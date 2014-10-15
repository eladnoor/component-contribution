#/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import logging

from matplotlib.font_manager import FontProperties
import types
import pulp

from python.thermodynamic_constants import default_RT
from python.kegg_reaction import KeggReaction

class Pathway(object):
    """Container for doing pathway-level thermodynamic analysis."""
   
    DEFAULT_FORMATION_LB = -1e6
    DEFAULT_FORMATION_UB = 1e6
    DEFAULT_REACTION_LB = -1e3
    DEFAULT_REACTION_UB = 0.0
    DEFAULT_C_RANGE = (1e-6, 0.1)
    DEFAULT_PHYSIOLOGICAL_CONC = 1e-3
   
    def __init__(self, S, dG0_r_prime, dG0_r_std=None, fluxes=None):
        """Create a pathway object.
       
        Args:
            S: Stoichiometric matrix of the pathway.
                Reactions are on the rows, compounds on the columns.
            dG0_r_prime: the change in Gibbs energy for the reactions
                in standard conditions, corrected for pH, ionic strength, etc.
                Should be a column vector in numpy.matrix format.
            dG0_r_std: (optional) the square root of the covariance matrix
                corresponding to the uncertainty in the dG0_r values.
            fluxes: the list of relative fluxes through each of the reactions.
                By default, all fluxes are 1.
        """
        self.S = S
        self.Nc, self.Nr = S.shape

        self.dG0_r_prime = dG0_r_prime
        if dG0_r_std is None:
            self.dG0_r_std = np.matrix(np.zeros((self.Nr, self.Nr)))
        else:
            self.dG0_r_std = dG0_r_std
       
        # make sure dG0_r' is the right size
        assert self.dG0_r_prime.shape[0] == self.Nr
        assert self.dG0_r_std.shape[0] == self.Nr
        assert self.dG0_r_std.shape[1] == self.Nr

        if fluxes is None:
            self.fluxes = np.matrix(np.ones((1, self.Nr)))
        elif type(fluxes) == types.ListType:
            self.fluxes = np.matrix(fluxes)
        else:
            self.fluxes = fluxes

        assert self.fluxes.shape[1] == self.Nr
        
        self.I_dir = np.matrix(np.diag(map(np.sign, self.fluxes.flat)))
        self.c_bounds = None
        self.r_bounds = None
        self.c_range = self.DEFAULT_C_RANGE

    def CalculateReactionEnergiesUsingConcentrations(self, concentrations):
        log_conc = np.log(concentrations)
        if np.isnan(self.dG0_r_prime).any():
            dG_r_prime = self.dG0_r_prime.copy()
            for r in xrange(self.Nr):
                reactants = list(self.S[:, r].nonzero()[0].flat)
                dG_r_prime[0, r] += default_RT * log_conc[reactants, 0].T * self.S[reactants, r]
            return dG_r_prime
        else:
            return self.dG0_r_prime + default_RT * self.S.T * log_conc

    def GetPhysiologicalConcentrations(self, bounds=None):
        conc = np.matrix(np.ones((self.Nc, 1))) * self.DEFAULT_PHYSIOLOGICAL_CONC
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                if lb is not None and ub is not None:
                    if not (lb < conc[i, 0] < ub):
                        conc[i, 0] = np.sqrt(lb * ub)
       
        return conc

    def _MakeLnConcentratonBounds(self):
        """Make bounds on logarithmic concentrations."""
        c_lower, c_upper = self.c_range or self.DEFAULT_C_RANGE
        ln_conc_lb = np.matrix(np.ones((self.Nc, 1)) * np.log(c_lower))
        ln_conc_ub = np.matrix(np.ones((self.Nc, 1)) * np.log(c_upper))
       
        if self.c_bounds:
            for i, bound in enumerate(self.c_bounds):
                lb, ub = bound
                log_lb = np.log(lb or c_lower)
                log_ub = np.log(ub or c_upper)
                if log_lb > log_ub:
                    raise Exception("Lower bound is greater than upper bound: "
                                    "%d > %d" % (log_lb, log_ub))
                elif abs(log_lb - log_ub) < 1e-2:
                    log_lb = log_ub - 1e-2
                   
                ln_conc_lb[i, 0] = log_lb
                ln_conc_ub[i, 0] = log_ub

        return ln_conc_lb, ln_conc_ub

    def _MakeDrivingForceConstraints(self, ln_conc_lb, ln_conc_ub):
        """
            Generates the A matrix and b & c vectors that can be used in a 
            standard form linear problem:
                max          c'x
                subject to   Ax <= b
                            
            x is the vector of (y | log-conc | B)
            where y dG'0 are the reaction Gibbs energy variables, log-conc
            are the natural log of the concentrations of metabolites, and
            B is the max-min driving force variable which is being maximized
            by the LP
        """
        
        
       
        A = np.matrix(np.vstack([np.hstack([self.dG0_r_std / default_RT, self.I_dir * self.S.T, np.ones((self.Nr, 1)) ]),
                                 np.hstack([ np.eye(self.Nr),              np.zeros((self.Nr, self.Nc+1))]),
                                 np.hstack([-np.eye(self.Nr),              np.zeros((self.Nr, self.Nc+1))]),
                                 np.hstack([np.zeros((self.Nc, self.Nr)),  np.eye(self.Nc), np.zeros((self.Nc, 1))]),
                                 np.hstack([np.zeros((self.Nc, self.Nr)), -np.eye(self.Nc), np.zeros((self.Nc, 1))])]))
        b = np.matrix(np.vstack([-self.I_dir * self.dG0_r_prime / default_RT,
                                 np.ones((self.Nr, 1)),
                                 np.ones((self.Nr, 1)),
                                  ln_conc_ub,
                                 -ln_conc_lb]))
        c = np.matrix(np.vstack([np.zeros((self.Nr + self.Nc, 1)),
                                 np.ones((1, 1))]))

        # change the constaints such that reaction that have an explicit
        # r_bound will not be constrained by B, but will be constained by
        # their specific bounds. Note that we need to divide the bound
        # by R*T since the variables in the LP are not in kJ/mol but in units
        # of R*T.
        if self.r_bounds:
            for i, r_ub in enumerate(self.r_bounds):
                if r_ub is not None:
                    A[i, -1] = 0
                    b[i] += r_ub/default_RT + 1e-5 
        
        return A, b, c
   
    def _MakeBasicProblem(self, mdf):
        # Define and apply the constraints on the concentrations
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds()

        # Create the driving force variable and add the relevant constraints
        A, b, c = self._MakeDrivingForceConstraints(ln_conc_lb, ln_conc_ub)
       
        # the dG'0 covariance eigenvariables        
        y = pulp.LpVariable.dicts("y", ["%d" % i for i in xrange(self.Nr)])
        y = [y["%d" % i] for i in xrange(self.Nr)]

        # ln-concentration variables
        l = pulp.LpVariable.dicts("l", ["%d" % i for i in xrange(self.Nc)])
        l = [l["%d" % i] for i in xrange(self.Nc)]

        x = y + l + [mdf]
        return A, b, c, x, y, l
   
    def _GetTotalEnergyProblem(self,
                               min_driving_force=0.0,
                               objective=pulp.LpMinimize):
        
        total_g = pulp.LpVariable("g_tot")
        A, b, _c, x, y, l = self._MakeBasicProblem(min_driving_force)

        lp = pulp.LpProblem("MDF", objective)
        
        for j in xrange(3*self.Nr + 2 * self.Nc):
            row = [A[j, i] * x[i] for i in xrange(self.Nr + self.Nc + 1)]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j
        
        total_g0 = float(self.fluxes * self.dG0_r_prime)
        total_reaction = self.S * self.fluxes.T
        row = [total_reaction[i, 0] * x[i] for i in xrange(self.Nc)]
        lp += (total_g == total_g0 + pulp.lpSum(row)), "Total G"

        lp.setObjective(total_g)
        
        #lp.writeLP("res/total_g.lp")
        
        return lp, total_g
           
    def _MakeMDFProblem(self):
        """Create a CVXOPT problem for finding the Maximal Thermodynamic
        Driving Force (MDF).
       
        Does not set the objective function... leaves that to the caller.
       
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
       
        Returns:
            A tuple (dgf_var, motive_force_var, problem_object).
        """
        B = pulp.LpVariable("mdf")
        A, b, c, x, y, l = self._MakeBasicProblem(B)

        lp = pulp.LpProblem("MDF", pulp.LpMaximize)
        
        for j in xrange(3*self.Nr + 2 * self.Nc):
            row = [A[j, i] * x[i] for i in xrange(self.Nr + self.Nc + 1)]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j
        
        objective = pulp.lpSum([c[i] * x[i] for i in xrange(A.shape[1])])
        lp.setObjective(objective)
        
        #lp.writeLP("res/mdf_primal.lp")
        
        return lp, y, l, B

    def _MakeMDFProblemDual(self):
        """Create a CVXOPT problem for finding the Maximal Thermodynamic
        Driving Force (MDF).
       
        Does not set the objective function... leaves that to the caller.
       
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
       
        Returns:
            A tuple (dgf_var, motive_force_var, problem_object).
        """
        # Define and apply the constraints on the concentrations
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds()

        # Create the driving force variable and add the relevant constraints
        A, b, c = self._MakeDrivingForceConstraints(ln_conc_lb, ln_conc_ub)
       
        lp = pulp.LpProblem("MDF", pulp.LpMinimize)
        
        w = pulp.LpVariable.dicts("w", 
                                  ["%d" % i for i in xrange(self.Nr)],
                                  lowBound=0)
        w = [w["%d" % i] for i in xrange(self.Nr)]

        g = pulp.LpVariable.dicts("g", 
                                  ["%d" % i for i in xrange(2*self.Nr)],
                                  lowBound=0)
        g = [g["%d" % i] for i in xrange(2*self.Nr)]

        z = pulp.LpVariable.dicts("z", 
                                  ["%d" % i for i in xrange(self.Nc)],
                                  lowBound=0)
        z = [z["%d" % i] for i in xrange(self.Nc)]

        u = pulp.LpVariable.dicts("u", 
                                  ["%d" % i for i in xrange(self.Nc)],
                                  lowBound=0)
        u = [u["%d" % i] for i in xrange(self.Nc)]
        
        x = w + g + z + u
        
        for i in xrange(A.shape[1]):
            row = [A[j, i] * x[j] for j in xrange(A.shape[0])]
            lp += (pulp.lpSum(row) == c[i, 0]), "dual_%02d" % i

        objective = pulp.lpSum([b[i] * x[i] for i in xrange(A.shape[0])])
        lp.setObjective(objective)
        
        #lp.writeLP("res/mdf_dual.lp")
        
        return lp, w, g, z, u
    
    def FindMDF(self, calculate_totals=True):
        """Find the MDF (Optimized Bottleneck Energetics).
       
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
       
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal mdf).
        """
        lp_primal, y, l, B = self._MakeMDFProblem()
        lp_primal.solve(pulp.CPLEX(msg=0))
        if lp_primal.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF primal")
        y = map(pulp.value, y)
        l = map(pulp.value, l)
        mdf = pulp.value(B)
        conc = np.matrix(map(np.exp, l)).T
        dG0_r_prime = self.dG0_r_prime + np.dot(self.dG0_r_std, y).T

        lp_dual, w, g, z, u = self._MakeMDFProblemDual()
        
        lp_dual.solve(pulp.CPLEX(msg=0))
        if lp_dual.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF dual")
        w = map(pulp.value, w)
        z = map(pulp.value, z)
        u = map(pulp.value, u)
        reaction_prices = np.matrix(w).T
        compound_prices = np.matrix(z).T - np.matrix(u).T

        params = {'MDF': mdf * default_RT,
                  'reaction energies': dG0_r_prime,
                  'concentrations' : conc,
                  'reaction prices' : reaction_prices,
                  'compound prices' : compound_prices}
                  
        if calculate_totals:
            # find the maximum and minimum total Gibbs energy of the pathway,
            # under the constraint that the driving force of each reaction is >= MDF
            lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-6, pulp.LpMinimize)
            lp_total.solve(pulp.CPLEX(msg=0))
            if lp_total.status != pulp.LpStatusOptimal:
                #raise pulp.solvers.PulpSolverError("cannot solve minimal total delta-G problem")
                logging.warning("cannot solve minimal total delta-G problem")
                min_tot_dg = np.nan
            else:
                min_tot_dg = pulp.value(total_dg)
    
            params['minimum total dG'] = min_tot_dg * default_RT
    
            lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-6, pulp.LpMaximize)
            lp_total.solve(pulp.CPLEX(msg=0))
            if lp_total.status != pulp.LpStatusOptimal:
                #raise pulp.solvers.PulpSolverError("cannot solve maximal total delta-G problem")
                logging.warning("cannot solve maximal total delta-G problem")
                max_tot_dg = np.nan
            else:
                max_tot_dg = pulp.value(total_dg)

            params['maximum total dG'] = max_tot_dg * default_RT
            
        dG_r_prime = self.CalculateReactionEnergiesUsingConcentrations(conc)
        params['gibbs energies raw'] = dG_r_prime + np.dot(self.dG0_r_std, y).T

        # adjust dG to flux directions
        dG_r_prime_adj = self.I_dir * params['gibbs energies raw']
        params['gibbs energies'] = dG_r_prime_adj
        
        return mdf * default_RT, params

class KeggPathway(Pathway):
   
    def __init__(self, S, rids, fluxes, cids,
                 dG0_r_prime, dG0_r_std=None,
                 rid2bounds=None,
                 cid2bounds=None,
                 c_range=None):
                     
        """
            S           - the stoichiometric matrix
            rid         - a list of names of the reactions in S
            fluxes      - a vector the relative fluxes in each reaction
            cids        - a list of names of the compounds in S
            formation_energies - the standard Gibbs energy of formation of the 
                                 compounds
            rid2bounds - a dictionary mapping rid to an upper bound on its dG'
                         if the value is None then the upper bound
                         is the B variable (corresponding to the MDF)
            reaction_energies - the standard Gibbs energies of the reactions
            cid2bounds - a dictionary mapping cid to a pair of lower/upper
                         bound on its concentration. if the value is (None, None)
                         the default bounds are used (i.e. c_range)
            c_range    - the default lower/upper bounds on the compound conc.
        """
        Pathway.__init__(self, S, dG0_r_prime=dG0_r_prime, dG0_r_std=dG0_r_std, fluxes=fluxes)
        assert len(cids) == self.Nc
        assert len(rids) == self.Nr
       
        self.rids = rids
        self.cids = cids
        if cid2bounds:
            self.c_bounds = [cid2bounds.get(cid, (None, None)) for cid in self.cids]
        else:
            self.c_bounds = None
        self.cid2bounds = cid2bounds
        
        if rid2bounds:
            self.r_bounds = [rid2bounds.get(rid, None) for rid in self.rids]
        else:
            self.r_bounds = None
        
        self.rid2bounds = rid2bounds
        self.c_range = c_range

    def GetConcentrationBounds(self, cid):
        lb, ub = None, None
        if cid in self.cid2bounds:
            lb, ub = self.cid2bounds[cid]
        lb = lb or self.c_range[0]
        ub = ub or self.c_range[1]
        return lb, ub

    def GetMillimolarConcentrations(self):
        conc = np.matrix(np.ones((self.Nc, 1))) * self.DEFAULT_PHYSIOLOGICAL_CONC
        try:
            i_h2o = self.cids.index('C00001')
            conc[i_h2o, 0] = 1
        except ValueError:
            pass
        return conc

    def GetReactionString(self, r, show_cids=False):
        rid = self.rids[r]
        sparse = dict([(self.cids[c], self.S[c, r])
                       for c in self.S[:, r].nonzero()[0].flat])
        if self.fluxes[0, r] >= 0:
            direction = '=>'
        else:
            direction = '<='
        reaction = KeggReaction(sparse, arrow=direction, rid=rid)
        return str(reaction)

    def GetTotalReactionString(self, show_cids=False):
        total_S = self.S * self.fluxes.T
        sparse = dict([(self.cids[c], total_S[c, 0])
                       for c in total_S.nonzero()[0].flat])
        reaction = KeggReaction(sparse, arrow="=>", rid="Total")
        return str(reaction)

    @staticmethod
    def _EnergyToString(dG):
        if np.isnan(dG):
            return "N/A"
        else:
            return "%.1f" % dG

    @staticmethod
    def _AddProfileToFigure(figure, dGs, fluxes, style, label):
        Nr = dGs.shape[0]
        dGs_adjusted = np.diag(fluxes.flat) * dGs
        cum_dG = np.cumsum([0] + list(dGs_adjusted.flat))
        plt.plot(np.arange(0.5, Nr + 1), cum_dG, style,
                 figure=figure, label=label)
   
    def PlotProfile(self, params, figure=None):
        phys_concentrations = self.GetMillimolarConcentrations()

        if figure is None:
            figure = plt.figure(figsize=(8,8), dpi=100)
        plt.title(r'Thermodynamic profile', figure=figure)
        plt.ylabel(r"cumulative $\Delta_r G'*$ [kJ/mol]", figure=figure)
        plt.xlabel(r'Reaction', figure=figure)
       
        nonzero_reactions = list(np.nonzero(self.fluxes)[1].flat)
        nonzero_fluxes = self.fluxes[0, nonzero_reactions]
       
        plt.xticks(np.arange(0, len(nonzero_reactions)) + 1,
                   [self.rids[i] for i in nonzero_reactions],
                   fontproperties=FontProperties(size=10), rotation=30)

        if np.isnan(self.dG0_r_prime).any():
            return figure
       
        KeggPathway._AddProfileToFigure(figure,
            self.dG0_r_prime[nonzero_reactions, 0], nonzero_fluxes, 'm--', r"$\Delta_r G'^\circ$")
       
        dGm_r_prime = self.CalculateReactionEnergiesUsingConcentrations(phys_concentrations)
        KeggPathway._AddProfileToFigure(figure,
            dGm_r_prime[nonzero_reactions, 0], nonzero_fluxes, 'g--', r"$\Delta_r G'^m$")

        dG_r_prime = params['gibbs energies raw']
        KeggPathway._AddProfileToFigure(figure,
            dG_r_prime[nonzero_reactions, 0], nonzero_fluxes, 'b-', r"$\Delta_r G'$")

        plt.legend(loc='lower left')
        return figure
       
    def PlotConcentrations(self, params, figure=None):
        concentrations = params['concentrations']

        if figure is None:
            figure = plt.figure()
        plt.xscale('log', figure=figure)
        plt.ylabel('Compound KEGG ID', figure=figure)
        plt.xlabel('Concentration [M]', figure=figure)
        plt.yticks(range(self.Nc, 0, -1), self.cids,
                   fontproperties=FontProperties(size=8))
        plt.plot(concentrations.T, range(self.Nc, 0, -1), '*b', figure=figure)

        x_min = concentrations.min() / 10
        x_max = concentrations.max() * 10
        y_min = 0
        y_max = self.Nc + 1
       
        for c, cid in enumerate(self.cids):
            plt.text(concentrations[0, c] * 1.1, self.Nc - c, cid, \
                       figure=figure, fontsize=6, rotation=0)
            b_low, b_up = self.GetConcentrationBounds(cid)
            plt.plot([b_low, b_up], [self.Nc - c, self.Nc - c], '-k', linewidth=0.4)

        if self.c_range is not None:
            plt.axvspan(self.c_range[0], self.c_range[1],
                          facecolor='r', alpha=0.3, figure=figure)
        plt.axis([x_min, x_max, y_min, y_max], figure=figure)
        return figure

    def WriteResultsToHtmlTables(self, html_writer, params):
        self.WriteProfileToHtmlTable(html_writer, params)
        self.WriteConcentrationsToHtmlTable(html_writer, params)

    def WriteProfileToHtmlTable(self, html_writer, params=None):
        phys_concentrations = self.GetMillimolarConcentrations()
        if params is not None:
            reaction_shadow_prices = params['reaction prices']
        else:
            reaction_shadow_prices = np.zeros((len(self.rids), 1))

        dG_r_prime_c     = self.CalculateReactionEnergiesUsingConcentrations(phys_concentrations)
        dG_r_prime_c_adj = self.I_dir * dG_r_prime_c # adjust dG to flux directions
        dG_r_prime       = params['gibbs energies raw']
        dG_r_prime_adj   = params['gibbs energies']
        headers=["reaction", 'formula', 'flux',
                 "&Delta;<sub>r</sub>G'<sup>m</sup> [kJ/mol] (%g M)" % self.DEFAULT_PHYSIOLOGICAL_CONC,
                 "&Delta;<sub>r</sub>G' [kJ/mol]", "shadow price"]

        dict_list = []
        for r, rid in enumerate(self.rids):
            d = {}
            d['reaction'] = rid
            d['flux'] = "%g" % abs(self.fluxes[0, r])
            d['formula'] = self.GetReactionString(r, show_cids=False)
            d[headers[3]] = dG_r_prime_c_adj[r, 0]
            d[headers[4]] = dG_r_prime_adj[r, 0]
            d[headers[5]] = '%.3g' % reaction_shadow_prices[r, 0]
                
            dict_list.append(d)

        d = {'reaction':'Total',
             'flux':'1',
             'formula':self.GetTotalReactionString(show_cids=False),
             headers[3]: float(self.fluxes * self.dG0_r_prime),
             headers[4]: float(self.fluxes * dG_r_prime),
             headers[5]: '%.3g' % np.sum(reaction_shadow_prices[:, 0])}
        dict_list.append(d)
        
        html_writer.write_table(dict_list, headers=headers, decimal=1)

    def WriteConcentrationsToHtmlTable(self, html_writer, params=None):
        if params is not None:
            concentrations = params['concentrations']
            compound_shadow_prices = params['compound prices']
        else:
            concentrations = self.GetMillimolarConcentrations()
            compound_shadow_prices = np.zeros((len(self.cids), 1))

        dict_list = []
        for c, cid in enumerate(self.cids):
            d = {}
            d['KEGG CID'] = cid
            lb, ub = self.GetConcentrationBounds(cid)
            d['Concentration LB [M]'] = '%.2e' % lb
            d['Concentration [M]'] = '%.2e' % concentrations[c, 0]
            d['Concentration UB [M]'] = '%.2e' % ub
            d['shadow price'] = '%.3g' % compound_shadow_prices[c, 0]
            dict_list.append(d)
        headers = ['KEGG CID', 'Concentration LB [M]',
                   'Concentration [M]', 'Concentration UB [M]', 'shadow price']
       
        html_writer.write_table(dict_list, headers=headers)
   
       
if __name__ == '__main__':
    dG0_r_prime = np.matrix([[-10, -10, -10]]).T
    dG0_r_std = 3.0 * np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    #dG0_r_std = None
    
    Nr = dG0_r_prime.shape[0]
    Nc = Nr+1
    S = np.matrix(np.vstack([-np.eye(Nr), np.zeros((1, Nr))]) + np.vstack([np.zeros((1, Nr)), np.eye(Nr)]))
    fluxes = np.matrix(np.ones((1, Nr)))
    rids = range(Nr)
    cids = range(Nc)
    
    print 'dG0_r_prime = ' + ', '.join(['%.2f' % i for i in dG0_r_prime.flat])
    #keggpath = KeggPathway(S, rids, fluxes, cids, dGs, c_range=(1e-6, 1e-3))
    #dGf, concentrations, protein_cost = keggpath.FindKineticOptimum()
    #print 'protein cost: %g protein units' % protein_cost
    #print 'concentrations: ', concentrations
    #print 'sum(concs): ', sum(concentrations), 'M'
   
    keggpath = KeggPathway(S, rids, fluxes, cids, dG0_r_prime, dG0_r_std, c_range=(1e-6, 1e-3))
    mdf, params = keggpath.FindMDF()
    print 'MDF: %.2f' % mdf
    print 'reaction shadow prices: ' + ', '.join(['%g' % i for i in params['reaction prices'].flat])
    print 'compound shadow prices: ' + ', '.join(['%g' % i for i in params['compound prices'].flat])
    print 'concentrations: ' + ', '.join(['%.2e' % i for i in params['concentrations'].flat])
    print 'minimal total dG: %.2f' % params['minimum total dG']
    print 'maximal total dG: %.2f' % params['maximum total dG']
    
    print 'dG_r_prime(optimal) = ' + ', '.join(['%.2f' % i for i in params['gibbs energies'].flat])
