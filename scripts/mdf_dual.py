#/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

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
   
    def __init__(self, S,
                 formation_energies=None, reaction_energies=None, fluxes=None):
        """Create a pathway object.
       
        Args:
            S: Stoichiometric matrix of the pathway.
                Reactions are on the rows, compounds on the columns.
            formation_energies: the Gibbs formation energy for the compounds
                in standard conditions, corrected for pH, ionic strength, etc.
                Should be a column vector in numpy.array format.
            reaction_energies: the change in Gibbs energy for the reactions
                in standard conditions, corrected for pH, ionic strength, etc.
                Should be a column vector in numpy.array format.
            fluxes: the list of relative fluxes through each of the reactions.
                By default, all fluxes are 1.
        """
        if formation_energies is None and reaction_energies is None:
            raise ValueError("In order to use 'Pathway' xeither formation xor "
                             "reaction energies must be provided.")
        if formation_energies is not None and reaction_energies is not None:
            raise ValueError("In order to use 'Pathway' xeither formation xor "
                             "reaction energies must be provided.")
       
        self.S = S
        self.dG0_f_prime = formation_energies
        self.dG0_r_prime = reaction_energies
        self.Nc, self.Nr = S.shape
       
        # make sure dG0_f' and dG0_r' are both 2D arrays of the right size
        if self.dG0_f_prime is not None:
            assert self.dG0_f_prime.shape[1] == self.Nc
            self.dG0_r_prime = self.CalculateReactionEnergies(self.dG0_f_prime)
        else:
            assert self.dG0_r_prime.shape[1] == self.Nr

        if fluxes is None:
            self.fluxes = np.matrix(np.ones((1, self.Nr)))
        elif type(fluxes) == types.ListType:
            self.fluxes = np.matrix(fluxes)
        else:
            self.fluxes = fluxes

        assert self.fluxes.shape[1] == self.Nr
           
        self.bounds = None
        self.c_range = self.DEFAULT_C_RANGE

    def CalculateReactionEnergies(self, dG_f):
        if np.isnan(dG_f).any():
            # if there are NaN values in dG_f, multiplying the matrices will not
            # work, since NumPy will not convert 0*NaN into 0 in the sum. Therefore,
            # the multiplication must be done explicitly and using only the nonzero
            # stoichiometric coefficients and their corresponding dG_f.
            dG_r = np.matrix(np.zeros((1, self.Nr)))
            for r in xrange(self.Nr):
                reactants = list(self.S[:, r].nonzero()[0].flat)
                dG_r[0, r] = dG_f[0, reactants] * self.S[reactants, r]
            return dG_r
        else:
            return dG_f * self.S

    def CalculateReactionEnergiesUsingConcentrations(self, concentrations):
        log_conc = np.log(concentrations)
        if np.isnan(self.dG0_r_prime).any(): # see CalculateReactionEnergies
            dG_r_prime = self.dG0_r_prime.copy()
            for r in xrange(self.Nr):
                reactants = list(self.S[:, r].nonzero()[0].flat)
                dG_r_prime[0, r] += default_RT * log_conc[reactants, 0].T * self.S[reactants, r]
            return dG_r_prime
        else:
            return self.dG0_r_prime + default_RT * log_conc.T * self.S

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
       
        if self.bounds:
            for i, bound in enumerate(self.bounds):
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
                subject to   Ax >= b
                             x >= 0
        """
        I_dir = np.matrix(np.diag([np.sign(x) for x in self.fluxes.flat]))
       
        A = np.matrix(np.vstack([np.hstack([I_dir * self.S.T, np.ones((self.Nr, 1)) ]),
                                 np.hstack([np.eye(self.Nc),  np.zeros((self.Nc, 1))]),
                                 np.hstack([-np.eye(self.Nc), np.zeros((self.Nc, 1))])]))
        b = np.matrix(np.vstack([-I_dir * self.dG0_r_prime.T / default_RT,
                                 ln_conc_ub,
                                 -ln_conc_lb]))
        c = np.matrix(np.vstack([np.zeros((self.Nc, 1)),
                                 np.ones((1, 1))]))
       
        return A, b, c
   
    def _GetTotalEnergyProblem(self, min_driving_force=0, objective=pulp.LpMinimize):
        
        # Define and apply the constraints on the concentrations
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds()

        # Create the driving force variable and add the relevant constraints
        A, b, _c = self._MakeDrivingForceConstraints(ln_conc_lb, ln_conc_ub)
       
        lp = pulp.LpProblem("MDF", objective)
        
        # ln-concentration variables
        l = pulp.LpVariable.dicts("l", ["%d" % i for i in xrange(self.Nc)])
        x = [l["%d" % i] for i in xrange(self.Nc)] + [min_driving_force]
        
        total_g = pulp.LpVariable("g_tot")
        
        for j in xrange(A.shape[0]):
            row = [A[j, i] * x[i] for i in xrange(A.shape[1])]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j
        
        total_g0 = float(self.dG0_r_prime * self.fluxes.T)
        total_reaction = self.S * self.fluxes.T
        row = [total_reaction[i, 0] * x[i] for i in xrange(self.Nc)]
        lp += (total_g == total_g0 + pulp.lpSum(row)), "Total G"

        lp.setObjective(total_g)
        
        #lp.writeLP("../res/total_g.lp")
        
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
        # Define and apply the constraints on the concentrations
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds()

        # Create the driving force variable and add the relevant constraints
        A, b, c = self._MakeDrivingForceConstraints(ln_conc_lb, ln_conc_ub)
       
        lp = pulp.LpProblem("MDF", pulp.LpMaximize)
        
        # ln-concentration variables
        l = pulp.LpVariable.dicts("l", ["%d" % i for i in xrange(self.Nc)])
        B = pulp.LpVariable("B")
        x = [l["%d" % i] for i in xrange(self.Nc)] + [B]
        
        for j in xrange(A.shape[0]):
            row = [A[j, i] * x[i] for i in xrange(A.shape[1])]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j
        
        objective = pulp.lpSum([c[i] * x[i] for i in xrange(A.shape[1])])
        lp.setObjective(objective)
        
        #lp.writeLP("../res/mdf_primal.lp")
        
        return lp, x, ln_conc_lb

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

        z = pulp.LpVariable.dicts("z", 
                                  ["%d" % i for i in xrange(self.Nc)],
                                  lowBound=0)

        u = pulp.LpVariable.dicts("u", 
                                  ["%d" % i for i in xrange(self.Nc)],
                                  lowBound=0)
        
        y = [w["%d" % i] for i in xrange(self.Nr)] + \
            [z["%d" % i] for i in xrange(self.Nc)] + \
            [u["%d" % i] for i in xrange(self.Nc)]
        
        for i in xrange(A.shape[1]):
            row = [A[j, i] * y[j] for j in xrange(A.shape[0])]
            lp += (pulp.lpSum(row) == c[i, 0]), "dual_%02d" % i

        objective = pulp.lpSum([b[i] * y[i] for i in xrange(A.shape[0])])
        lp.setObjective(objective)
        
        #lp.writeLP("../res/mdf_dual.lp")
        
        return lp, w, z, u
    
    def FindMDF(self):
        """Find the MDF (Optimized Bottleneck Energetics).
       
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
       
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal mdf).
        """
        lp_primal, x, _ = self._MakeMDFProblem()
        lp_primal.solve(pulp.CPLEX(msg=0))
        if lp_primal.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF primal")
            
        mdf = pulp.value(x[-1])
        conc = np.matrix([np.exp(pulp.value(x[j])) for j in xrange(self.Nc)]).T

        lp_dual, w, z, u = self._MakeMDFProblemDual()
        lp_dual.solve(pulp.CPLEX(msg=0))
        if lp_dual.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF dual")
        reaction_prices = np.matrix([pulp.value(w["%d" % i]) for i in xrange(self.Nr)]).T
        compound_prices = np.matrix([pulp.value(z["%d" % j]) for j in xrange(self.Nc)]).T - \
                          np.matrix([pulp.value(u["%d" % j]) for j in xrange(self.Nc)]).T
        
        # find the maximum and minimum total Gibbs energy of the pathway,
        # under the constraint that the driving force of each reaction is >= MDF
        lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-6, pulp.LpMinimize)
        lp_total.solve(pulp.CPLEX(msg=0))
        if lp_total.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve total delta-G problem")
        min_tot_dg = pulp.value(total_dg)

        lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-6, pulp.LpMaximize)
        lp_total.solve(pulp.CPLEX(msg=0))
        if lp_total.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve total delta-G problem")
        max_tot_dg = pulp.value(total_dg)
        
        params = {'MDF': mdf * default_RT,
                  'concentrations' : conc,
                  'reaction prices' : reaction_prices,
                  'compound prices' : compound_prices,
                  'maximum total dG' : max_tot_dg * default_RT,
                  'minimum total dG' : min_tot_dg * default_RT}
        return mdf * default_RT, params

class KeggPathway(Pathway):
   
    def __init__(self, S, rids, fluxes, cids, formation_energies=None,
                 reaction_energies=None, cid2bounds=None, c_range=None):
        Pathway.__init__(self, S, formation_energies=formation_energies,
                         reaction_energies=reaction_energies, fluxes=fluxes)
        assert len(cids) == self.Nc
        assert len(rids) == self.Nr
       
        self.rids = rids
        self.cids = cids
        if cid2bounds:
            self.bounds = [cid2bounds.get(cid, (None, None)) for cid in self.cids]
        else:
            self.bounds = None
        self.cid2bounds = cid2bounds
        self.c_range = c_range

    def GetConcentrationBounds(self, cid):
        lb, ub = None, None
        if cid in self.cid2bounds:
            lb, ub = self.cid2bounds[cid]
        lb = lb or self.c_range[0]
        ub = ub or self.c_range[1]
        return lb, ub

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
        Nr = dGs.shape[1]
        dGs_adjusted = np.multiply(dGs, fluxes)
        cum_dG = np.cumsum([0] + [dGs_adjusted[0, r] for r in xrange(Nr)])
        plt.plot(np.arange(0.5, Nr + 1), cum_dG, style,
                 figure=figure, label=label)
   
    def PlotProfile(self, concentrations, figure=None):
        if figure is None:
            figure = plt.figure(figsize=(6,6), dpi=100)
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
            self.dG0_r_prime[0, nonzero_reactions], nonzero_fluxes, 'm--', r"$\Delta_r G'^\circ$")
       
        phys_concentrations = self.GetPhysiologicalConcentrations(self.bounds)
        dGm_r_prime = self.CalculateReactionEnergiesUsingConcentrations(phys_concentrations)
        KeggPathway._AddProfileToFigure(figure,
            dGm_r_prime[0, nonzero_reactions], nonzero_fluxes, 'g--', r"$\Delta_r G'^m$")

        dG_r_prime = self.CalculateReactionEnergiesUsingConcentrations(concentrations)
        KeggPathway._AddProfileToFigure(figure,
            dG_r_prime[0, nonzero_reactions], nonzero_fluxes, 'b-', r"$\Delta_r G'$")

        plt.legend(loc='lower left')
        return figure
       
    def PlotConcentrations(self, concentrations, figure=None):
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

    def WriteResultsToHtmlTables(self, html_writer, concentrations,
                                 reaction_shadow_prices, compound_shadow_prices):
        
        self.WriteProfileToHtmlTable(html_writer, concentrations, reaction_shadow_prices)
        self.WriteConcentrationsToHtmlTable(html_writer, concentrations, compound_shadow_prices)

    def WriteConcentrationsToHtmlTable(self, html_writer, concentrations=None,
                                       compound_shadow_prices=None):
        phys_concentrations = self.GetPhysiologicalConcentrations(self.bounds)
        if concentrations is None:
            concentrations = phys_concentrations
        if compound_shadow_prices is None:
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
   
    def WriteProfileToHtmlTable(self, html_writer, concentrations=None,
                                reaction_shadow_prices=None):

        phys_concentrations = self.GetPhysiologicalConcentrations(self.bounds)
        if concentrations is None:
            concentrations = phys_concentrations
        if reaction_shadow_prices is None:
            reaction_shadow_prices = np.zeros((len(self.rids), 1))

        dG_r_prime_c = self.CalculateReactionEnergiesUsingConcentrations(phys_concentrations)
        dG_r_prime_c_adj = np.multiply(dG_r_prime_c, np.sign(self.fluxes)) # adjust dG to flux directions
        dG_r_prime = self.CalculateReactionEnergiesUsingConcentrations(concentrations)
        dG_r_prime_adj = np.multiply(dG_r_prime, np.sign(self.fluxes)) # adjust dG to flux directions
        headers=["reaction", 'formula', 'flux',
                 "&Delta;<sub>r</sub>G'<sup>c</sup> [kJ/mol] (%g M)" % self.DEFAULT_PHYSIOLOGICAL_CONC,
                 "&Delta;<sub>r</sub>G' [kJ/mol]", "shadow price"]

        dict_list = []
        for r, rid in enumerate(self.rids):
            d = {}
            d['reaction'] = rid
            d['flux'] = "%g" % abs(self.fluxes[0, r])
            d['formula'] = self.GetReactionString(r, show_cids=False)
            d[headers[3]] = dG_r_prime_c_adj[0, r]
            d[headers[4]] = dG_r_prime_adj[0, r]
            d[headers[5]] = '%.3g' % reaction_shadow_prices[r, 0]
                
            dict_list.append(d)

        d = {'reaction':'Total',
             'flux':'1',
             'formula':self.GetTotalReactionString(show_cids=False),
             headers[3]: float(self.dG0_r_prime * self.fluxes.T),
             headers[4]: float(dG_r_prime * self.fluxes.T),
             headers[5]: '%.3g' % np.sum(reaction_shadow_prices[:, 0])}
        dict_list.append(d)
        
        html_writer.write_table(dict_list, headers=headers, decimal=1)
       
if __name__ == '__main__':
    dGs = np.matrix([[0, -50, -30, -100]])
    n = dGs.shape[1]
    S = np.matrix(np.vstack([-np.eye(n-1), np.zeros((1, n-1))]) + np.vstack([np.zeros((1, n-1)), np.eye(n-1)]))
    print 'S = %s' % str(S)
    fluxes = np.matrix(np.ones((1, n-1)))
    rids = range(n-1)
    cids = range(n)
    
    print 'G0 = %s' % str(dGs * S)
    #keggpath = KeggPathway(S, rids, fluxes, cids, dGs, c_range=(1e-6, 1e-3))
    #dGf, concentrations, protein_cost = keggpath.FindKineticOptimum()
    #print 'protein cost: %g protein units' % protein_cost
    #print 'concentrations: ', concentrations
    #print 'sum(concs): ', sum(concentrations), 'M'
   
    keggpath = KeggPathway(S, rids, fluxes, cids, dGs, c_range=(1e-6, 1e-3))
    mdf, params = keggpath.FindMDF()
    print 'MDF: %g' % mdf
    print 'reaction shadow prices: ' + ', '.join(['%g' % i for i in params['reaction prices'].flat])
    print 'compound shadow prices: ' + ', '.join(['%g' % i for i in params['compound prices'].flat])
    print 'concentrations: ' + ', '.join(['%.2e' % i for i in params['concentrations'].flat])
    print 'minimal total dG: %.1f' % params['minimum total dG']
    print 'maximal total dG: %.1f' % params['maximum total dG']
    
    G_opt = (dGs + default_RT * np.log(params['concentrations'].T))*S
    print 'G_opt = %s' % str(G_opt)
