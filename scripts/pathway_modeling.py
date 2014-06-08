# -*- coding: utf-8 -*-
"""
Created on Wed May 28 13:22:51 2014

@author: eladn
"""
import numpy as np
import types
import cvxpy
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from python.thermodynamic_constants import default_RT, symbol_dr_G0_prime, \
                                           symbol_dr_G_prime, symbol_dr_Gc_prime
                                           
from python.kegg_reaction import KeggReaction

class DeltaGNormalization:
    TIMES_FLUX = 1 # motivation is having the most uniform entropy production
    DIVIDE_BY_FLUX = 2 # motivation is requiring more force in reaction that have more flux
    SIGN_FLUX = 3 # motivation is putting limits on the allowed backward/forward fluxes
    
    DEFAULT = SIGN_FLUX

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
        else:
            assert fluxes.shape[1] == self.Nr
            self.fluxes = fluxes
            
        self.normalization = DeltaGNormalization.DEFAULT

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
                dG_r_prime[0, r] += default_RT * log_conc[0, reactants] * self.S[reactants, r]
            return dG_r_prime
        else:
            return self.dG0_r_prime + default_RT * log_conc * self.S

    def GetPhysiologicalConcentrations(self, bounds=None):
        conc = np.matrix(np.ones((1, self.Nc))) * self.DEFAULT_PHYSIOLOGICAL_CONC
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                if lb is not None and ub is not None:
                    if not (lb < conc[0, i] < ub):
                        conc[0, i] = np.sqrt(lb * ub)
        
        return conc

    def _MakeLnConcentratonBounds(self, ln_conc, bounds=None, c_range=None):
        """Make bounds on logarithmic concentrations."""
        c_lower, c_upper = c_range or self.DEFAULT_C_RANGE
        ln_conc_lb = np.matrix(np.ones((1, self.Nc))) * np.log(c_lower)
        ln_conc_ub = np.matrix(np.ones((1, self.Nc))) * np.log(c_upper)
        
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                log_lb = np.log(lb or c_lower)
                log_ub = np.log(ub or c_upper)
                if log_lb > log_ub:
                    raise Exception("Lower bound is greater than upper bound: "
                                    "%d > %d" % (log_lb, log_ub))
                elif abs(log_lb - log_ub) < 1e-2:
                    log_lb = log_ub - 1e-2
                    
                ln_conc_lb[0, i] = log_lb
                ln_conc_ub[0, i] = log_ub
        
        return [cvxpy.geq(ln_conc, cvxpy.matrix(ln_conc_lb)) + \
                cvxpy.leq(ln_conc, cvxpy.matrix(ln_conc_ub))]

    def _MakeDrivingForceConstraints(self, ln_conc, driving_force_lb=0):
        """
            driving_force_lb can either be a cvxpy variable use later in the optimization
            or a scalar, which sets it as a constraint. By default the lower bound is 0.
        """
        constraints = []
        S = cvxpy.matrix(self.S)
        dg0r_primes = cvxpy.matrix(self.dG0_r_prime)
        for i in xrange(self.Nr):
            # if the dG0 is unknown, this reaction imposes no new constraints
            if np.isnan(self.dG0_r_prime[0, i]):
                continue
            
            curr_dgr = dg0r_primes[0, i] + default_RT * ln_conc * S[:, i]
            if self.fluxes[0, i] == 0:
                constraints += cvxpy.eq(curr_dgr, 0)
            else:
                if self.normalization == DeltaGNormalization.DIVIDE_BY_FLUX:
                    motive_force = -curr_dgr * (1.0 / self.fluxes[0, i])
                elif self.normalization == DeltaGNormalization.TIMES_FLUX:
                    motive_force = -curr_dgr * self.fluxes[0, i]
                elif self.normalization == DeltaGNormalization.SIGN_FLUX:
                    motive_force = -curr_dgr * np.sign(self.fluxes[0, i])
                else:
                    raise ValueError("bad value for normalization method: "
                                     + str(self.normalization))

                constraints += [cvxpy.geq(motive_force, driving_force_lb)]

        return constraints

    def _GetTotalReactionEnergy(self, c_range=(1e-6, 1e-2), bounds=None, min_driving_force=0):
        constraints = []
        
        # Define and apply the constraints on the concentrations
        ln_conc = cvxpy.variable(1, self.Nc, name='lnC')
        constraints += self._MakeLnConcentratonBounds(ln_conc, bounds=bounds,
                                                      c_range=c_range)
        
        # find the row vector describing the overall stoichiometry
        S = cvxpy.matrix(self.S)
        f = cvxpy.matrix(self.fluxes)
        g0 = cvxpy.matrix(self.dG0_r_prime)
        g = g0 + default_RT * ln_conc * S
        total_g = f * g.T

        constraints += self._MakeDrivingForceConstraints(ln_conc, min_driving_force)
        
        return ln_conc, constraints, total_g
            
    def _MakeMDFProblem(self, c_range=(1e-6, 1e-2), bounds=None):
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
        constraints = []
        
        # Define and apply the constraints on the concentrations
        ln_conc = cvxpy.variable(1, self.Nc, name='lnC')
        constraints += self._MakeLnConcentratonBounds(ln_conc, bounds=bounds,
                                                      c_range=c_range)

        # Create the driving force variable and add the relevant constraints
        driving_force_lb = cvxpy.variable(name='B')
        constraints += self._MakeDrivingForceConstraints(ln_conc, driving_force_lb)
        return ln_conc, driving_force_lb, constraints

    def _FindMDF(self, c_range=(1e-6, 1e-2), bounds=None,
                  normalization=DeltaGNormalization.DEFAULT):
        """Find the MDF (Maximal Thermodynamic Driving Force).
        
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal mdf).
        """
        ln_conc, motive_force_lb, constraints = self._MakeMDFProblem(
                                                c_range, bounds)
        program = cvxpy.program(cvxpy.maximize(motive_force_lb), constraints)
        program.solve(quiet=True)
        return ln_conc.value, program.objective.value

    def FindMDF_Regularized(self, c_range=(1e-6, 1e-2), bounds=None,
                            c_mid=1e-3,
                            min_mdf=None,
                            max_mdf=None):
        """Find the MDF (Max-min Driving Force).
        
        Uses l2 regularization to minimize the log difference of 
        concentrations from c_mid.
        
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            c_mid: the defined midpoint concentration.
            max_mdf: the maximum value for the motive force.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal mdf).
        """
        ln_conc, motive_force_lb, constraints = self._MakeMDFProblem(c_range, bounds)
        
        # Set the objective and solve.
        norm2_resid = cvxpy.norm2(ln_conc - np.log(c_mid))
        if max_mdf is not None and min_mdf is not None:
            constraints.append(cvxpy.leq(motive_force_lb, max_mdf))
            constraints.append(cvxpy.geq(motive_force_lb, min_mdf))
            objective = cvxpy.minimize(norm2_resid)
        elif max_mdf is not None:
            constraints.append(cvxpy.leq(motive_force_lb, max_mdf))
            objective = cvxpy.minimize(norm2_resid)
        elif min_mdf is not None:
            constraints.append(cvxpy.geq(motive_force_lb, min_mdf))
            objective = cvxpy.minimize(norm2_resid)
        else:
            objective = cvxpy.minimize(motive_force_lb + norm2_resid)

        program = cvxpy.program(objective, constraints)
        program.solve(quiet=True)
        return ln_conc.value, program.objective.value

    def FindMDF_OptimizeConcentrations(self, c_range=1e-3,
                                       bounds=None, c_mid=1e-3):
        """Optimize concentrations at optimal pCr.
        
        Runs two rounds of optimization to find "optimal" concentrations
        at the optimal MDF. First finds the globally optimal MDF.
        Then minimizes the l2 norm of deviations of log concentrations
        from c_mid given the optimal MDF.

        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
            c_mid: the median concentration.
 
        Returns:
            A 3 tuple (dGfs, concentrations, MDF value).
        """
        _, opt_mdf = self._FindMDF(c_range, bounds)
        return self.FindMDF_Regularized(c_range, bounds, c_mid,
                                        max_mdf=opt_mdf)

    def _MakeMinimalFeasbileConcentrationProblem(self, bounds=None, c_range=(1e-6, 1e-2)):
        # Define and apply the constraints on the concentrations
        constraints = []
        
        # Define and apply the constraints on the concentrations
        ln_conc = cvxpy.variable(1, self.Nc, name='lnC')
        constraints += self._MakeLnConcentratonBounds(ln_conc, bounds=bounds,
                                                      c_range=c_range)

        # find the row vector describing the overall stoichiometry
        S = cvxpy.matrix(self.S)
        dg0r_primes = cvxpy.matrix(self.dG0_r_prime)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr in the direction of the flux.
        # Reactions with a flux of 0 must be in equilibrium. 
        for i in xrange(self.Nr):
            # if the dG0 is unknown, this reaction imposes no new constraints
            if np.isnan(self.dG0_r_prime[0, i]):
                continue
            
            curr_dgr = dg0r_primes[0, i] + default_RT * ln_conc * S[:, i]

            if self.fluxes[0, i] != 0:
                constraints.append(cvxpy.leq(curr_dgr * np.sign(self.fluxes[0, i]),
                                             self.DEFAULT_REACTION_UB))
                constraints.append(cvxpy.geq(curr_dgr * np.sign(self.fluxes[0, i]),
                                             self.DEFAULT_REACTION_LB))
            else:
                constraints.append(cvxpy.eq(curr_dgr, 0))
        
        # Set the constraints
        return ln_conc, constraints

    def FindMinimalFeasibleConcentration(self, index_to_minimize,
                                         bounds=None, c_range=(1e-6, 1e-2)):
        """
            Compute the smallest ratio between two concentrations which makes the pathway feasible.
            All other compounds except these two are constrained by 'bounds' or unconstrained at all.
        
            Arguments:
                index_to_minimize - the column index of the compound whose concentration 
                                    is to be minimized
        
            Returns:
                dGs, concentrations, target-concentration
        """
        ln_conc, constraints = self._MakeMinimalFeasbileConcentrationProblem(bounds, c_range)
        objective = cvxpy.minimize(ln_conc[index_to_minimize])
        program = cvxpy.program(objective, constraints)
        program.solve(quiet=True) 
        return ln_conc.value, program.objective.value

    def _MakeMinimumFeasbileConcentrationsProblem(self, bounds=None,
                                                  c_range=(1e-6, 1e-2)):
        """Creates the CVXOPT problem for finding minimum total concentrations.
        
        Returns:
            Two tuple (ln_concentrations var, problem).
        """
        assert self.dG0_f_prime is not None
        
        constraints = []
        
        # Define and apply the constraints on the concentrations
        ln_conc = cvxpy.variable(1, self.Nc, name='lnC')
        constraints += self._MakeLnConcentratonBounds(ln_conc, bounds=bounds,
                                                      c_range=c_range)
        
        # Make the objective and problem.
        S = cvxpy.matrix(self.S)
        
        # Make flux-based constraints on reaction free energies.
        # All reactions must have negative dGr in the direction of the flux.
        # Reactions with a flux of 0 must be in equilibrium.
        dgf_primes = default_RT * ln_conc + cvxpy.matrix(self.dG0_f_prime)
        for i in xrange(self.Nr):
            if self.fluxes[0, i] > 0:
                constraints.append(cvxpy.leq(S[i, :] * dgf_primes,
                                             self.DEFAULT_REACTION_UB))
                constraints.append(cvxpy.geq(S[i, :] * dgf_primes,
                                             self.DEFAULT_REACTION_LB))
            elif self.fluxes[0, i] == 0:
                constraints.append(cvxpy.eq(S[i, :] * dgf_primes,
                                            0))
            else:
                constraints.append(cvxpy.geq(S[i, :] * dgf_primes,
                                             -self.DEFAULT_REACTION_UB))
                constraints.append(cvxpy.leq(S[i, :] * dgf_primes,
                                             -self.DEFAULT_REACTION_LB))
        
        return ln_conc, constraints

    def FindMinimumFeasibleConcentrations(self, bounds=None):
        """Use the power of convex optimization!
        
        minimize sum (concentrations)
        
        we can do this by using ln(concentration) as variables and leveraging 
        the convexity of exponentials. 
        
        min sum (exp(ln(concentrations)))
        """
        assert self.dG0_f_prime is not None
        
        ln_conc, constraints = self._MakeMinimumFeasbileConcentrationsProblem(bounds=bounds)
        total_conc = cvxpy.sum(cvxpy.exp(ln_conc))
        program = cvxpy.program(cvxpy.minimize(total_conc), constraints)
        program.solve(quiet=True)
        return ln_conc.value, total_conc.value

    def FindKineticOptimum(self, bounds=None):
        """Use the power of convex optimization!
        
        minimize sum (protein cost)
        
        we can do this by using ln(concentration) as variables and leveraging 
        the convexity of exponentials.         
        """
        assert self.dG0_f_prime is not None
        
        ln_conc, constraints = self._MakeMinimumFeasbileConcentrationsProblem()
        total_inv_conc = cvxpy.sum(cvxpy.exp(-ln_conc))
        program = cvxpy.program(cvxpy.minimize(total_inv_conc), constraints)
        program.solve(quiet=True)
        return ln_conc.value, total_inv_conc.value

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
        reaction = KeggReaction(sparse, arrow='=>')
        return str(reaction)

    def FindMDF(self, normalization=DeltaGNormalization.DEFAULT):
        """Find the MDF (Maximal Thermodynamic Driving Force).
        
        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.
        
        Returns:
            A 3 tuple (optimal dGfs, optimal concentrations, optimal mdf).
        """
        return self._FindMDF(self.c_range, self.bounds, normalization)

    def GetTotalReactionEnergy(self, min_driving_force=0, maximize=True):
        """
            Maximizes the total pathway dG' (i.e. minimize energetic cost).
            Arguments:
                min_driving_force - the lower limit on each reaction's driving force
                                    (it is common to provide the optimize driving force
                                    in order to find the concentrations that minimize the
                                    cost, without affecting the MDF).
                maximize          - if True then finds the maximal total dG.
                                    if False then finds the minimal total dG.
        """
        ln_conc, constraints, total_g = self._GetTotalReactionEnergy(
                                self.c_range, self.bounds, min_driving_force)
        
        if maximize:
            objective = cvxpy.maximize(total_g)
        else:
            objective = cvxpy.minimize(total_g)
        
        program = cvxpy.program(objective, constraints)
        program.solve(quiet=True)
        return ln_conc.value, program.objective.value
       
    def FindMinimalFeasibleConcentration(self, cid_to_minimize):
        """
            Compute the smallest ratio between two concentrations which makes the pathway feasible.
            All other compounds except these two are constrained by 'bounds' or unconstrained at all.
        min_conc
            Arguments:
                cid - the CID of the compound whose concentration should be minimized
        
            Returns:
                dGs, concentrations, target-concentration
        """
        index = self.cids.index(cid_to_minimize)
        return Pathway.FindMinimalFeasibleConcentration(self, index,
                                                    self.bounds, self.c_range)

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
        plt.yticks(range(self.Nc, 0, -1), self.cids, fontproperties=FontProperties(size=8))
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

    def WriteResultsToHtmlTables(self, html_writer, concentrations):
        self.WriteConcentrationsToHtmlTable(html_writer, concentrations)
        self.WriteProfileToHtmlTables(html_writer, concentrations)

    def WriteConcentrationsToHtmlTable(self, html_writer, concentrations=None):
        #html_writer.write('<b>Compound Concentrations</b></br>\n')
        dict_list = []
        for c, cid in enumerate(self.cids):
            d = {}
            d['KEGG CID'] = cid
            lb, ub = self.GetConcentrationBounds(cid)
            d['Concentration LB [M]'] = '%.2e' % lb
            if concentrations is not None:
                d['Concentration [M]'] = '%.2e' % concentrations[0, c]
            d['Concentration UB [M]'] = '%.2e' % ub
            dict_list.append(d)
        if concentrations is not None:
            headers = ['KEGG CID', 'Concentration LB [M]',
                         'Concentration [M]', 'Concentration UB [M]']
        else:
            headers = ['KEGG CID', 'Concentration LB [M]',
                       'Concentration UB [M]']
        
        html_writer.write_table(dict_list, headers=headers)
    
    def WriteProfileToHtmlTable(self, html_writer, concentrations=None):
        #html_writer.write('<b>Biochemical Reaction Energies</b></br>\n')
        phys_concentrations = np.matrix(np.ones((1, len(self.cids)))) * self.DEFAULT_PHYSIOLOGICAL_CONC
        if 1 in self.cids:
            # C00001 (water) is an exception, its concentration is always set to 1
            phys_concentrations[0, self.cids.index(1)] = 1
        
        dG_r_prime_c = self.CalculateReactionEnergiesUsingConcentrations(phys_concentrations)
        headers=["reaction", 'formula', 'flux', 
                 symbol_dr_G0_prime + " [kJ/mol]",
                 symbol_dr_Gc_prime + \
                 " [kJ/mol] (%g M)" % self.DEFAULT_PHYSIOLOGICAL_CONC] 

        if concentrations is not None:
            dG_r_prime = self.CalculateReactionEnergiesUsingConcentrations(concentrations)
            headers.append(symbol_dr_G_prime + " [kJ/mol]")
            dG_mat = np.hstack([self.dG0_r_prime.T, dG_r_prime_c.T, dG_r_prime.T])
        else:
            dG_mat = np.hstack([self.dG0_r_prime.T, dG_r_prime_c.T])
        
        dict_list = []
        for r, rid in enumerate(self.rids):
            d = {}
            d['reaction'] = rid
            d['flux'] = "%g" % abs(self.fluxes[0, r])
            d['formula'] = self.GetReactionString(r, show_cids=False)

            flux_sign = np.sign(self.fluxes[0, r])
            d[headers[3]] = dG_mat[r, 0] * flux_sign
            d[headers[4]] = dG_mat[r, 1] * flux_sign
            if concentrations is not None:
                d[headers[5]] = dG_mat[r, 2] * flux_sign
            dict_list.append(d)
            
        d = {'reaction':'Total',
             'flux':'1',
             'formula':self.GetTotalReactionString(show_cids=False)}
        d[headers[3]] = float(self.dG0_r_prime * self.fluxes.T)
        d[headers[4]] = float(dG_r_prime_c * self.fluxes.T)
        if concentrations is not None:
            d[headers[5]] = float(dG_r_prime * self.fluxes.T)
        dict_list.append(d)
        html_writer.write_table(dict_list, headers=headers, decimal=1)
        
        html_writer.insert_toggle(start_here=True, label='Show CSV')
        for r, rid in enumerate(self.rids):
            html_writer.write('%g,' % self.fluxes[0, r] + 
                              ','.join(['%.1f' % x for x in dG_mat[r, :].flat]) + 
                              '</br>\n')
        html_writer.div_end()
        
if __name__ == '__main__':
    S = np.matrix("-1, 0, 0; 1, -1, 0; 0, 1, 1; 0, 0, -1")
    dGs = np.matrix([[0, 10, 12, 2]])
    fluxes = np.matrix([[1, 1, -1]])
    rids = [1, 2, 3]
    cids = [1, 2, 3, 4]
    #keggpath = KeggPathway(S, rids, fluxes, cids, dGs, c_range=(1e-6, 1e-3))
    #dGf, concentrations, protein_cost = keggpath.FindKineticOptimum()
    #print 'protein cost: %g protein units' % protein_cost
    #print 'concentrations: ', concentrations
    #print 'sum(concs): ', sum(concentrations), 'M'
    
    keggpath = KeggPathway(S, rids, fluxes, cids, dGs, c_range=(1e-6, 1e-3))
    concentrations, mdf = keggpath.FindMDF()
    print 'MDF: %g' % mdf
    print 'concentrations:', concentrations
    