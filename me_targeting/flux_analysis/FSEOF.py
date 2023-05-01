'''
Created on 2014. 6. 4.

@author: user
'''

from pandas import DataFrame
from pandas.io.parsers import read_csv
import numpy as np
from scipy import stats

from me_targeting.flux_analysis import Simulator
# import Simulator


class FSEOF(Simulator.Simulator):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.target_info_dic = {}
        self.wild_opt_flux = {}

    def run_FSEOF(self, biomass_rxn, target_rxn):
        moma_obj = Simulator.Simulator()
        moma_obj.load_cobra_model(self.cobra_model)
        constdic = {}

        fba_obj = Simulator.Simulator()
        fba_obj.load_cobra_model(self.cobra_model)
        a, b, flux_dic = fba_obj.run_FBA(new_objective=biomass_rxn, internal_flux_minimization=True)
        wild_opt_flux = flux_dic

        biomass_max = flux_dic[biomass_rxn]

        new_objective_target = target_rxn

        a, b, flux_dic = fba_obj.run_FBA(new_objective=new_objective_target, flux_constraints=constdic, mode='max', internal_flux_minimization=True)
        max_const = flux_dic[new_objective_target]

        a, b, flux_dic = fba_obj.run_FBA(new_objective=new_objective_target, flux_constraints=constdic, mode='min', internal_flux_minimization=True)
        min_const = flux_dic[new_objective_target]

        print('Target minimum flux : %s. Target maximum flux: %s.'%(min_const, max_const))

        flux_dist_dic_set = {}
        count = 1
        for each_flux_const in np.linspace(min_const, max_const*0.8, 10):  # No. of steps
        #for each_flux_const in np.linspace(min_const, 0.01713, 10):  # No. of steps
            constdic = {}
            count += 1
            #constdic[biomass_rxn] = [biomass_max, 1000.0]
            constdic[new_objective_target] = [each_flux_const, each_flux_const]

            stat, b, flux_dic = moma_obj.run_FBA(new_objective=biomass_rxn, flux_constraints=constdic, internal_flux_minimization=True)

            print('Status : %s\t Target production : %s \t Biomass max : %s' % (
                a, flux_dic[new_objective_target], flux_dic[biomass_rxn]))
            flux_dist_dic_set[each_flux_const] = flux_dic


        df = DataFrame.from_dict(flux_dist_dic_set)
        self.fseof_df = df
        return df

    def result_summary(self):
        rm_reactions=[]
        all_reactions = [reaction.id for reaction in self.cobra_model.reactions]
        for each_reaction in self.cobra_model.reactions:
            all_metabolite_compartments = [met.compartment for met in each_reaction.reactants+each_reaction.products]
            all_metabolite_compartments = list(set(all_metabolite_compartments))
            #if len(all_metabolite_compartments) > 1:
            #    rm_reactions.append(each_reaction.id)
            #if len(all_metabolite_compartments) == 1:
            #    if all_metabolite_compartments[0] == 'e':
            #        rm_reactions.append(each_reaction.id)
            #if len(each_reaction.get_gene()) == 0:
            #    rm_reactions.append(each_reaction.id)

        selected_reactions = set(all_reactions).difference(set(rm_reactions))

        df = self.fseof_df.T
        df.to_csv('temp.csv')
        df = read_csv('temp.csv')
        corr = df.abs().corr()
        cov = df.abs().cov()

        p_corr_df = corr[corr['Unnamed: 0'] >= 0.1]
        n_corr_df = corr[corr['Unnamed: 0'] <= -0.1]

        new_cov = cov.sort_values('Unnamed: 0', ascending=False)
        top_new_cov = new_cov['Unnamed: 0']
        new_cov = cov.sort_values('Unnamed: 0', ascending=True)
        down_new_cov = new_cov['Unnamed: 0']

        up_target_reactions = list(set(p_corr_df.index) & set(top_new_cov.index) & set(selected_reactions))
        down_target_reactions = list(set(n_corr_df.index) & set(down_new_cov.index) & set(selected_reactions))


        summary_target_info = {}
        for each_reaction in up_target_reactions:
            summary_target_info[each_reaction]={}
            target_corr = corr[each_reaction].loc['Unnamed: 0']
            target_cov = cov[each_reaction].loc['Unnamed: 0']
            summary_target_info[each_reaction]['pearson'] = target_corr
            summary_target_info[each_reaction]['covariance'] = target_cov
            summary_target_info[each_reaction]['status'] = 'UP'

        for each_reaction in down_target_reactions:
            summary_target_info[each_reaction]={}
            target_corr = corr[each_reaction].loc['Unnamed: 0']
            target_cov = cov[each_reaction].loc['Unnamed: 0']
            summary_target_info[each_reaction]['pearson'] = target_corr
            summary_target_info[each_reaction]['covariance'] = target_cov
            summary_target_info[each_reaction]['status'] = 'DOWN'
        df = DataFrame.from_dict(summary_target_info)
        return df.T
