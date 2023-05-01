import copy
from itertools import combinations
from multiprocessing import Process, Queue

import numpy as np
from pandas import DataFrame
from pandas.io.parsers import read_csv
from scipy import stats

# import Simulator
from me_targeting.flux_analysis import Simulator
from gurobipy import *

class FVSEOF(Simulator.Simulator):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.fvseof_df = None

    def run_FVSEOF(self, biomass_rxn, target_rxn, flux_constraints={}, cpu_num=1):
        cobra_model = self.cobra_model
        all_reactions = [each_rxn.id for each_rxn in cobra_model.reactions]

        const_dic = flux_constraints

        a,b,c = self.run_FBA(new_objective=target_rxn, flux_constraints=const_dic, mode='min')
        r1_min_value = c[target_rxn]

        a,b,c = self.run_FBA(new_objective=target_rxn, flux_constraints=const_dic, mode='max')
        r1_max_value = c[target_rxn]

        steps = 10
        manipulation_ratio_list = np.linspace(r1_min_value, r1_max_value, steps)
        n_step = 1

        cnt = 1
        flux_info = {}
        for each_target_flux in manipulation_ratio_list:
            flux_const = {}
            flux_const = copy.deepcopy(const_dic)
            flux_const[target_rxn] = [each_target_flux, each_target_flux]

            a, b, c = self.run_FBA(new_objective=biomass_rxn, flux_constraints=flux_const, mode='max')
            growth_rate = c[biomass_rxn]

            flux_const2 = {}
            flux_const2 = copy.deepcopy(flux_const)
            flux_const2[biomass_rxn] = [growth_rate * 0.95, growth_rate * 0.95]

            reaction_sets = np.array_split(all_reactions, cpu_num)
            procs = []
            queues = {}
            fvseof_results = {}
            count = 0
            for reaction_set in reaction_sets:
                queue = Queue()
                proc = Process(target=self.fast_optimize, args=(reaction_set, target_rxn, each_target_flux, biomass_rxn, flux_const2, n_step, queue))
                queues[count] = queue
                procs.append(proc)
                count+=1
            for proc in procs:
                proc.start()
            for num in queues:
                fvseof_results[num] = queues[num].get()
                queues[num].close()
            for proc in procs:
                proc.join()

            for num in fvseof_results:
                tmp_result = fvseof_results[num]
                for each_rxn in tmp_result:
                    flux_info[cnt] = tmp_result[each_rxn]
                    cnt += 1
            n_step += 1  # n_step means each enforcement step

        df = DataFrame.from_dict(flux_info)
        self.fvseof_df = copy.deepcopy(df)

        return df.T


    def fast_optimize(self, reaction_set, target_rxn, each_target_flux, biomass_rxn, flux_const, n_step, queue):
        flux_info = {}
        for each_rxn in reaction_set:
            a,b,c = self.run_FBA(new_objective=each_rxn, flux_constraints=flux_const, mode='min')
            if c:
                min_value = c[each_rxn]
            else:
                min_value = -999999
            a,b,c = self.run_FBA(new_objective=each_rxn, flux_constraints=flux_const, mode='max')
            if c:
                max_value = c[each_rxn]
                biomass = c[biomass_rxn]
            else:
                max_value = 999999
                biomass = 999999

            each_rxn_min_value = min_value
            each_rxn_max_value = max_value

            v_avg = (each_rxn_max_value + each_rxn_min_value) / 2.0
            l_sol = abs(each_rxn_max_value - each_rxn_min_value)

            b = {}
            b['TargetRxn'] = target_rxn
            b['Step'] = n_step
            b['TargetFlux'] = each_target_flux
            b['Biomass'] = biomass
            b['Reaction'] = each_rxn
            b['MinFlux'] = each_rxn_min_value
            b['MaxFlux'] = each_rxn_max_value
            b['AvgFlux'] = v_avg
            b['SolFlux'] = l_sol
            flux_info[each_rxn] = b
        queue.put(flux_info)
        return


            


    def result_summary(self):
        df = self.fvseof_df.T
        target_reaction = list(df['TargetRxn'])[0]

        result_summarty_info = {}

        for each_reaction in list(df['Reaction']):
            result_summarty_info[each_reaction]={}
            temp_info = {}

            avg_flux = list(df[df['Reaction']==each_reaction]['AvgFlux'])
            target_flux = list(df[df['Reaction']==each_reaction]['TargetFlux'])
            min_flux = list(df[df['Reaction']==each_reaction]['MinFlux'])
            max_flux = list(df[df['Reaction']==each_reaction]['MaxFlux'])
            average_sol_range = list(df[df['Reaction']==each_reaction]['SolFlux'])

            mean_sol = np.mean(average_sol_range)

            slope_avg_flux, avg_intercept, avg_r_value, avg_p_value, avg_std_err = stats.linregress(target_flux, np.abs(avg_flux))
            slope_min_flux, min_intercept, min_r_value, min_p_value, min_std_err = stats.linregress(target_flux, np.abs(min_flux))
            slope_max_flux, max_intercept, max_r_value, max_p_value, max_std_err = stats.linregress(target_flux, np.abs(max_flux))


            temp_info['TARGET_REACTION'] = target_reaction
            # temp_info['REACTION'] = each_reaction
            temp_info['AVG_SLOPE_FLUX'] = slope_avg_flux
            temp_info['AVG_INTERCEPT'] = avg_intercept
            temp_info['AVG_R'] = avg_r_value
            temp_info['AVG_P'] = avg_p_value
            temp_info['AVG_STDERR'] = avg_std_err

            temp_info['MIN_SLOPE_FLUX'] = slope_min_flux
            temp_info['MIN_INTERCEPT'] = min_intercept
            temp_info['MIN_R'] = min_r_value
            temp_info['MIN_P'] = min_p_value
            temp_info['MIN_STDERR'] = min_std_err

            temp_info['MAX_SLOPE_FLUX'] = slope_max_flux
            temp_info['MAX_INTERCEPT'] = max_intercept
            temp_info['MAX_R'] = max_r_value
            temp_info['MAX_P'] = max_p_value
            temp_info['MAX_STDERR'] = max_std_err
            temp_info['MEAN_SOLUTION_RANGE'] = mean_sol
            result_summarty_info[each_reaction]=temp_info

        df = DataFrame.from_dict(result_summarty_info)
        self.sum_df = df
        return df.T

    def get_target_info(self,df_file):
        df = read_csv(df_file, index_col=0)
        new_df = df[df['MIN_P']<0.05]
        new_df = new_df[new_df['MIN_R']>0.9]
        new_df = new_df[new_df['AVG_P']<0.05]
        new_df = new_df[new_df['AVG_R']>0.9]
        new_df = new_df[new_df['MAX_P']<0.05]
        new_df = new_df[new_df['MAX_R']>0.9]

        up_df = new_df

        # new_df.to_csv('./Up_Target.csv')

        df = read_csv(df_file, index_col=0)
        new_df = df[df['MIN_P']<0.05]
        new_df = new_df[new_df['MIN_R']<-0.9]
        new_df = new_df[new_df['AVG_P']<0.05]
        new_df = new_df[new_df['AVG_R']<-0.9]
        new_df = new_df[new_df['MAX_P']<0.05]
        new_df = new_df[new_df['MAX_R']<-0.9]

        # new_df.to_csv('./Down_Target.csv')
        down_df = new_df

        return up_df, down_df