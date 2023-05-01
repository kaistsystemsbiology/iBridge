'''
Created on 2014. 6. 4.

@author: user
'''

import copy
import os

from gurobipy import *

# import Simulator
from me_targeting.flux_analysis import Simulator
import cobra.io as io

class FVA(Simulator.Simulator):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''

    def run_FVA(self, reaction_list = [], flux_constraints={}, inf_flag = False):
        print('Start FVA simulation ... ')

        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        fva_result_dic = {}
        for each_reaction in reaction_list:
            fva_result_dic[ each_reaction ] = {'min':0.0, 'max':0.0, 'model_stat':0.0}

        for each_obj_reaction in reaction_list:
            objective = each_obj_reaction

            Smatrix = self.Smatrix
            LowerFlux=self.Lower_Boundary_Constraint
            UpperFlux=self.Upper_Boundary_Constraint

            if inf_flag == False:
                for key in LowerFlux.keys():
                    if LowerFlux[key] == float("-inf"):
                        LowerFlux[key] = -1000.0

                for key in UpperFlux.keys():
                    if UpperFlux[key] == float("inf"):
                        UpperFlux[key] = 1000.0

            pairs, coffvalue = multidict( Smatrix )
            pairs = tuplelist(pairs)

            m = Model('FVA')
            m.setParam('OutputFlag', 0)
            m.reset()

            # create variables
            v = {}
            for each_reaction in model_reactions:
                if each_reaction in flux_constraints.keys():
                    v[each_reaction] = m.addVar(lb = flux_constraints[each_reaction][0], ub = flux_constraints[each_reaction][1] , name=each_reaction)
                else:
                    v[each_reaction] = m.addVar(lb = LowerFlux[each_reaction], ub = UpperFlux[each_reaction] , name=each_reaction)

            m.update()

            # Add constraints
            for each_metabolite in model_metabolites:
                if len( pairs.select(each_metabolite, '*') ) == 0:
                    continue
                m.addConstr( quicksum( v[reaction]*coffvalue[metabolite,reaction] for metabolite, reaction in pairs.select(each_metabolite, '*') ) == 0 )

            m.update()

            m.setObjective(v[objective], GRB.MAXIMIZE)
            m.optimize()

            max_model_stat = m.status
            print('Model stat %s'%m.status)
            max_flux = m.ObjVal

            m.setObjective(v[objective], GRB.MINIMIZE)
            m.optimize()

            min_model_stat = m.status
            if m.status != 2:
                print(max_model_stat, min_model_stat, each_obj_reaction)
                #raw_input('temp')


            if max_model_stat == 2 and min_model_stat == 2:
                min_flux = m.ObjVal
                fva_result_dic[ each_obj_reaction ] = {'min':min_flux, 'max':max_flux, 'model_stat':2}

        return fva_result_dic


if __name__ == '__main__':
    import time
    obj = FVA( )
    #obj.readModel( './Recon2_Human_Model.xml')
    (model_metabolites, model_reactions, Smatrix, Lower_Boundary_Constraint, Upper_Boundary_Constraint, objective_reaction) = obj.readModel( './iAF1260_SBML_export.xml')
    fva_result_dic = obj.run_FVA( model_reactions[0:3]  )
    print(fva_result_dic)
