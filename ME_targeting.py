#Copyright 2014-2016 BioInformatics Research Center, KAIST
import copy
import glob
import logging
import os
import sys
import time
import warnings
from itertools import combinations
from multiprocessing import Process, Queue

import numpy as np
import pandas as pd

from cobra.io import read_sbml_model, write_sbml_model

from me_targeting import utils
from me_targeting.flux_analysis import Simulator
from me_targeting.flux_analysis import FSEOF
from me_targeting.flux_analysis import FVSEOF


def do_fba_simulation(cobra_model, biomass_reaction, 
                       target_reaction, carbon_reaction,
                       minimum_production_rate, constraints, 
                       each_reaction_set, targeting_mode, 
                       gene_rxn_dict, queue):
    
    

    flux_constraints = constraints

    obj = Simulator.Simulator()
    obj.load_cobra_model(cobra_model)
    model_status, objective, flux = obj.run_FBA(new_objective=target_reaction, mode='max')
    max_target_production_rate = flux[target_reaction]
    lb_target_production_rate = max_target_production_rate * minimum_production_rate
    flux_constraints[target_reaction] = [lb_target_production_rate, 1000.0]

    result_list = []
    for each_set in each_reaction_set:
        key_string = ';'.join(each_set)
        simulation_flux_constraints = copy.deepcopy(flux_constraints)
        for each_reaction in each_set:
            if targeting_mode=='reaction':
                simulation_flux_constraints[each_reaction] = [0.0, 0.0]
            elif targeting_mode=='gene':
                for unit_reaction in gene_rxn_dict[each_reaction]:
                    simulation_flux_constraints[unit_reaction] = [0.0, 0.0]
        try:
            _, _, perturbed_flux = obj.run_FBA(new_objective=biomass_reaction,
                                                flux_constraints=simulation_flux_constraints, 
                                                mode='max', 
                                                internal_flux_minimization=False)
                                                     
  
            biomass_flux = perturbed_flux[biomass_reaction]
            target_flux = perturbed_flux[target_reaction]
            c_source_flux = perturbed_flux[carbon_reaction]
            result_list.append([key_string, biomass_flux, target_flux, c_source_flux, simulation_flux_constraints])
        except:
            result_list.append([key_string, 'error', 'error', 'error', simulation_flux_constraints])
            
    queue.put(result_list)
    return


def run_FBA_targeting(output_dir, cobra_model, biomass_reaction, 
                       target_reaction, carbon_reaction, constraints, 
                       minimum_production_rate, targeting_mode='reaction', 
                       target_num=1, cpu_num=1):

    candidate_reactions = []
    candidate_genes = []
    gene_rxn_dict = {}
    for each_reaction in cobra_model.reactions:
        metabolites = each_reaction.reactants + each_reaction.products
        compartments = [each_metabolite.compartment for each_metabolite in metabolites]
        compartments = list(set(compartments))
        if len(compartments) == 1:
            if len(each_reaction.genes) > 0:
                candidate_reactions.append(each_reaction.id)
        for gene in each_reaction.genes:
            candidate_genes.append(gene.id)
            if gene.id not in gene_rxn_dict:
                gene_rxn_dict[gene.id] = [rxn.id for rxn in gene.reactions]
            else:
                for rxn in gene.reactions:
                    gene_rxn_dict[gene.id].append(rxn.id)
    candidate_genes = list(set(candidate_genes))


    combinatorial_sets = []
    if targeting_mode=='reaction':
        for each_reaction_set in combinations(candidate_reactions, target_num):
            combinatorial_sets.append(each_reaction_set)
    elif targeting_mode=='gene':
        for each_gene_set in combinations(candidate_genes, target_num):
            combinatorial_sets.append(each_gene_set)
    else:
        raise NotImplementedError

    reaction_sets = np.array_split(combinatorial_sets, cpu_num)
    procs = []
    queues = {}
    fba_results = {}
    count = 0
    for reaction_set in reaction_sets:
        queue = Queue()
        proc = Process(target=do_fba_simulation, 
                       args=(cobra_model, biomass_reaction, 
                             target_reaction, carbon_reaction,
                             minimum_production_rate, 
                             constraints, reaction_set,
                             targeting_mode, gene_rxn_dict,
                             queue))
        queues[count] = queue
        procs.append(proc)
        count += 1
    for proc in procs:
        proc.start()
    for num in queues:
        fba_results[num] = queues[num].get()
        queues[num].close()
    for proc in procs:
        proc.join()

    simulation_results = []
    for num in fba_results:
        for item in fba_results[num]:
            simulation_results.append(item)
    
    with open(output_dir + '/fba_knockout_result.txt', 'w') as fp:
        fp.write('reaction(s)\tbiomass flux\ttarget flux\tcarbon uptake flux\n')
        for item in simulation_results:
            fp.write('%s\t%s\t%s\t%s\n'%(item[0], str(item[1]), str(item[2]), str(item[3])))
    return



def run_FSEOF_targeting(output_dir, cobra_model, biomass_reaction, target_reaction):
    obj = FSEOF.FSEOF()
    obj.load_cobra_model(cobra_model)
    df = obj.run_FSEOF(biomass_reaction, target_reaction)
    df.to_csv(output_dir+'/fseof_result_df.csv')
    result_df = obj.result_summary()
    result_df.to_csv(output_dir+'/fseof_summary_result.csv')
    return


def run_FVSEOF_targeting(output_dir, cobra_model, biomass_reaction, target_reaction, cpu_num=1):
    obj = FVSEOF.FVSEOF()
    obj.load_cobra_model(cobra_model)
    ResultDF = obj.run_FVSEOF(biomass_reaction, target_reaction, {}, cpu_num=cpu_num)
    ResultDF.to_csv(output_dir+'/fvseof_result_df.csv')
    df = obj.result_summary()
    df.to_csv(output_dir+'/fvseof_result_summary.csv')
    up_df, down_df = obj.get_target_info(output_dir+'/fvseof_result_summary.csv')
    up_df.to_csv(output_dir+'/fvseof_up_targets.csv')
    down_df.to_csv(output_dir+'/fvseof_down_targets.csv')
    return


def main():
    start = time.time()
    warnings.filterwarnings("ignore")

    parser = utils.argument_parser()
    options = parser.parse_args()

    output_dir = options.output_dir
    model_file = options.model_input
    target_reaction = options.target_rxn
    bio_rxn = options.biomass_rxn
    carbon_reaction = options.carbon_rxn
    method = options.method
    cpu_num = options.cpu_num
    target_num = options.target_num
    target_mode = options.target_mode
    minimum_production_rate = options.minimum_production_rate


    if not os.path.exists(output_dir): # If the directory does not exist, make it
        os.makedirs(output_dir)

    cobra_model = read_sbml_model(model_file)

    constraints = {}

    #constraints[bio_rxn] = (0.0, 1000.0)
    #constraints['PYK'] = (0, 0)

    print('Simulation start: %s'%model_file)

    if method == 'fba':
        run_FBA_targeting(output_dir, cobra_model, bio_rxn, 
                           target_reaction, carbon_reaction, constraints,
                           minimum_production_rate=minimum_production_rate,
                           targeting_mode=target_mode, 
                           target_num=target_num, cpu_num=cpu_num)
    elif method == 'fvseof':
        run_FVSEOF_targeting(output_dir, cobra_model, bio_rxn, 
                             target_reaction, cpu_num)
    elif method == 'fseof':
        run_FSEOF_targeting(output_dir, cobra_model, bio_rxn, 
                            target_reaction)

    else:
        print('Wrong method input')
        raise NotImplementedError
    
    logging.info(time.strftime("Elapsed time %H:%M:%S", time.gmtime(time.time() - start)))

if __name__ == '__main__':
    main()
