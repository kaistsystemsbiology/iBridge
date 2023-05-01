from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
import os

def run_FVA(model, target_rxn, target_const_rate, flux_constraints={}):
    with model as m:
        m.objective = target_rxn
        m.objective_direction = 'max'
        for rxn_id in flux_constraints:
            m.reactions.get_by_id(rxn_id).bounds = flux_constraints[rxn_id]
        target_max = m.slim_optimize()

    fva_result = {}
    for rxn in model.reactions:
        with model as m:
            tmp = target_const_rate * target_max
            m.reactions.get_by_id(target_rxn).bounds = (tmp, tmp)

            for rxn_id in flux_constraints:
                m.reactions.get_by_id(rxn_id).bounds = flux_constraints[rxn_id]

            m.objective = rxn.id

            m.objective_direction = 'min'
            min_flux = m.slim_optimize()
            m.objective_direction = 'max'
            max_flux = m.slim_optimize()
            fva_result[rxn.id] = (min_flux, max_flux)
    return fva_result


def compare_FV_range(bio_result, target_result):
    down_targets = {}
    up_targets = {}
    for rxn_id in target_result:
        target_lb, target_ub = target_result[rxn_id]
        bio_lb, bio_ub = bio_result[rxn_id]
        if target_lb > bio_ub:
            up_targets[rxn_id] = {'bio':(bio_lb, bio_ub),
                                  'target':(target_lb, target_ub)}
        elif target_ub < bio_lb:
            down_targets[rxn_id] = {'bio':(bio_lb, bio_ub),
                                    'target':(target_lb, target_ub)}
    return down_targets, up_targets


def recordResults(target_rxn, down_targets, up_targets):
    if not os.path.exists('./output'): # If the directory does not exist, make it
        os.makedirs('./output')
        
    with open(f'./output/OptForce_{target_rxn}.txt', 'w') as fp:
        fp.write(f'Target reaction: {target_rxn}\n')
        fp.write('\nDown targets\tBio lb\tBio ub\tTarget lb\tTarget ub\n')
        for rxn_id in down_targets:
            bio_lb, bio_ub = down_targets[rxn_id]['bio']
            target_lb, target_ub = down_targets[rxn_id]['target']
            fp.write(f'{rxn_id}\t{bio_lb:4f}\t{bio_ub:4f}\t{target_lb:4f}\t{target_ub:4f}\n')
        fp.write('\nUp targets\tBio lb\tBio ub\tTarget lb\tTarget ub\n')
        for rxn_id in up_targets:
            bio_lb, bio_ub = up_targets[rxn_id]['bio']
            target_lb, target_ub = up_targets[rxn_id]['target']
            fp.write(f'{rxn_id}\t{bio_lb:4f}\t{bio_ub:4f}\t{target_lb:4f}\t{target_ub:4f}\n')


        
if __name__ == '__main__':
    model_dir = './input/ijo1366_IRR_indirubin.xml'
    biomass_rxn = 'Ec_biomass_iJO1366_core_53p95M'
    target_rxn = 'EX_INDIRUBIN_LPAREN_e_RPAREN_'

    biomass_const = 0.9
    target_const = 0.9

    model = read_sbml_model(model_dir)
    with model as m:
        bio_max = m.slim_optimize()

    bio_fva_result = run_FVA(model, biomass_rxn, biomass_const)

    ###
    flux_constraint = {biomass_rxn:(0.1*bio_max, 0.1*bio_max)}
    ###

    target_fva_result = run_FVA(model, target_rxn, target_const, flux_constraints=flux_constraint)

    down_targets, up_targets = compare_FV_range(bio_fva_result, target_fva_result)
    recordResults(target_rxn, down_targets, up_targets)
