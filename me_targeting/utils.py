import os
import argparse


def argument_parser(version=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--model_input', required=True, help='Input model file directory')
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    parser.add_argument('-t', '--target_rxn', required=True, help='Target reaction ID')
    parser.add_argument('-b', '--biomass_rxn', default='BIOMASS_Ec_iML1515_core_75p37M',
                        required=False, help='BIOMASS reaction ID')
    parser.add_argument('-c', '--carbon_rxn', default='EX_glc__D_e', 
                        required=False, help='Carbon source reaction ID')
    parser.add_argument('-n', '--cpu_num', type=int, default=1, \
                        required=False, help='Number of cpus to use (default: 1')
    parser.add_argument('-m', '--method', required=True, help='Simulation method: moma or fvseof')
    parser.add_argument('-p', '--percentage', required=False, help='Knockdown percentage')
    parser.add_argument('-target_num', '--target_num', type=int, default=1, \
                        required=False, help='Number of knockout targets to use (default: 1')
    parser.add_argument('-target_mode', '--target_mode', default='reaction', \
                        required=False, help='MOMA mode (reaction or gene)')
    parser.add_argument('-min_rate', '--minimum_production_rate', type=float, default=0.1, \
                        required=False, help='Required minimum production rate for MOMA')
    return parser
