# iBridge

#iBridge #Project This project is to identify target reactions or metabolites that will enhance production of commercial chemicals through a metabolite-centric approach


##Example

1. Clone the repository

        git clone https://github.com/kaistsystemsbiology/iBridge.git

2. Change the directory

        cd ibridge

3. Create and activate virtual environment

        conda env create -f environment.yml
        conda activate ibridge

4. Install the virtual environment kerenl into the jupyter

        python -m ipykernel install —user —name ibridge —display-name "ibridge"
        
        
        
# FSEOF, FVSEOF, OptForce simulation

1. Create and activate virtual environment

        conda env create -f environment2.yml
        conda activate metool
        
2. Run FVSEOF

        python ME_targeting.py -i ./input/ijo1366_IRR_indirubin.xml -o ./output/fvseof -t EX_INDIRUBIN_LPAREN_e_RPAREN_ -b Ec_biomass_iJO1366_core_53p95M -m fvseof
        
3. Run FSEOF

        python ME_targeting.py -i ./input/ijo1366_IRR_indirubin.xml -o ./output/fvseof -t EX_INDIRUBIN_LPAREN_e_RPAREN_ -b Ec_biomass_iJO1366_core_53p95M -m fseof
        

4. Run OptForce

        python OptForce.py
