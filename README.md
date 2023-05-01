# iBridge

#iBridge #Project This project is to identify target reactions or metabolites that will enhance production of commercial chemicals through a metabolite-centric approach


##Example

1. Clone the repository

        git clone https://bitbucket.org/kaistsystemsbiology/ibridge.git

2. Create and activate virtual environment

        virtualenv venv
        source venv/bin/activate


3. Install gurobipy

    In our case, we installed gurobipy in the root of a server, and created its symbolic link in venv:

        ln -s /usr/local/lib/python2.7/dist-packages/gurobipy/ ./venv/lib/python2.7/site-packages/


4. Change the directory

        cd ibridge

5. Install packages

    If you created the environment by conda, you can skip this step

        pip install pip --upgrade
        pip install -r requirements.txt


6. Or, install conda environment

        conda env create -f environment.yml
        conda activate ibridge
        
        
        
##FSEOF, FVSEOF, OptForce simulation

1. Create and activate virtual environment

        conda env create -f environment2.yml
        conda activate metool
        
2. Run FVSEOF

        python ME_targeting.py -i ./input/ijo1366_IRR_indirubin.xml -o ./output/fvseof -t EX_INDIRUBIN_LPAREN_e_RPAREN_ -b Ec_biomass_iJO1366_core_53p95M -m fvseof
        
3. Run FSEOF

        python ME_targeting.py -i ./input/ijo1366_IRR_indirubin.xml -o ./output/fvseof -t EX_INDIRUBIN_LPAREN_e_RPAREN_ -b Ec_biomass_iJO1366_core_53p95M -m fseof
        

4. Run OptForce

        python OptForce.py