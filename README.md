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