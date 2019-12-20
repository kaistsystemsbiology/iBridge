# iBridge

#iBridge #Project This project is to identify target reactions or metabolites that will enhance production of commercial chemicals through a metabolite-centric approach


1. Clone the repository

        git clone https://anlito@bitbucket.org/anlito/ibridge.git

2. Create and activate virtual environment

        virtualenv venv
        source venv/bin/activate

    or you can use pre-built anaconda environemnt

        conda env create -f environment.yml
        conda activate engem_env

3. Install gurobipy

    In our case, we installed gurobipy in the root of a server, and created its symbolic link in venv:

        ln -s /usr/local/lib/python2.7/dist-packages/gurobipy/ ./venv/lib/python2.7/site-packages/


4. Change the directory

        cd engem

5. Install packages

    If you created the environment by conda, you can skip this step

        pip install pip --upgrade
        pip install -r requirements.txt


##Example

- Run modeling code with retrieving relavant data

        python MetScore.py