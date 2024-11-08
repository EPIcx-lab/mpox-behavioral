# -------------------------------------------------------------------------------------------------------------------------
# ** REMARK ABOUT REM_LIST ****************************************************************************************
# -------------------------------------------------------------------------------------------------------------------------
# if analysis_type == "Targeted_Mixed_Removal" or == "Random_Mixed_Removal", then the values in rem_list
# correspond to the % of nodes that spontaneously change behavior. 
# if instead analysis_type == "Contactbased_link_removal_general", then the values in rem_list correspond
# to the % of nodes who got in contact when a detected in the previous back_in_time days that change behavior
# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------

import subprocess
import sys
import numpy as np
from multiprocessing import Pool, cpu_count
import os
import time
import pandas as pd
from datetime import timedelta

# --------------------------------------------------------------------------------------------------------------
# ** PARAMETERS DEFS *********************************************************************************** 
# --------------------------------------------------------------------------------------------------------------
analysis_type = "Targeted_Mixed_Removal"
subfolder = analysis_type

net_list = ["N1","N2","N3","N4","N5"]
parentdir = '/'.join(os.path.abspath('./').split(os.sep)[:-6])
maindir = f'{parentdir}/Dropbox/DM/INSERM/MONKEYPOX'
net_folder = 'FIVE_NETWORKS_threshold'
exe_file_name = "mine.exe"

start_simulation_date_default = pd.to_datetime('5-7-2022') 
end_simulation_date = pd.to_datetime('8-31-2022') 
start_vaccines_date = pd.to_datetime('5-27-2022')
end_vaccines_date = pd.to_datetime('7-10-2022')
start_firstdose_date = pd.to_datetime('7-11-2022')
interrupt_reference_date = pd.to_datetime('6-30-2022')

# float ------------------------------------------------------
mu_1_before_start = 1/8.82
mu_1_May = 1/8.82 
mu_1_June = 1/6.71
mu_1_after_June = 1/6.71 
mu  = 1.0/14.0
epsilon = 1.0/8.0
VES_pep = 0.89
VEI_pep = 0.0
VEE_pep = 0.0
VES_smallpox = 0.71
VEI_smallpox = 0.0
VEE_smallpox = 0.0
VES_firstdose = 0.78
VEI_firstdose = 0.0
VEE_firstdose = 0.0
first_age_to_immunize = 43
efficacy_delay_pep = 14   
efficacy_delay_prep = 14   
exposure_vaccination_delay = 0
T_data = 180
n_runs = 50
n_initial_I = 10
degree_vaccination_threshold = 10
prep_vaccination = 1   # 1 is true
msm_population = 65_000
verbose = 1  
save_state = 1         # 0 is false

class No_removal_parameters:
    def __init__(self):
        # zipped
        self.p_detection_list = [.2,.5,.8]
        self.beta_q_list = [.42,.32,.28]
        self.beta_list = list(np.array(self.beta_q_list) * 1)
        self.start_simulation_date_list = pd.to_datetime(['5-6-2022','5-8-2022','5-7-2022'])

        # non zipped 
        self.total_doses_to_be_given_list = [10_000]

    def default_parameters(self, analysis_type):
        if analysis_type == "PY_ages_runs_vaccines_CORRECTED":
            self.back_in_time = -100
            self.prem = 0
            self.rem = 0
            self.start_degree_date = pd.to_datetime('6-30-2022')
            self.end_degree_date = pd.to_datetime('7-1-2022')
            self.start_rem_date = pd.to_datetime('6-30-2022') 
            self.end_rem_date = pd.to_datetime('7-1-2022')   
            self.only_non_vaccinated_change_behavior = -1 # 1 is true. Here no b. changes -> default=-1
            self.save_weights = 1

    def debugging(self):
        if(len(self.beta_q_list) != len(self.p_detection_list)):
            sys.exit("wrong zipped params length")
        if(len(self.beta_q_list) != len(self.start_simulation_date_list)):
             sys.exit("wrong zipped params length")
        if(len(self.beta_q_list) != len(self.beta_list)):
             sys.exit("wrong zipped params length")

class Removal_parameters:
    def __init__(self):
        self.total_doses_to_be_given = 10_000
        self.settings = [2]    
        self.back_in_time = 14
        self.save_weights = 1     #0 is false. weights file is very heavy
    
    def default_parameters(self, analysis_type):
        if analysis_type == "Targeted_Mixed_Removal" or analysis_type == "Random_Mixed_Removal":
            self.back_in_time = -100

class Setting_Parameters:
    def __init__(self, setting):
        self.setting = setting
    
    def set_parameters(self, start_simulation_date_default):
        if self.setting == 1:
            self.p_detection = 0.2
            self.beta_q = 0.42
            self.beta = self.beta_q * 1         
            self.lag = -1
            self.prem_list = [77]
            self.rem_list = [18]
            self.start_simulation_date = start_simulation_date_default + timedelta(days=self.lag)
            self.start_degree_date = self.start_simulation_date
            self.start_rem_date_list = pd.to_datetime(['6-18-2022']) #pd.to_datetime(['6-15-2022','6-18-2022','6-21-2022']) 
            self.end_rem_date_list = pd.to_datetime(['7-18-2022']) #pd.to_datetime(['7-15-2022','7-18-2022','7-21-2022'])   
            self.end_degree_date_list = self.end_rem_date_list
            self.only_non_vaccinated_change_behavior = 0
            
        elif self.setting == 2:
            self.p_detection = 0.5
            self.beta_q = 0.29  
            self.beta = self.beta_q * 1    # you should have beta < beta_q (non detected are less infectious)            
            self.lag = 0
            self.prem_list = [95]
            self.rem_list = [48]
            self.start_simulation_date = start_simulation_date_default + timedelta(days=self.lag)
            self.start_degree_date = self.start_simulation_date
            self.start_rem_date_list = pd.to_datetime(['6-18-2022']) #pd.to_datetime(['6-15-2022','6-18-2022','6-21-2022']) 
            self.end_rem_date_list = pd.to_datetime(['7-16-2022']) #pd.to_datetime(['7-15-2022','7-18-2022','7-21-2022'])   
            self.end_degree_date_list = self.end_rem_date_list
            self.only_non_vaccinated_change_behavior = 0
            
        elif self.setting == 3:
            self.p_detection = 0.8
            self.beta_q = 0.28     
            self.beta = self.beta_q * 1           
            self.lag = 0
            self.prem_list = [62]
            self.rem_list = [17]
            self.start_simulation_date = start_simulation_date_default + timedelta(days=self.lag)
            self.start_degree_date = self.start_simulation_date
            self.start_rem_date_list = pd.to_datetime(['6-21-2022']) #pd.to_datetime(['6-15-2022','6-18-2022','6-21-2022']) 
            self.end_rem_date_list = pd.to_datetime(['7-21-2022']) #pd.to_datetime(['7-15-2022','7-18-2022','7-21-2022'])   
            self.end_degree_date_list = self.end_rem_date_list
            self.only_non_vaccinated_change_behavior = 0

        else:
            sys.exit("ERROR: Non-valid setting \n")

# --------------------------------------------------------------------------------------------------------------
# ** FUNCTIONS DEFS *********************************************************************************** 
# --------------------------------------------------------------------------------------------------------------
def define_analysis_code(analysis_type):
    if analysis_type == "PY_ages_runs_vaccines_CORRECTED":
        analysis_code = 1
    elif analysis_type == "Targeted_Mixed_Removal":
        analysis_code = 2
    elif analysis_type == "Random_Mixed_Removal":
        analysis_code = 3
    elif analysis_type == "Contactbased_link_removal_general":
        analysis_code = 4
    else:
        sys.exit("Got invalid analysis type \n")
    return analysis_code

def define_paths(maindir, net_folder, net, subfolder):
    ids_filename = f"IDS_{net}.txt"
    ages_filename = f"ages_{net}.txt"
    net_filename = f'NET_{net}.txt'

    indir = maindir + f"/DATA/Networks/{net_folder}/"
    outdir_results = maindir + f"/OUTPUTS/{net_folder}/{net}/{subfolder}/results/"
    outdir_states = maindir + f"/OUTPUTS/{net_folder}/{net}/{subfolder}/states/"
    outdir_weights = maindir + f"/OUTPUTS/{net_folder}/{net}/{subfolder}/weights/"
    
    return ids_filename, ages_filename, net_filename, indir, outdir_results, outdir_states, outdir_weights

def Time_Related_Quantities(start_simulation_date, end_simulation_date, mu_1_May, mu_1_June, mu_1_after_June):
    May_days = (pd.to_datetime('5-31-2022') - start_simulation_date).days + 1
    After_June_days = (end_simulation_date - pd.to_datetime('6-30-2022')).days
    mu_1_vector = [mu_1_May]*May_days + [mu_1_June]*30 + [mu_1_after_June]*After_June_days
    T_simul = May_days + 30 + After_June_days
    return May_days, After_June_days, mu_1_vector, T_simul

def From_date_to_day(date, start_simulation_date):
    return ((date - start_simulation_date).days + 1)
