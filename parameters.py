# -------------------------------------------------------------------------------------------------------------------------
# ** DETAILS ****************************************************************************************
# -------------------------------------------------------------------------------------------------------------------------
'''
Here detailed information about some of the input parameters are provided

- Analysis type
Possible values: "PY_ages_runs_vaccines_CORRECTED" -> no behavioral changes
                 "Random_Mixed_Removal" -> Widespread behavioral changes. Each MSM has the same probability of changing behavior
                 "Targeted_Mixed_Removal" -> Behavioral changes occurring in high-risk MSM preferentially. Probability of changing behavior proportional to node's degree
                 "Contactbased_link_removal_general" Cases' contacts change behavior

- Back_in_time: Applicable only for "Contactbased_link_removal_general". If >0, when an S individual is infected, some of its contacts occured in the 
                previous 'back_in_time' days change behavior. If=-1, contacts are searched since the beginning of the simulation.

- start_degree_date, end_degree_date: Nodes' degree is computed aggregating the temporal network from start_degree_date to end_degree_date

'''
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
analysis_type = "Targeted_Mixed_Removal"                             # see details
subfolder = analysis_type                                            # leave it as is

net_list = ["N1","N2","N3","N4","N5"]                                #name of input network files (default extension is .txt)
parentdir = '/'.join(os.path.abspath('./').split(os.sep)[:-6])
maindir = f'{parentdir}/Dropbox/DM/INSERM/MONKEYPOX'
net_folder = 'FIVE_NETWORKS_threshold'
exe_file_name = "mine.exe"

start_simulation_date_default = pd.to_datetime('5-7-2022')           # starting date of the simulation
end_simulation_date = pd.to_datetime('8-31-2022')                    # ending date of the simulation
start_vaccines_date = pd.to_datetime('5-27-2022')                    # starting date of the PEP vaccination
end_vaccines_date = pd.to_datetime('7-10-2022')                      # ending date of the PEP vaccination
start_firstdose_date = pd.to_datetime('7-11-2022')                   # starting date of the PrEP vaccination
interrupt_reference_date = pd.to_datetime('6-30-2022')               # ending date of the PrEP vaccination

# float ------------------------------------------------------
mu_1_before_start = 1/8.82                                           # Inverse of the onset-to-testing period, before May
mu_1_May = 1/8.82                                                    # Inverse of the onset-to-testing period, in May
mu_1_June = 1/6.71                                                   # Inverse of the onset-to-testing period, in June
mu_1_after_June = 1/6.71                                             # Inverse of the onset-to-testing period, after June
mu  = 1.0/14.0                                                       # Inverse of the infectious period
epsilon = 1.0/8.0                                                    # Inverse of the incubation period
VES_pep = 0.89                                                       # PEP effectiveness against infection
VEI_pep = 0.0                                                        # PEP effectiveness against transmission
VEE_pep = 0.0                                                        # PEP effectiveness against onset of symptoms (never used, leave 0)
VES_smallpox = 0.71                                                  # Smallpox vaccine effectiveness against infection
VEI_smallpox = 0.0                                                   # Smallpox vaccine effectiveness against transmission
VEE_smallpox = 0.0                                                   # Smallpox vaccine effectiveness against onset of symptoms (never used, leave 0)
VES_firstdose = 0.78                                                 # PrEP effectiveness against infection
VEI_firstdose = 0.0                                                  # PrEP effectiveness against transmission
VEE_firstdose = 0.0                                                  # PrEP effectiveness against onset of symptoms (never used, leave 0)
first_age_to_immunize = 43                                           # minimum age of MSM to whom assign smallpox vaccination                  
efficacy_delay_pep = 14                                              # delay between PEP administration and PEP becoming effective (effectiveness is 0 in the meanwhile)            
efficacy_delay_prep = 14                                             # delay between PrEP administration and PEP becoming effective (effectiveness is 0 in the meanwhile)
exposure_vaccination_delay = 0                                       # Delay between contact at-risk and vaccine administration (applicable only for Contactbased_link_removal_general)
T_data = 180                                                         # Duration of the input temporal networks (days)
n_runs = 50                                                          # Number of stochastic runs per network
n_initial_I = 10                                                     # Number of infected seeds at time 0
degree_vaccination_threshold = 0                                     # Minimum degree to be eligible to receive PrEP vaccination
prep_vaccination = 1                                                 # 1 is true, 0 is false (no PrEP vaccination is given)
msm_population = 65_000                                              # MSM population in the Paris region (not in the networks)
verbose = 1                                                         
save_state = 1                                                       # 0 is false, 1 is true (not save/save state file. State files are heavy.)

# -- applies for "PY_ages_runs_vaccines_CORRECTED" only --
class No_removal_parameters:
    def __init__(self):
        # these lists are zipped and must have the same length
        self.p_detection_list = [.2,.5,.8]                                                       # detection probabilities
        self.beta_q_list = [.42,.32,.28]                                                         # transmission rate beta of Id MSM (see SI)
        self.beta_list = list(np.array(self.beta_q_list) * 1)                                    # transimssion rate beta of I MSM (see SI)
        self.start_simulation_date_list = pd.to_datetime(['5-6-2022','5-8-2022','5-7-2022'])     # start simulation date

        # non zipped 
        self.total_doses_to_be_given_list = [10_000]                                             # Total number of PEP doses to be administrated

    # -- do not edit. These are parameters useful only for behavioral changes. Here some defaults are assigned. --
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
            
            # changeable
            self.save_weights = 1                                                                # 1: save weights file, 0 does not. Weights file is heavy.

    def debugging(self):
        if(len(self.beta_q_list) != len(self.p_detection_list)):
            sys.exit("wrong zipped params length")
        if(len(self.beta_q_list) != len(self.start_simulation_date_list)):
             sys.exit("wrong zipped params length")
        if(len(self.beta_q_list) != len(self.beta_list)):
             sys.exit("wrong zipped params length")

class Removal_parameters:                                                       
    def __init__(self):
        self.total_doses_to_be_given = 10_000                              # Total number of PEP doses to be administrated
        self.settings = [2]                                                # Settings to be launched. Settings are defined below.                                     
        self.back_in_time = 14                                             # see details
        self.save_weights = 1                                              # 0 is false. weights file is very heavy.

    # -- default, to leave as is --
    def default_parameters(self, analysis_type):
        if analysis_type == "Targeted_Mixed_Removal" or analysis_type == "Random_Mixed_Removal":
            self.back_in_time = -100

class Setting_Parameters:
    def __init__(self, setting):
        self.setting = setting
        
    # -- three possible settings are defined here below. User can define its own setting --
    def set_parameters(self, start_simulation_date_default):
        if self.setting == 1:
            self.p_detection = 0.2                                                                                                 # detection probability
            self.beta_q = 0.42                                                                                                     # transmission rate beta of Id MSM (see SI)
            self.beta = self.beta_q * 1                                                                                            # transmission rate beta of I MSM (see SI)
            self.lag = -1                                                                                                          # first date in surveillance data - start simulation date
            self.prem_list = [77]                                                                                                  # probability (%) of averting a contact if having changed behavior
            self.rem_list = [18]                                                                                                   # cumulative total % of MSM that are changing behavior
            self.start_simulation_date = start_simulation_date_default + timedelta(days=self.lag)                                  # start simulation date
            self.start_degree_date = self.start_simulation_date                                                                    # see details
            self.start_rem_date_list = pd.to_datetime(['6-18-2022']) #pd.to_datetime(['6-15-2022','6-18-2022','6-21-2022'])        # starting date of behavioral changes
            self.end_rem_date_list = pd.to_datetime(['7-18-2022']) #pd.to_datetime(['7-15-2022','7-18-2022','7-21-2022'])          # ending date of behavioral changes
            self.end_degree_date_list = self.end_rem_date_list                                                                     # see details
            self.only_non_vaccinated_change_behavior = 1                                                                           # 1: Only non-smallpox vaccinated MSM can change behavior. 0: All MSM can
            
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
