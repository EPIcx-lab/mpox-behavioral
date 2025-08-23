'''
READ OUTPUT FILES
Output files can be read more easily in the form of pandas dataframes by using the Python function: 'Unified_Code_Load_Results_threshold'. This function assigns names to the pandas dataframe columns. 

-- INPUTS --
- folder must be the path to the folder containing the three folders: state, results, and weights. Example: folder=f'./OUTPUTS/precision/{network}/{analysis_type}'
- analysis_code is an integer that must be defined according to analysis_type. It can be done through the function 'define_analysis_code(analysis_type)'.
- also_state=1 if you want to load the state file also, 0 otherwise
- also_weights=1 if you want to load the weights file also, 0 otherwise

Notice that the number of outputs of the function depends on also_state and also_weights
For each file described below, the number corresponds to the column (first column is numbered as 0), and the string to the name assigned to the pandas dataframe column by the function

## RESULTS FILE DESCRIPTION ##
0:'time':                        day
1:'S'                            # of Susceptible MSM
2:'E'                            # of Exposed MSM
3:'IQ'                           # of ID MSM (see SI)
4:'I'                            # of I MSM (see SI)
5:'Q'                            # of isolated MSM
6:'R'                            # of removed MSM (previously I)
7:'RQ'                           # of removed MSM (previously IQ)
8:'All_new_I_cases'              # of new infections that day (both I and IQ)
9:'new_I1_cases'                 # of new cases that day (only I)
10:'behavior'                    # of MSM changing behavior
11:'eligible_contacts'           # of contacts eligible to change behavior (see SI and manuscript, applicable only if analysis_type=='Contactbased_link_removal_general') 
12:'non_eligible_contacts'       # of contacts non eligible to change behavior (see SI and manuscript, applicable only if analysis_type=='Contactbased_link_removal_general') 
13:'averted_eligible_contacts'   # of averted contacts among the eligible ones (applicable only if analysis_type=='Contactbased_link_removal_general') 
14:'run'                         stochastic run number

## STATE FILE DESCRIPTION ##
The state file provides some information about each node in the network. Information always refers to the end of the simulation. For example, the column 'compartment' refers to the compartment of the MSM at the end of the epidemic.

0:'id':                                   ID of the node
1:'compartment'                           compartment number. Coding: 0:S, 1:E, 2:I, 3:ID, 4:Q, 5:R, 6:RD
2:'vacc_status'                           see details on vacc_status
3:'infecting_id'                          ID of the infecting node
4:'generation',                           infection generation (If I am infected by a seed, generation=1)
5:'detected'                              0: non-detected, 1:detected, -1: never infected
6:'inf_time'                              time when becoming infectious (entering either I or ID compartment)             
7:'gen_time'                              generation time
8:'vacc_time'                             time when vaccinated
9:'exposure_time'                         time when entering the E compartment
10:'age'                                  age
11:'behavior'                             1: changed behavior, 0: did not change behavior
12:'change_time'                          time when changing behavior
13:'compartment_when_vaccinated' 
14:'compartment_when_changing_behavior'
15:'run'                                  stochastic run number

-- details on vacc_status --
0: never vaccinated
1: PEP vaccinated (effective)
2: PrEP vaccinated (effective)
3: smallpox vaccinated
-200: PrEP received from vacc_status=0, not effective yet
5: PrEP received from vacc_status=3, not effective yet 
-1:listed to receive PEP, not received yet
-100: PEP received (not effective yet)

## WEIGHTS FILE DESCRIPTION ##
The weights file copies the input temporal network and specifies if each link was removed or not (due to behavioral changes).

0:'run'                             stochastic run number
1:'time',                           day
2:'weight'                          1: link not eliminated, 0: link eliminated

'''

def Unified_Code_Load_Results_threshold(folder, start_simulation_date, beta_q, beta, epsilon, mu_1_before_start, mu_1_May, mu_1_June, p_detection, mu, \
                              VES_pep, VEI_pep, VEE_pep, VES_smallpox, VEI_smallpox, VEE_smallpox,  VES_firstdose, VEI_firstdose, VEE_firstdose, T_data, T_simul, n_runs, n_initial_I, degree_vaccination_threshold, start_day_vaccines, \
                              end_day_vaccines, start_day_firstdose, exposure_vaccination_delay, total_doses_to_be_given, first_age_to_immunize, \
                              efficacy_delay_pep, efficacy_delay_prep, msm_population, prep_vaccination, analysis_code, also_weights, also_state, interrupt_reference_day, prem=None, rem=None,
                              only_non_vaccinated_change_behavior = None, start_removal_time=None, end_removal_time=None, back_in_time=None, start_degree_time=None, end_degree_time=None):
    
    if analysis_code == 1:
        back_in_time = -100
        prem = 0
        rem = 0
        start_degree_time = (pd.to_datetime('6-30-2022') - start_simulation_date).days + 1
        end_degree_time = (pd.to_datetime('7-1-2022') - start_simulation_date).days + 1
        start_removal_time = (pd.to_datetime('6-30-2022')  - start_simulation_date).days + 1
        end_removal_time = (pd.to_datetime('7-1-2022') - start_simulation_date).days + 1
        only_non_vaccinated_change_behavior = -1
        
    elif analysis_code == 2 or analysis_code == 3:
        back_in_time = -100
    elif analysis_code == 4:
        pass
    else:
        sys.exit("ERROR: got invalid analysis code")
        
    filename = 19*['%.3f'] + 23*['%d']                                      
    filename = '_'.join(filename)
    filename = filename%(beta_q, beta, epsilon, mu_1_before_start, mu_1_May, mu_1_June, p_detection, mu, \
                        VES_pep, VEI_pep, VEE_pep, VES_smallpox, VEI_smallpox, VEE_smallpox, VES_firstdose, VEI_firstdose, VEE_firstdose, \
                        prem, rem, T_data, T_simul, n_runs, n_initial_I, degree_vaccination_threshold, \
                        start_day_vaccines, end_day_vaccines, start_day_firstdose, exposure_vaccination_delay, total_doses_to_be_given, \
                        first_age_to_immunize, efficacy_delay_pep, efficacy_delay_prep, start_removal_time, end_removal_time, back_in_time, interrupt_reference_day, \
                        start_degree_time, end_degree_time, msm_population, prep_vaccination, only_non_vaccinated_change_behavior, analysis_code)
    
    # results
    results = pd.read_csv(f'{folder}/results/{filename}.csv', header=None)
    results = results.rename(columns={0:'time',1:'S',2:'E',3:'IQ',4:'I',5:'Q',6:'R',7:'RQ',8:'All_new_I_cases',9:'new_I1_cases',
                                      10:'behavior', 11:'eligible_contacts', 12:'non_eligible_contacts', 13:'averted_eligible_contacts',
                                      14:'run'})
    
    if results.loc[:,['S','E','IQ','I','Q','R','RQ']].sum(axis=1).var() != 0:
        sys.exit('ERROR: Population not conserved')
        
    runs = results['run'].unique()
    results = results[results["time"]>0]
    results['time'] = results['time'].apply(lambda x: x-1)
    results['newItot'] = list(itertools.chain.from_iterable(list(map(lambda run: \
            np.concatenate(([0],-np.diff(results.loc[results['run']==run,['S','E']].sum(axis=1)))),runs))))
    
    #state
    if also_state:
        state = pd.read_csv(f'{folder}/states/{filename}.csv', header=None)
        state = state.rename(columns={0:'id', 1:'compartment', 2:'vacc_status', 3:'infecting_id', 4:'generation',
                                    5:'detected', 6:'inf_time', 7:'gen_time', 8:'vacc_time',9:'exposure_time',
                                    10:'age',11:'behavior',12:'change_time',13:'compartment_when_vaccinated',
                                    14:'compartment_when_changing_behavior', 15:'run'})

        if len(state[(state['vacc_status']==-1) & (state['vacc_time']>0)]) > 0:
            sys.exit("ERROR: Got people waiting for the vaccine with vaccination time >0")
        if len(state[(state['inf_time'] < state['exposure_time']) & (state['vacc_status']==0)]) > 0:
            sys.exit("Got people non-vaccinated with infection time < exposure time")
        if len(state[(state['vacc_time']>0) & (state['vacc_time']-state['exposure_time']<0)]) > 0:
            sys.exit("ERROR: some nodes have exposure_time > vacc_time ")
   
    # weights
    if also_weights:
        weights = pd.read_csv(f'{folder}/weights/{filename}.csv', header=None)
        weights = weights.rename(columns={0:'run', 1:'time', 2:'weight'})
    
    if also_state:
        if also_weights:
            return filename, results, state, weights
        else:
            return filename, results, state, None
    else:
        if also_weights:
            return filename, results, None, weights
        else:
            return filename, results, None, None

# -----------------------------------------------------------------------------------------------------------
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
        sys.exit(f"Got invalid analysis type {analysis_type} \n")
    return analysis_code
