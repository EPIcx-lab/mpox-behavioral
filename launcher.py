from parameters import *

if os.system(f"g++ -std=c++11 mine.cpp functions.cpp -o {exe_file_name}") == 0:
    print ("Script compiled")
else:
    print ("Failed: script not compiled")
    exit()

start = time.perf_counter()
analysis_code = define_analysis_code(analysis_type)

if analysis_type in ['Targeted_Mixed_Removal','Random_Mixed_Removal',"Contactbased_link_removal_general"]:
    # -- removal parameters --------------------------------------------------------------------------
    removal_parameters = Removal_parameters()
    removal_parameters.default_parameters(analysis_type)

    settings = removal_parameters.settings
    total_doses_to_be_given = removal_parameters.total_doses_to_be_given
    back_in_time = removal_parameters.back_in_time
    save_weights = removal_parameters.save_weights

    for count_net, net in enumerate(net_list, start=1):
        for count_setting, setting in enumerate(settings, start=1):
            
             # -- settings related quantities --------------------------------------------------------
            params_from_setting = Setting_Parameters(setting)
            params_from_setting.set_parameters(start_simulation_date_default)
            
            p_detection = params_from_setting.p_detection
            beta = params_from_setting.beta
            beta_q = params_from_setting.beta_q
            prem_list = params_from_setting.prem_list
            rem_list = params_from_setting.rem_list
            start_simulation_date = params_from_setting.start_simulation_date
            start_degree_date = params_from_setting.start_degree_date
            start_rem_date_list = params_from_setting.start_rem_date_list
            end_rem_date_list = params_from_setting.end_rem_date_list
            end_degree_date_list = params_from_setting.end_degree_date_list
            only_non_vaccinated_change_behavior = params_from_setting.only_non_vaccinated_change_behavior

            for count_rem_date, (start_rem_date, end_rem_date, end_degree_date) in enumerate(zip(start_rem_date_list, end_rem_date_list, end_degree_date_list), start=1):
                # -- strings -----------------------------------------------------
                ids_filename, ages_filename, net_filename, indir,  outdir_results, outdir_states, outdir_weights = \
                define_paths(maindir, net_folder, net, subfolder)

                # -- time related quantities ------------------------------------------------------------
                May_days, After_June_days, mu_1_vector, T_simul = Time_Related_Quantities(start_simulation_date, 
                end_simulation_date, mu_1_May, mu_1_June, mu_1_after_June)

                start_day_vaccines, end_day_vaccines, start_day_firstdose, start_day_degree, end_day_degree, interrupt_reference_day,\
                day_start_rem, day_end_rem = list(map(lambda xx: From_date_to_day(xx, start_simulation_date), \
                    [start_vaccines_date, end_vaccines_date, start_firstdose_date, start_degree_date, end_degree_date, interrupt_reference_date,\
                    start_rem_date, end_rem_date]))
            
                # -- debugging ----------------------------------------------------
                if(T_simul > T_data):
                    sys.exit('ERROR: Got T_simul > T_data')
                if(len(mu_1_vector) != T_simul):
                    sys.exit('ERROR: Got len(mu_1_vector) != T_simul')
                if(len(start_rem_date_list) != len(end_rem_date_list)):
                    sys.exit("ERROR: start_rem_date_list and end_rem_date_list have different len")

                for count_prem, prem in enumerate(prem_list, start=1):
                    for count_rem, rem in enumerate(rem_list, start=1):
                        
                        # -- print ------------------------------------------------------
                        print(f"******** Net {count_net}/{len(net_list)}; Setting {count_setting}/{len(settings)}; RemDate {count_rem_date}/{len(start_rem_date_list)};"
                            f"Prem {count_prem}/{len(prem_list)}; Rem {count_rem}/{len(rem_list)} ********* ")

                        filename = 19*['%.3f'] + 23*['%d']                                      
                        filename = '_'.join(filename)
                        filename = filename%(beta_q, beta, epsilon, mu_1_before_start, mu_1_May, mu_1_June, p_detection, mu, \
                                            VES_pep, VEI_pep, VEE_pep, VES_smallpox, VEI_smallpox, VEE_smallpox, VES_firstdose, VEI_firstdose, VEE_firstdose, \
                                            prem, rem, T_data, T_simul, n_runs, n_initial_I, degree_vaccination_threshold, \
                                            start_day_vaccines, end_day_vaccines, start_day_firstdose, exposure_vaccination_delay, total_doses_to_be_given, \
                                            first_age_to_immunize, efficacy_delay_pep, efficacy_delay_prep, day_start_rem, day_end_rem, back_in_time, interrupt_reference_day, \
                                            start_day_degree, end_day_degree, msm_population, prep_vaccination, only_non_vaccinated_change_behavior, analysis_code)

                        numb_list = [f'./{exe_file_name}'] + \
                                    ["%f" % ii for ii in [beta_q, beta, epsilon, p_detection, mu, VES_pep, VEI_pep, VEE_pep, VES_smallpox, VEI_smallpox, VEE_smallpox, VES_firstdose, VEI_firstdose, VEE_firstdose, prem, rem]] + \
                                    ["%d" % ii for ii in [T_data, T_simul, n_runs, n_initial_I, degree_vaccination_threshold, verbose, start_day_vaccines, end_day_vaccines, start_day_firstdose, exposure_vaccination_delay, \
                                                        total_doses_to_be_given, first_age_to_immunize, efficacy_delay_pep, efficacy_delay_prep, day_start_rem, day_end_rem, back_in_time, \
                                                        interrupt_reference_day, start_day_degree, end_day_degree, save_state, save_weights, msm_population, prep_vaccination, only_non_vaccinated_change_behavior, analysis_code]] + \
                                    ["%f" % ii for ii in mu_1_vector] + \
                                    [outdir_results, outdir_states, outdir_weights, indir, net_filename, ids_filename, ages_filename, filename]
                        str_launch = ' '.join(numb_list)

                        # execute ----------------------------------------------------
                        terminal_input = str_launch
                        subprocess.call([terminal_input], shell=True)
                        print(f"Time: {(time.perf_counter()-start)/60} min")

elif analysis_type == "PY_ages_runs_vaccines_CORRECTED":
    # -- removal parameters ------------------------------------------------------
    no_removal_parameters = No_removal_parameters()
    no_removal_parameters.default_parameters(analysis_type)
    no_removal_parameters.debugging()

    total_doses_to_be_given_list = no_removal_parameters.total_doses_to_be_given_list
    p_detection_list =  no_removal_parameters.p_detection_list
    beta_list =  no_removal_parameters.beta_list
    beta_q_list =  no_removal_parameters.beta_q_list
    start_simulation_date_list = no_removal_parameters.start_simulation_date_list

    back_in_time = no_removal_parameters.back_in_time
    prem = no_removal_parameters.prem
    rem = no_removal_parameters.rem
    start_degree_date = no_removal_parameters.start_degree_date
    end_degree_date = no_removal_parameters.end_degree_date
    start_rem_date = no_removal_parameters.start_rem_date
    end_rem_date = no_removal_parameters.end_rem_date
    save_weights = no_removal_parameters.save_weights
    only_non_vaccinated_change_behavior = no_removal_parameters.only_non_vaccinated_change_behavior

    for count_net, net in enumerate(net_list, start=1):
        for count_doses, total_doses_to_be_given in enumerate(total_doses_to_be_given_list, start=1):
            for count, (beta_q, beta, p_detection, start_simulation_date) in enumerate(zip(beta_q_list, beta_list, p_detection_list, start_simulation_date_list), start=1):
                # -- strings -----------------------------------------------------
                ids_filename, ages_filename, net_filename, indir,  outdir_results, outdir_states, outdir_weights = \
                define_paths(maindir, net_folder, net, subfolder)

                # -- time related quantities ------------------------------------------------------------
                May_days, After_June_days, mu_1_vector, T_simul = Time_Related_Quantities(start_simulation_date, \
                end_simulation_date, mu_1_May, mu_1_June, mu_1_after_June)

                start_day_vaccines, end_day_vaccines, start_day_firstdose, start_day_degree, end_day_degree, interrupt_reference_day,\
                day_start_rem, day_end_rem = list(map(lambda xx: From_date_to_day(xx, start_simulation_date), \
                    [start_vaccines_date, end_vaccines_date, start_firstdose_date, start_degree_date, end_degree_date, interrupt_reference_date,\
                    start_rem_date, end_rem_date]))
                
                # -- debugging ----------------------------------------------------
                if(T_simul > T_data):
                    sys.exit('ERROR: Got T_simul > T_data')
                if(len(mu_1_vector) != T_simul):
                    sys.exit(f'ERROR: Got len(mu_1_vector) != T_simul: {len(mu_1_vector)} vs {T_simul}')

                print(f"******** Net {count_net}/{len(net_list)}; Doses {count_doses}/{len(total_doses_to_be_given_list)};"
                                f"Scenario {count}/{len(beta_list)}; ********* ")

                filename = 19*['%.3f'] + 23*['%d']                                      
                filename = '_'.join(filename)
                filename = filename%(beta_q, beta, epsilon, mu_1_before_start, mu_1_May, mu_1_June, p_detection, mu, \
                                    VES_pep, VEI_pep, VEE_pep, VES_smallpox, VEI_smallpox, VEE_smallpox, VES_firstdose, VEI_firstdose, VEE_firstdose, \
                                    prem, rem, T_data, T_simul, n_runs, n_initial_I, degree_vaccination_threshold, \
                                    start_day_vaccines, end_day_vaccines, start_day_firstdose, exposure_vaccination_delay, total_doses_to_be_given, \
                                    first_age_to_immunize, efficacy_delay_pep, efficacy_delay_prep, day_start_rem, day_end_rem, back_in_time, interrupt_reference_day, \
                                    start_day_degree, end_day_degree, msm_population, prep_vaccination, only_non_vaccinated_change_behavior, analysis_code)

                numb_list = [f'./{exe_file_name}'] + \
                            ["%f" % ii for ii in [beta_q, beta, epsilon, p_detection, mu, VES_pep, VEI_pep, VEE_pep, VES_smallpox, VEI_smallpox, VEE_smallpox, VES_firstdose, VEI_firstdose, VEE_firstdose, prem, rem]] + \
                            ["%d" % ii for ii in [T_data, T_simul, n_runs, n_initial_I, degree_vaccination_threshold, verbose, start_day_vaccines, end_day_vaccines, start_day_firstdose, exposure_vaccination_delay, \
                                                total_doses_to_be_given, first_age_to_immunize, efficacy_delay_pep, efficacy_delay_prep, day_start_rem, day_end_rem, back_in_time, \
                                                interrupt_reference_day, start_day_degree, end_day_degree, save_state, save_weights, msm_population, prep_vaccination, only_non_vaccinated_change_behavior, analysis_code]] + \
                            ["%f" % ii for ii in mu_1_vector] + \
                            [outdir_results, outdir_states, outdir_weights, indir, net_filename, ids_filename, ages_filename, filename]
                str_launch = ' '.join(numb_list)

                # execute ----------------------------------------------------
                terminal_input = str_launch
                subprocess.call([terminal_input], shell=True)
                print(f"Time: {(time.perf_counter()-start)/60} min")
else:
    sys.exit("Got invalid analysis type")

print(f"Total time: {(time.perf_counter()-start)/60} min")


# compile ----------------------------------------------------
# if os.system(f"g++ -std=c++17 mine.cpp functions.cpp -o {exe_file_name}") == 0:
#     print ("Script compiled")
# else:
#     print ("Failed: script not compiled")
#     exit()