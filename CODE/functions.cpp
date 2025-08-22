// Libraries *********************************************************************************
#include "functions.h"
#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <algorithm>
#include <vector>
#include <iterator>
#include <ctime>
#include <fstream>
#include <numeric>
#include <utility>
#include <tuple>

// *********************************************************************************************************
std::ofstream output;
std::ifstream input;
std::random_device rd;

// ############################################################################################################
// ############################################################################################################

void Non_Negative_Debugging(int var, std::string message)
{
    if (var < 0)
    {
        std::cout << "Error: got negative variable in: " << message << '\n';
        exit(EXIT_FAILURE);
    }
}
void Rates_Debugging(float rate, std::string rate_name)
{
    if (rate < 0)
    {
        std::cout << "ERROR: Negative " << rate_name << '\n';
        exit(EXIT_FAILURE); // wonderful command to interrupt the script
    }
    else if (rate == 0)
    {
        std::cout << "WARNING: " << rate_name << " is = 0 \n";
    }
}

void Strong_Rates_Debugging(float rate, std::string rate_name)
{
    if (rate < 0)
    {
        std::cout << "ERROR: Negative " << rate_name << '\n';
        exit(EXIT_FAILURE); // wonderful command to interrupt the script
    }
    if (rate > 1)
    {
        std::cout << "ERROR: Got larger than 1 " << rate_name << '\n';
        exit(EXIT_FAILURE); // wonderful command to interrupt the script
    }
}

void Weak_Rates_Debugging(float rate, std::string rate_name)
{
    if (rate < 0)
    {
        std::cout << "ERROR: Negative " << rate_name << '\n';
        exit(EXIT_FAILURE);
    }
}

void MayorEqual_Debugging(float arg1, float arg2, std::string arg1_name, std::string arg2_name){
    if(arg1 >= arg2){
        std::cout<<"ERROR: got " << arg1_name << " >= " << arg2_name << '\n';
        exit(EXIT_FAILURE);
    }
}

void Mayor_Debugging(float arg1, float arg2, std::string arg1_name, std::string arg2_name){
    if(arg1 > arg2){
        std::cout<<"ERROR: got " << arg1_name << " > " << arg2_name << '\n';
        exit(EXIT_FAILURE);
    }
}

void Solve_Vaccine_Queue(int current_time, STATE &state, std::vector<std::vector<int> > &Queue_Vaccines){
    for(auto node : Queue_Vaccines[current_time-1]){
        if(state.vacc_status[node] == -1){
            if(state.vacc_time[node] == 0){
                state.vacc_status[node] = state.vacc_status[node]*100;
                state.vacc_time[node] = current_time;
            }else{
                std::cout<<"ERROR: Got node with non 0 vacc_time in the Solve_Vaccine_Queue function \n";
                exit(EXIT_FAILURE);
            } 
        }else{
            std::cout<<"ERROR: Got node with >= 0 vacc_status in the Solve_Vaccine_Queue function \n";
            exit(EXIT_FAILURE);
        }
    }       
}

void Activate_PEP_doses(int current_time, unsigned int efficacy_delay_pep, STATE &state){
    for(auto id : state.ids){
        if(((state.vacc_time[id]+efficacy_delay_pep) == current_time) && (state.vacc_status[id]==-100)){
            state.vacc_status[id] = 1;
        }
        if(((state.vacc_time[id]+efficacy_delay_pep) < current_time) && (state.vacc_status[id]==-100)){
            std::cout<<"sum: " << state.vacc_time[id]+efficacy_delay_pep << '\n';
            std::cout<<"time: " << current_time << " vacc_time: " << state.vacc_time[id] << " vacc status: " << state.vacc_status[id] << 
            " efficacy delay pep: " << efficacy_delay_pep <<'\n';
            std::cout<<"ERROR: agent with status -100 should already have vacc_status=1";
            exit(EXIT_FAILURE);
        }
    }
}

void Activate_prep_doses(int current_time, unsigned int efficacy_delay_prep, STATE &state){
    for(auto id : state.ids){
        if(((state.vacc_time[id]+efficacy_delay_prep) == current_time) && ((state.vacc_status[id]==-200) || (state.vacc_status[id]==5))){
            state.vacc_status[id] = 2;
        }
        if(((state.vacc_time[id]+efficacy_delay_prep) < current_time) && ((state.vacc_status[id]==-200) || (state.vacc_status[id]==5))){
            std::cout<<"sum: " << state.vacc_time[id]+efficacy_delay_prep << '\n';
            std::cout<<"time: " << current_time << " vacc_time: " << state.vacc_time[id] << " vacc status: " << state.vacc_status[id] << 
            " efficacy delay pep: " << efficacy_delay_prep <<'\n';
            std::cout<<"ERROR: agent with status -200 or 5 should already have vacc_status=2";
            exit(EXIT_FAILURE);
        }
    }
}

unsigned int Extract_number_nodes_to_vaccines_from_custom_distribution()
{
    float p0 = 30.0/64.0;   //30
    float p1 = 31.0/64.0;   //31
    float p2 = 2.0/64.0;    //2
    float p3 = 0.0;         //0
    float p4 = 1.0/64.0;    //1

    std::discrete_distribution<unsigned int> custom_distribution( {p0, p1, p2, p3, p4} );
    unsigned int number_nodes_to_vaccine = custom_distribution(rd);
    return number_nodes_to_vaccine;
}

unsigned int PEP_Vaccination(std::vector<std::vector<CONTACT> > &temp_network, STATE& state, std::vector<std::vector<int> > &Queue_Vaccines, int current_time,
                    unsigned int infected_node, unsigned int number_nodes_to_vaccine, int exposure_vaccination_delay, bool debugging){

    // debugging
    Weak_Rates_Debugging(exposure_vaccination_delay, "exposure vaccination delay");
    Weak_Rates_Debugging(number_nodes_to_vaccine, "number nodes to vaccine");
    
    // no vaccines case
    if(number_nodes_to_vaccine==0){
        return 0;
    }
   
    // definition of start and end time (of vaccination)
    int start_search_time = state.inf_time[infected_node] + 1;   //it starts infecting new people starting from the day following the infection day!
    int end_search_time = current_time;
    if(start_search_time > end_search_time){
        std::cout<<"ERROR: Got start_search_time > end_search_time in reactive_vaccination function \n";
        exit(EXIT_FAILURE);
    }
     
    // Prepare the candidates for the vaccination
    int total_eligible = 0;
    std::vector<ID_AND_EXPOSURE_TIME> Candidates_to_Vaccine;
    for(auto tt=start_search_time; tt<=end_search_time; tt++){

        // debugging. retrieve contact list at tt-1
        if(tt<1){
            std::cout<<"ERROR in the Reactive Vaccination function got tt<1 \n";
            exit(EXIT_FAILURE);
        }
        
        std::vector<CONTACT> vacc_contactList = temp_network[tt - 1];

        if(vacc_contactList.size()==0){
            std::cout<<"ERROR: Got vacc_contactList of size 0. Probably you are checking the list from time 0 and not 1 \n";
            exit(EXIT_FAILURE);
        }
         
        // iterate on the contact list. Non-vaccinated old contacts are stored in a list
        for(auto contact_iterator = vacc_contactList.begin(); contact_iterator != vacc_contactList.end(); contact_iterator++){
            unsigned int i = (*contact_iterator).i;
            unsigned int j = (*contact_iterator).j;
            float w = (*contact_iterator).w;
            
            if((i == infected_node) && (state.vacc_status[j]==0) && ((state.compartment[j]==0) || (state.compartment[j]==1) || 
                (state.compartment[j]==3) || (state.compartment[j]==2) || (state.compartment[j]==5))){
                ID_AND_EXPOSURE_TIME id_and_exposure_time;
                id_and_exposure_time.id = j;
                id_and_exposure_time.exposure_time = tt;
                Candidates_to_Vaccine.push_back(id_and_exposure_time);
            }
            if((j == infected_node) && (state.vacc_status[i]==0 ) && ((state.compartment[i]==0) || (state.compartment[i]==1) || 
                (state.compartment[i]==3) || (state.compartment[j]==2) || (state.compartment[j]==5))){
                ID_AND_EXPOSURE_TIME id_and_exposure_time;
                id_and_exposure_time.id = i;
                id_and_exposure_time.exposure_time = tt;
                Candidates_to_Vaccine.push_back(id_and_exposure_time);
            }
            if((i == infected_node) || (j == infected_node)){
                total_eligible += 1;
            }
        }
    }

    // std::cout<< "Total eligible:" << total_eligible << "    Len: " << Candidates_to_Vaccine.size() << "  time: " << end_search_time-start_search_time << '\n';
    // Actually vaccine them
    std::vector<int> vaccinated_nodes_this_loop;
    //std::random_device rd;
    std::shuffle(Candidates_to_Vaccine.begin(), Candidates_to_Vaccine.end(), rd);
     
    for(auto vv : Candidates_to_Vaccine){
        
        // If: the node must have not been already vaccinated in this loop
        if (std::find(vaccinated_nodes_this_loop.begin(), vaccinated_nodes_this_loop.end(), vv.id) == vaccinated_nodes_this_loop.end()){
            
            vaccinated_nodes_this_loop.push_back(vv.id);
            
            //If: the node will be either vaccinated now or at exposure_time + exposure_vaccination_delay
            if((vv.exposure_time + exposure_vaccination_delay) > current_time){
                
                state.vacc_status[vv.id] = -1;
                state.exposure_time[vv.id] = vv.exposure_time;

                // If: obv, vaccine can't be put in the queue if vv.exposure_time+exposure_vaccination_delay > T_simul ! (recall: Queue_Vaccines.size)
                if((vv.exposure_time + exposure_vaccination_delay) < Queue_Vaccines.size()){
                    Queue_Vaccines[vv.exposure_time + exposure_vaccination_delay].push_back(vv.id);
                }
                
            }else{
                
                state.vacc_status[vv.id] = -100;
                state.vacc_time[vv.id] = current_time;
                state.exposure_time[vv.id] = vv.exposure_time;
            }
            state.compartment_when_vaccinated = state.compartment;   //niente queue, ecc. Copio subito il compart.
        }
        // If: exit the loop if we have vaccinated number_nodes_to_vaccine nodes. 
        if(vaccinated_nodes_this_loop.size() == number_nodes_to_vaccine){
            break;
        }
    }
    // debugging
    if(debugging){
        std::cout<<"Candidates to vaccine and actually vaccinated: " << Candidates_to_Vaccine.size() << "  " << vaccinated_nodes_this_loop.size() << '\n';
    }

    unsigned int newly_vaccinated = vaccinated_nodes_this_loop.size();
    return newly_vaccinated;
    
    //Candidates_to_Vaccine.clear();
    //vaccinated_nodes_this_loop.clear();
}

void Contact_Based_Behavioral_Changes(std::vector<std::vector<CONTACT> > &temp_network, STATE& state, int& eligible_contacts,
                                     int& non_eligible_contacts, int& averted_eligible_contacts,
                                     int back_in_time, int current_time, unsigned int detected_node, float ratio_how_many_change_behavior, int only_non_vaccinated_change_behavior){

    bool eligibility = false;
    std::uniform_real_distribution<float> unif(0.0, 1.0);
    std::vector<int> compartments_susceptible_to_change_behavior = {0,1,2,3,5}; //S,E,I,IQ,R
    int start_search_time;
    int end_search_time = current_time;
    if (back_in_time < 0){
        start_search_time = 1;   //convention: if got negative back in time, start searching from the beginning of the network (time 1)
    }else{
        if((current_time - back_in_time) >= 1){
             start_search_time = current_time - back_in_time;
        }else{
            start_search_time = 1;
        }
    }
    if(start_search_time > end_search_time){
        std::cout<<"ERROR: Got start_search_time > end_search_time in reactive_vaccination function \n";
        exit(EXIT_FAILURE);
    }
   
    for(auto tt=start_search_time; tt<=end_search_time; tt++){
        // debugging. retrieve contact list at tt-1
        if(tt<1){
            std::cout<<"ERROR in the Reactive Vaccination function got tt<1 \n";
            exit(EXIT_FAILURE);
        }
        
        std::vector<CONTACT> behav_contactList = temp_network[tt - 1];

        if(behav_contactList.size()==0){
            std::cout<<"ERROR: Got behav_contactList of size 0. Probably you are checking the list from time 0 and not 1 \n";
            exit(EXIT_FAILURE);
        }
         
        // iterate on the contact list. Non-vaccinated old contacts are stored in a list
        for(auto contact_iterator = behav_contactList.begin(); contact_iterator != behav_contactList.end(); contact_iterator++){
            unsigned int i = (*contact_iterator).i;
            unsigned int j = (*contact_iterator).j;

            // found variable is true if node j is in a compartment susceptible of behavioral changes
            bool found_j = (std::find(compartments_susceptible_to_change_behavior.begin(), 
                    compartments_susceptible_to_change_behavior.end(), state.compartment[j]) != compartments_susceptible_to_change_behavior.end());
            bool found_i = (std::find(compartments_susceptible_to_change_behavior.begin(), 
                    compartments_susceptible_to_change_behavior.end(), state.compartment[i]) != compartments_susceptible_to_change_behavior.end());

            if((i == detected_node) && (found_j)){
                if(only_non_vaccinated_change_behavior){
                    eligibility = (state.vacc_status[j] == 0);
                }else{
                    eligibility = ((state.vacc_status[j] == 0) || (state.vacc_status[j] == 3));
                }
                if(eligibility){
                    eligible_contacts += 1;
                    float ww = unif(rd);
                    if(ww < ratio_how_many_change_behavior){
                        state.behavior[j] = 1; 
                        state.change_time[j] = current_time;
                        state.compartment_when_changing_behavior[j] = state.compartment[j];
                        averted_eligible_contacts += 1;
                    }else{
                       ;
                    }  
                }else{
                    non_eligible_contacts += 1;
                }
            }
            if((j == detected_node) && (found_i)){
                if(only_non_vaccinated_change_behavior){
                    eligibility = (state.vacc_status[i] == 0);
                }else{
                    eligibility = ((state.vacc_status[i] == 0) || (state.vacc_status[i] == 3));
                }
                if(eligibility){
                    eligible_contacts += 1;
                    float ww = unif(rd);
                    if(ww < ratio_how_many_change_behavior){
                        state.behavior[i] = 1; 
                        state.change_time[i] = current_time;
                        state.compartment_when_changing_behavior[i] = state.compartment[i];
                        averted_eligible_contacts += 1;
                    }else{
                       ;
                    }  
                }else{
                    non_eligible_contacts += 1;
                }
            }
            // old protocol
            // if((i == detected_node) && (found_j)){
            //     if(state.vacc_status[j] == 0){
            //         non_eligible_contacts += 1;
            //         float ww = unif(rd);
            //         if(ww < ratio_how_many_change_behavior){
            //             state.behavior[j] = 1; 
            //             state.change_time[j] = current_time;
            //             state.compartment_when_changing_behavior[j] = state.compartment[j];
            //             averted_non_eligible_contacts += 1;
            //         }else{
            //            ;
            //         }  
            //     }else{
            //         eligible_contacts += 1;
            //     }
            // }
        }
    }
}

void Immunize_From_Age(unsigned int first_age_to_immunize, STATE &state){
    for(int jj = 0; jj < state.compartment.size(); jj++){
        if(state.age[jj] >= first_age_to_immunize){
            state.compartment[jj] = 5;
        }
    }
}

std::vector<int> Get_degree_vector(std::vector<std::vector<CONTACT>> &temp_network, unsigned int start_day_degree, unsigned int end_day_degree, unsigned int size){
   
    std::vector<int> inf_degree(size);                       // store the degree of each node
    
    for (auto tt = start_day_degree; tt <= end_day_degree; tt++)
    {
        std::vector<CONTACT> contactList = temp_network[tt - 1]; // t-1!!! but dependent on the dataset // loop over time. contact list at time tt

        // get degrees for all nodes
        for (auto contact_iterator = contactList.begin(); contact_iterator != contactList.end(); contact_iterator++)
        { // iterate on the contact list
            unsigned int i = (*contact_iterator).i;
            unsigned int j = (*contact_iterator).j;
            float w = (*contact_iterator).w;
            inf_degree[i]++;
            inf_degree[j]++;
        }
    } // end of the first run to compute all degrees
    return inf_degree;
}

void Runtime_targeted_Behavioral_Changes(STATE& state, std::vector<int>& inf_degree, int nrem, int tt, int only_non_vaccinated_change_behavior){
    Weak_Rates_Debugging(nrem, "nrem");
    Rates_Debugging(inf_degree.size(), "shuffled ids size");
    //Rates_Debugging(std::accumulate(inf_degree.begin(), inf_degree.end(),0), "all inf degree are 0");
    if(nrem ==0){
        return;
    }else{
        std::random_device rd;
        std::mt19937 gen(rd());
        
        for(int removed_nodes=1; removed_nodes<=nrem; removed_nodes++){ //se parto da zero ne rimuove nrem+1

            if(std::accumulate(inf_degree.begin(), inf_degree.end(), 0) == 0){
                //std::cout<<"WARNING: all the nodes removed!";
                return;
            }
            //This line must be INTO the loop, the distribution must be updated every time with the updated inf_degree
            std::discrete_distribution<unsigned int> d(inf_degree.begin(), inf_degree.end());
            unsigned int node = d(gen);
            
            if(only_non_vaccinated_change_behavior){
                if(((state.compartment[node] == 0) || (state.compartment[node] == 1) || (state.compartment[node] == 2)
                 || (state.compartment[node] == 3) || (state.compartment[node] == 5)) && (state.vacc_status[node]==0)){
                inf_degree[node] = 0;
                state.change_time[node] = tt;
                state.behavior[node] = 1;
                state.compartment_when_changing_behavior[node] = state.compartment[node];
                }else{
                    std::cout<<"got unexpected error in behav changes \n";
                }
            }else{
                if(((state.compartment[node] == 0) || (state.compartment[node] == 1) || (state.compartment[node] == 2)
                 || (state.compartment[node] == 3) || (state.compartment[node] == 5)) && ((state.vacc_status[node]==0) || 
                    (state.vacc_status[node]==3))){
                inf_degree[node] = 0;
                state.change_time[node] = tt;
                state.behavior[node] = 1;
                state.compartment_when_changing_behavior[node] = state.compartment[node];
                }else{
                    std::cout<<"got unexpected error in behav changes \n";
                }
            }
        }
    }
}

void give_prep_first_dose(STATE& state, std::vector<int>& weekly_first_doses, std::vector<int>& firstdose_unnormalized_probabilities, unsigned int start_day_firstdose, unsigned int msm_population, int tt){
    int doses_given_that_week_index = floor((tt - start_day_firstdose)/7);
    int doses_given_that_day = weekly_first_doses[doses_given_that_week_index]/7*state.ids.size()/msm_population;
   
    if(doses_given_that_day > 0){
        std::random_device rd;
        std::mt19937 gen(rd());
        
        for(int nn=1; nn<=doses_given_that_day; nn++){
            if(std::accumulate(firstdose_unnormalized_probabilities.begin(), firstdose_unnormalized_probabilities.end(), 0) == 0){
                std::cout<<"WARNING: all the nodes vaccinated!";
                return;
            }
            //This line must be INTO the loop, the distribution must be updated every time with the updated inf_degree
            std::discrete_distribution<unsigned int> d(firstdose_unnormalized_probabilities.begin(), firstdose_unnormalized_probabilities.end());
            unsigned int node = d(gen);
            if(state.vacc_status[node]==0){
                state.vacc_status[node] = -200;
                state.vacc_time[node] = tt;
                firstdose_unnormalized_probabilities[node] = 0;
            }else if(state.vacc_status[node]==3){
                state.vacc_status[node] = 5;
                state.vacc_time[node] = tt;
                firstdose_unnormalized_probabilities[node] = 0;
            }
            else{
                ;//std::cout<<"got unexpected error in prep first dose vaccination \n";
            }
        }
        return; 
    }else{
        return;
    }  
}

void Runtime_random_Behavioral_Changes(STATE& state, std::vector<int>& inf_degree_FITTIZIO, int nrem, int tt, int only_non_vaccinated_change_behavior){
    Weak_Rates_Debugging(nrem, "nrem");
    Rates_Debugging(inf_degree_FITTIZIO.size(), "shuffled ids size");
    //Rates_Debugging(std::accumulate(inf_degree_FITTIZIO.begin(), inf_degree_FITTIZIO.end(),0), "all inf degree are 0");
    if(nrem ==0){
        return;
    }else{
        std::random_device rd;
        std::mt19937 gen(rd());
        
        for(int removed_nodes=1; removed_nodes<=nrem; removed_nodes++){ //se parto da zero ne rimuove nrem+1

            if(std::accumulate(inf_degree_FITTIZIO.begin(), inf_degree_FITTIZIO.end(), 0) == 0){
                //std::cout<<"WARNING: all the nodes removed!";
                return;
            }
            //This line must be INTO the loop, the distribution must be updated every time with the updated inf_degree
            std::discrete_distribution<unsigned int> d(inf_degree_FITTIZIO.begin(), inf_degree_FITTIZIO.end());
            unsigned int node = d(gen);
            
            if(only_non_vaccinated_change_behavior){
                if(((state.compartment[node] == 0) || (state.compartment[node] == 1) || (state.compartment[node] == 2)
                 || (state.compartment[node] == 3) || (state.compartment[node] == 5)) && (state.vacc_status[node]==0)){
                inf_degree_FITTIZIO[node] = 0;
                state.change_time[node] = tt;
                state.behavior[node] = 1;
                state.compartment_when_changing_behavior[node] = state.compartment[node];
                }else{
                    std::cout<<"got unexpected error in behav changes \n";
                }
            }else{
                if(((state.compartment[node] == 0) || (state.compartment[node] == 1) || (state.compartment[node] == 2)
                 || (state.compartment[node] == 3) || (state.compartment[node] == 5)) && ((state.vacc_status[node]==0) || 
                    (state.vacc_status[node]==3))){
                inf_degree_FITTIZIO[node] = 0;
                state.change_time[node] = tt;
                state.behavior[node] = 1;
                state.compartment_when_changing_behavior[node] = state.compartment[node];
                }else{
                    std::cout<<"got unexpected error in behav changes \n";
                }
            }
        }
    }
}

void Update_inf_degree(STATE& state, std::vector<int>& inf_degree, int only_non_vaccinated_change_behavior){
    for(auto node : state.ids){
        if(!((state.compartment[node]==0) || (state.compartment[node]==1) || (state.compartment[node]==2) ||
        (state.compartment[node]==3) || (state.compartment[node]==5))){
            inf_degree[node] = 0.0;
        }
        if(only_non_vaccinated_change_behavior){
            if(state.vacc_status[node] !=0 ){
                inf_degree[node] = 0.0;
            }
        }else{
             if((state.vacc_status[node] !=0 ) & (state.vacc_status[node] !=3 )){
                inf_degree[node] = 0.0;
            }
        }    
    }
}

void Update_inf_degree_random(STATE& state, std::vector<int> &inf_degree_FITTIZIO, int only_non_vaccinated_change_behavior){
    for(auto node : state.ids){
        if(!((state.compartment[node]==0) || (state.compartment[node]==1) || (state.compartment[node]==2) ||
        (state.compartment[node]==3) || (state.compartment[node]==5))){
            inf_degree_FITTIZIO[node] = 0.0;
        }
        if(only_non_vaccinated_change_behavior){
            if(state.vacc_status[node] !=0 ){
                inf_degree_FITTIZIO[node] = 0.0;
            }
        }else{
             if((state.vacc_status[node] !=0 ) & (state.vacc_status[node] !=3 )){
                inf_degree_FITTIZIO[node] = 0.0;
            }
        }    
    }
}

void Update_firstdose_unnormalized_probabilities(STATE& state, std::vector<int>& firstdose_unnormalized_probabilities){
    for(auto node : state.ids){
        if((state.vacc_status[node] != 0) && (state.vacc_status[node] != 3)){
            firstdose_unnormalized_probabilities[node] = 0;
        }
        if(!((state.compartment[node]==0) || (state.compartment[node]==1) || (state.compartment[node]==2) ||
        (state.compartment[node]==3) || (state.compartment[node]==5))){
            firstdose_unnormalized_probabilities[node] = 0;
        }
    }
    //int acc = std::accumulate(firstdose_unnormalized_probabilities.begin(), firstdose_unnormalized_probabilities.end(), 0);
    //std::cout<< acc << '\n'; 
}

EPIDEMIC_RESULT Run_The_Epidemic(std::vector<std::vector<CONTACT> > &temp_network, STATE &state, float beta_q, float beta, float sigma, std::vector<float> delta_vector,
                                               float p_detection, float mu, float VES_pep, float VEI_pep, float VEE_pep, float VES_smallpox, float VEI_smallpox, float VEE_smallpox, 
                                               float VES_firstdose, float VEI_firstdose, float VEE_firstdose, int T_simul, int T_data, std::vector<unsigned int> n_initial_vec, int degree_vaccination_threshold, 
                                               unsigned int verbose, unsigned int run, 
                                               int exposure_vaccination_delay, unsigned int start_day_vaccines, unsigned int end_day_vaccines, unsigned int start_day_firstdose,
                                               int total_doses_to_be_given, unsigned int efficacy_delay_pep, unsigned int efficacy_delay_prep, unsigned int day_start_rem, unsigned int day_end_rem, 
                                               float prem, float rem, int back_in_time, int interrupt_reference_day, unsigned int start_day_degree, unsigned int end_day_degree, unsigned int msm_population,
                                               unsigned int prep_vaccination, int only_non_vaccinated_change_behavior, int analysis_code){
    // debugging 
    Rates_Debugging(beta, "beta");
    Rates_Debugging(beta_q, "beta_q");
    Rates_Debugging(sigma, "sigma");
    for (auto delta : delta_vector)
    {
        Rates_Debugging(delta, "delta");
    }
    Rates_Debugging(mu, "mu");
    Strong_Rates_Debugging(VES_pep, "VES_pep");
    Strong_Rates_Debugging(VEI_pep, "VEI_pep");
    Strong_Rates_Debugging(VEE_pep, "VEE_pep");
    Strong_Rates_Debugging(VES_smallpox, "VES_smallpox");
    Strong_Rates_Debugging(VEI_smallpox, "VEI_smallpox");
    Strong_Rates_Debugging(VEE_smallpox, "VEE_smallpox");
    Strong_Rates_Debugging(VES_firstdose, "VES_firstdose");
    Strong_Rates_Debugging(VEI_firstdose, "VEI_firstdose");
    Strong_Rates_Debugging(VEE_firstdose, "VEE_firstdose");
    Rates_Debugging(T_simul, "T_simul");
    Rates_Debugging(n_initial_vec[0], "n_initial_I1");
    Rates_Debugging(n_initial_vec[1], "n_initial_I2");
    Rates_Debugging(run, "run");
    Weak_Rates_Debugging(p_detection, "p_detection");
    Rates_Debugging(start_day_vaccines, "start day vaccines");
    Rates_Debugging(end_day_vaccines, "end day vaccines");
    Rates_Debugging(start_day_firstdose, "start_day_firstdose");
    MayorEqual_Debugging(start_day_vaccines, end_day_vaccines, "start_day_vaccines", "end_day_vaccines");
    MayorEqual_Debugging(start_day_vaccines, T_simul, "start_day_vaccines", "T simul");
    Mayor_Debugging(end_day_vaccines, T_simul, "end_day_vaccines", "T simul");
    Rates_Debugging(total_doses_to_be_given, "total doses to be given");
    if(analysis_code != 1){
        Weak_Rates_Debugging(prem, "prem");
        Mayor_Debugging(day_start_rem, T_simul, "start_day_removal", "T simul");
        Mayor_Debugging(day_end_rem, T_simul, "end_day_removal", "T simul");
        MayorEqual_Debugging(day_start_rem, day_end_rem, "start_day_removal", "end_day_removal");
        Rates_Debugging(efficacy_delay_pep, "efficacy delay pep");
    }
    Mayor_Debugging(beta, 1.0, "beta", "1");
    Mayor_Debugging(beta_q, 1.0, "beta_q", "1");
    Mayor_Debugging(beta, beta_q, "beta", "beta_q");

    // setup the quantities
    EPIDEMIC_RESULT epidemic_result;
    std::vector<INF_CONTACT> atrisk_contacts;
    std::uniform_real_distribution<float> unif(0.0, 1.0);                    // uniform distribution between 0 and 1
    std::vector<unsigned int> new_exposed;
    std::vector<unsigned int> infecting_nodes;
    INF_CONTACT atrisk_contact;
    std::vector<std::vector<CONTACT> > temp_network_copy = temp_network;       // needed to keep trace of removed links
    std::vector<std::vector<int> > results(T_simul + 1, std::vector<int>(15)); // last result will be at column 367!
    int new_I1_cases = 0;
    int new_cases = 0;
    int given_plus_queue_doses = 0;
    std::vector<std::vector<int> > Queue_Vaccines(T_simul);                   //for every day, it will be filled with the nodes to vaccine that day
    std::vector<float> VES_vector{0.0, VES_pep, VES_firstdose, VES_smallpox, 0.0, VES_smallpox};
    std::vector<float> VEI_vector{0.0, VEI_pep, VEI_firstdose, VEI_smallpox, 0.0, VEI_smallpox};
    std::vector<float> VEE_vector{0.0, VEE_pep, VEE_firstdose, VEE_smallpox, 0.0, VEE_smallpox};
    std::vector<int> inf_degree = Get_degree_vector(temp_network, start_day_degree, end_day_degree, state.ids.size()); 
    std::vector<int> inf_degree_six_months = Get_degree_vector(temp_network, 1, T_data, state.ids.size()); 
    std::vector<int> inf_degree_FITTIZIO(state.ids.size(), 1); 
    std::vector<int> firstdose_unnormalized_probabilities(state.ids.size(), 1); 
    std::vector<int> weekly_first_doses{1190, 2705, 4469, 5461, 6446, 6088, 8356, 6046, 5613, 3361, 0, 2783, 2294};

    //setup for prep multipartner vaccination. Only MSM with >=2 partners are eligible for prep
    if(prep_vaccination==1){
        for(auto node: state.ids){
            if(inf_degree_six_months[node] <= degree_vaccination_threshold){
                firstdose_unnormalized_probabilities[node] = 0;
            }
        }
    }
   
    // Initialize results matrix at time 0
    results[0][0] = 0;                                                                // time
    results[0][1] = state.compartment.size() - n_initial_vec[0] - n_initial_vec[1];   // S
    results[0][2] = 0;                                                                // E
    results[0][3] = n_initial_vec[0];                                                 // I1 (IQ)
    results[0][4] = n_initial_vec[1];                                                 // I2 (I)
    results[0][5] = 0;                                                                // Q
    results[0][6] = 0;                                                                // R
    results[0][7] = 0;                                                                // RQ
    results[0][8] = 0;                                                                // All_new_I_cases
    results[0][9] = 0;                                                                // new_I1_cases
    results[0][10] = 0;                                                               // behavior == 1
    results[0][11] = 0;                                                               // eligible contacts
    results[0][12] = 0;                                                               // non eligible contacts
    results[0][13] = 0;                                                               // averted eligible contacts
    results[0][14] = run;                                                             // run

 
    int nrem = round(rem/100*inf_degree.size()/(day_end_rem - day_start_rem + 1));   //se alla riga sotto usassi < come mettia, dovrei togliere il +1
    

    for (auto tt = 1; tt <= T_simul; tt++){ 

       // setup quantities
       int eligible_contacts = 0;
       int non_eligible_contacts = 0;
       int averted_eligible_contacts = 0;
       
       // delta
        float delta = delta_vector[tt - 1];

        // compute mu_q
        if(delta <= mu){
            std::cout<<"ERROR: got delta (muIQ) <= mu \n";
            exit(EXIT_FAILURE);
        }
        float mu_q = mu*delta/(delta - mu);

        // Update Vaccinal status
        Solve_Vaccine_Queue(tt, state, Queue_Vaccines);
        Activate_PEP_doses(tt, efficacy_delay_pep, state);
        Activate_prep_doses(tt, efficacy_delay_prep, state);

        // prep vaccination
        if((tt >= start_day_firstdose) && prep_vaccination==1){
            Update_firstdose_unnormalized_probabilities(state, firstdose_unnormalized_probabilities);
            give_prep_first_dose(state, weekly_first_doses, firstdose_unnormalized_probabilities, start_day_firstdose, msm_population, tt);
        }
        
        // behavioral changes
        if((tt>=day_start_rem) && (tt<=day_end_rem)){
            if(analysis_code == 2){
                Update_inf_degree(state, inf_degree, only_non_vaccinated_change_behavior);
                Runtime_targeted_Behavioral_Changes(state, inf_degree, nrem, tt, only_non_vaccinated_change_behavior);
            }
            if(analysis_code == 3){
                Update_inf_degree_random(state, inf_degree_FITTIZIO, only_non_vaccinated_change_behavior);
                Runtime_random_Behavioral_Changes(state, inf_degree_FITTIZIO, nrem, tt, only_non_vaccinated_change_behavior);
            }
        }

        // targeted links removal
        // we don't impose a temporal constraint here. Imposed only on when nodes change behavior
        std::random_device randd;
        
        if(std::accumulate(state.behavior.begin(), state.behavior.end(), 0) == 0){
                ;    //do nothing: no nodes are changing behavior yet
        }else{
            for (auto link_iterator = temp_network_copy[tt - 1].begin(); link_iterator != temp_network_copy[tt - 1].end(); link_iterator++){
                unsigned int i = (*link_iterator).i;
                unsigned int j = (*link_iterator).j;
                float w = (*link_iterator).w;
                if((state.behavior[i]==0) & (state.behavior[j]==0)){
                    ;
                }else if ((state.behavior[i]==1) & (state.behavior[j]==1))
                {
                    float uu = unif(rd);
                    if(uu <= (2*prem/100 - pow(prem/100, 2))){
                        (*link_iterator).w = 0.0;    //removed link: we set its weight to 0
                    }else{
                        ;
                    }
                }else{
                    float uu = unif(rd);
                    if(uu <= prem/100){
                        (*link_iterator).w = 0.0;    //removed link: we set its weight to 0
                    }else{
                        ;
                    }
                }  
            }
        }
    
        // Find contacts at risk
        for (auto contact_iterator = temp_network_copy[tt - 1].begin(); contact_iterator != temp_network_copy[tt - 1].end(); contact_iterator++)
        { // iterate on the contact list
            unsigned int i = (*contact_iterator).i;
            unsigned int j = (*contact_iterator).j;
            float w = (*contact_iterator).w;

            if(state.compartment[j] == 0){
                if(state.compartment[i] == 2){   //detected. they transmit more
                    atrisk_contact.ns = j;
                    atrisk_contact.ni = i;
                    atrisk_contact.w = w*beta_q;
                    atrisk_contacts.push_back(atrisk_contact);
                } else if(state.compartment[i] == 3){     // non detected. they transmit less
                    atrisk_contact.ns = j;
                    atrisk_contact.ni = i;
                    atrisk_contact.w = w*beta;
                    atrisk_contacts.push_back(atrisk_contact);
                }
            }
            if(state.compartment[i] == 0){
                if(state.compartment[j] == 2){   //detected. they transmit more
                    atrisk_contact.ns = i;
                    atrisk_contact.ni = j;
                    atrisk_contact.w = w*beta_q;
                    atrisk_contacts.push_back(atrisk_contact);
                } else if(state.compartment[j] == 3){     // non detected. they transmit less
                    atrisk_contact.ns = i;
                    atrisk_contact.ni = j;
                    atrisk_contact.w = w*beta;
                    atrisk_contacts.push_back(atrisk_contact);
                }
            }
        }

        // PREPARE TRANSITIONS S->E. I save in a vector which WILL do this transition
        for (auto inf_contact_iterator = atrisk_contacts.begin(); inf_contact_iterator != atrisk_contacts.end(); inf_contact_iterator++)
        {
            if (state.compartment[(*inf_contact_iterator).ns] == 0)
            { // this condition is needed in order to not count twice a node that gets infected by two infected neighbours, so we have to ckeck if it is still S
                // float vacc_eff = Map_Vaccinal_Status(state.vacc_status[(*inf_contact_iterator).ns], VE0, VE1, VE2);  //old one
                // VACCINE_EFFICACY_VALUES vaccine_efficacy_values = Map_Vaccinal_status_to_struct(state.vacc_status[(*inf_contact_iterator).ns], VES, VEI);
                float veff_s = Map_Vaccinal_Status(state.vacc_status[(*inf_contact_iterator).ns], VES_vector);
                float veff_i = Map_Vaccinal_Status(state.vacc_status[(*inf_contact_iterator).ni], VEI_vector);
                float x = unif(rd);
                float weighted_prob = (*inf_contact_iterator).w * (1.0 - veff_s) * (1.0 - veff_i);
                if (x <= weighted_prob)
                {
                    unsigned int new_exposed_node = (*inf_contact_iterator).ns;
                    unsigned int infecting_node = (*inf_contact_iterator).ni;
                    new_exposed.push_back(new_exposed_node);
                    infecting_nodes.push_back(infecting_node);
                }
            }
        }

        //  -- SPONTANEOUS REACTIONS ---------------------------------------------------
        for (int ii = 0; ii < state.compartment.size(); ii++)
        {
            float uu = unif(rd);
            switch (state.compartment[ii])
            {
            case 0:
                break; // I have to do nothing, now, if the state is S
            case 1:    // E -> I1   or  E -> I2
                if (uu <= sigma)
                {
                    float yy = unif(rd);
                    float veff_e = Map_Vaccinal_Status(state.vacc_status[ii], VEE_vector);
                    if (yy >= ((1-veff_e)*p_detection))     // E -> I2 (non detected) (I)  
                    {
                        state.compartment[ii] = 3;
                        state.detected[ii] = 0;
                    }
                    else if (yy < ((1-veff_e)*p_detection)) // E -> I1 (detected) (IQ) //now I have to be detected AND the vaccine must fail
                    {
                        state.compartment[ii] = 2;
                        state.detected[ii] = 1;
                        new_I1_cases += 1;
                    }
                    state.inf_time[ii] = tt;
                    state.gen_time[ii] = tt - state.inf_time[state.infecting_id[ii]];
                    state.generation[ii] = state.generation[state.infecting_id[ii]] + 1;
                    new_cases += 1;
                }
                break;
            case 2: // I1->Q
                if (uu <= delta)
                {
                    state.compartment[ii] = 4;

                    // PEP to the contacts
                    if((tt>= start_day_vaccines) && (tt<=end_day_vaccines) && (given_plus_queue_doses < total_doses_to_be_given)){ //vaccination happens in this time interval

                        //extract the number of nodes to vaccine
                        unsigned int number_nodes_to_vaccine = Extract_number_nodes_to_vaccines_from_custom_distribution();

                        // every single time that someone is vaccinated, I have to update given_plus_queue_doses
                        unsigned int newly_vaccinated = PEP_Vaccination(temp_network, state, Queue_Vaccines, tt, ii, number_nodes_to_vaccine, exposure_vaccination_delay, false);
                        given_plus_queue_doses += newly_vaccinated;
                    }

                    // contact-based behavioral changes
                    if(analysis_code == 4){
                        if((tt>=day_start_rem) && (tt<=day_end_rem)){
                        Contact_Based_Behavioral_Changes(temp_network, state, eligible_contacts, non_eligible_contacts, averted_eligible_contacts, back_in_time, tt, ii, rem/100, 
                        only_non_vaccinated_change_behavior);
                        }
                    }
                }
                break;
            case 3: // I2 -> R
                if (uu <= mu)
                {
                    state.compartment[ii] = 5;
                }
                break;
            case 4: // Q -> RQ
                if (uu <= mu_q)
                {
                    state.compartment[ii] = 6;
                }
                break;
            case 5: // I have to do nothing if state is R
                break;
            case 6: // I have to do nothing if state is RQ
                break;
            default:
                std::cout << "ERROR: entered in the default case of the spontaneous precesses: \n";
                std::cout << "Got compartment: " << state.compartment[ii] << '\n';
                break;
            }
        }

        // do transitions S->E
        if (new_exposed.size() != infecting_nodes.size())
        {
            std::cout << "ERROR: new_exposed size different from infecting_nodes size \n";
            break;
        }

        for (int ii = 0; ii < new_exposed.size(); ii++)
        {
            state.compartment[new_exposed[ii]] = 1;
            state.infecting_id[new_exposed[ii]] = infecting_nodes[ii];
        }

        // UPDATE RESULTS MATRIX. Done easily by summing up stuff of the states matrix
        for (int ii = 1; ii <= 8; ii++)
        {
            results[tt][ii] = std::count(state.compartment.begin(), state.compartment.end(), ii-1); //pay attention to that ii and to that ii-1
        }

        // Initialize results matrix at time tt. Store current time. Store current run.
        results[tt][0] = tt;
        results[tt][8] = new_cases;
        results[tt][9] = new_I1_cases;
        results[tt][10] = std::count(state.behavior.begin(), state.behavior.end(), 1);
        results[tt][11] = eligible_contacts;                                                      // eligible contacts
        results[tt][12] = non_eligible_contacts;                                                  // non eligible contacts
        results[tt][13] = averted_eligible_contacts;                                              // averted eligible contacts
        results[tt][14] = run;

        //given_plus_queue_doses = results[tt][8] + std::count(state.vacc_status.begin(), state.vacc_status.end(), -1);//si puo' fare piu' generale contando gli zeri e togliendoli

        // Interrupt simulation if infected and exposed are finished
        if (((results[tt][2] + results[tt][3] + results[tt][4]) == 0) && (tt<interrupt_reference_day))
        {
            epidemic_result.result = results;
            epidemic_result.interrupted = true;
            epidemic_result.temp_network_copy = temp_network_copy;
            return epidemic_result;
        }

        // RESET QUANTITIES
        new_I1_cases = 0;
        new_exposed.clear();
        atrisk_contacts.clear();
        infecting_nodes.clear();

        // DEBUGGING
        for (int ii = 0; ii <= 10; ii++)
        {
            std::string message = std::to_string(ii);
            message = "column " + message + "of results \n";
            Non_Negative_Debugging(results[tt][ii], message);
        }
    }

    epidemic_result.result = results;
    epidemic_result.interrupted = false;
    epidemic_result.temp_network_copy = temp_network_copy;
    return epidemic_result;
}

std::vector<unsigned int> Generalized_Assign_seeds_and_Ages(STATE &state, unsigned int &n_initial_I, float p_detection, unsigned int first_age_to_immunize, bool debugging = false)
{
    // first, recover aged people. These WON'T be chosen as SEEDS
    for(int jj = 0; jj < state.compartment.size(); jj++){
        if(state.age[jj] >= first_age_to_immunize){
            state.vacc_status[jj] = 3;
        }
    }

    // now shuffle to choose the seeds
    std::uniform_real_distribution<float> unif(0.0, 1.0); // uniform distribution between 0 and 1
    std::vector<unsigned int> copy_of_ids = state.ids;
    std::shuffle(copy_of_ids.begin(), copy_of_ids.end(), rd);
    
    int assigned = 0;
    int ii = 1;

    while(assigned < n_initial_I){
        unsigned int seed = copy_of_ids[ii]; // choose a random node
        if(state.vacc_status[seed] == 0){    // this way, already vaccinated can't be chosen
            float yy = unif(rd);
            if (yy >= p_detection)
            {
                state.compartment[seed] = 3; // put that node into the I2 compartment (not detected)
                state.detected[seed] = 0;
            }
            else if (yy < p_detection)
            {
                state.compartment[seed] = 2; // put that node into the I1 compartment (detected)
                state.detected[seed] = 1;
            }
            state.inf_time[seed] = 0.0;      // set its inf_time t0 0: he got the infection at time 0!
            assigned += 1;
            ii += 1;
        }else{
            ii += 1;                         // node not chosen to be seed. Let's try with the next one
        }
        if(ii > state.compartment.size()){
            std::cout<<"ERROR: could not assign " << n_initial_I << " seed. Probably due to the fact that there are too many R due to age \n";
            exit(EXIT_FAILURE);
        }
    }

    unsigned int n_initial_I1 = count(state.compartment.begin(),state.compartment.end(),2);
    unsigned int n_initial_I2 = count(state.compartment.begin(),state.compartment.end(),3);
    std::vector<unsigned int> n_initial_vec {n_initial_I1, n_initial_I2};
    return n_initial_vec;
}

std::vector<unsigned int> Assign_seeds_and_Ages(STATE &state, unsigned int &n_initial_I, float p_detection, unsigned int first_age_to_immunize, bool debugging = false)
{
    // first, recover aged people. These WON'T be chosen as SEEDS
    for(int jj = 0; jj < state.compartment.size(); jj++){
        if(state.age[jj] >= first_age_to_immunize){
            state.compartment[jj] = 5;
        }
    }

    // now shuffle to choose the seeds
    std::uniform_real_distribution<float> unif(0.0, 1.0); // uniform distribution between 0 and 1
    std::vector<unsigned int> copy_of_ids = state.ids;
    std::shuffle(copy_of_ids.begin(), copy_of_ids.end(), rd);
    
    int assigned = 0;
    int ii = 1;

    while(assigned < n_initial_I){
        unsigned int seed = copy_of_ids[ii]; // choose a random node
        if(state.compartment[seed] == 0){    // this way, already R can't be chosen
            float yy = unif(rd);
            if (yy >= p_detection)
            {
                state.compartment[seed] = 3; // put that node into the I2 compartment (not detected)
                state.detected[seed] = 0;
            }
            else if (yy < p_detection)
            {
                state.compartment[seed] = 2; // put that node into the I1 compartment (detected)
                state.detected[seed] = 1;
            }
            state.inf_time[seed] = 0.0;      // set its inf_time t0 0: he got the infection at time 0!
            assigned += 1;
            ii += 1;
        }else{
            ii += 1;                         // node not chosen to be seed. Let's try with the next one
        }
        if(ii > state.compartment.size()){
            std::cout<<"ERROR: could not assign " << n_initial_I << " seed. Probably due to the fact that there are too many R due to age \n";
            exit(EXIT_FAILURE);
        }
    }

    unsigned int n_initial_I1 = count(state.compartment.begin(),state.compartment.end(),2);
    unsigned int n_initial_I2 = count(state.compartment.begin(),state.compartment.end(),3);
    std::vector<unsigned int> n_initial_vec{n_initial_I1, n_initial_I2};
    return n_initial_vec;
}

STATE Initialize_State_Structure(std::vector<unsigned int> ids, std::vector<unsigned int> ages, unsigned int run)
{
    unsigned int NN = ids.size();
    STATE state;                      // definisco l'oggetto state come una STRUCT di tipo STATE
    state.compartment.assign(NN, 0);  // everyone is Susceptible at time 0
    state.vacc_status.assign(NN, 0);  // everyone with 0 doses
    state.inf_time.assign(NN, 0.0);   // noone is infected
    state.gen_time.assign(NN, 0.0);   // noone is infected
    state.infecting_id.assign(NN, 0); // noone is infected
    state.ids = ids;
    state.run.assign(NN, run);        // specify the run
    state.generation.assign(NN, 0);
    state.detected.assign(NN, -1);    //-1:never infected, 0:non detected, 1:detected
    state.vacc_time.assign(NN, 0.0);
    state.exposure_time.assign(NN, 0.0);
    state.age = ages;
    state.behavior.assign(NN, 0);
    state.change_time.assign(NN, 0);
    state.compartment_when_vaccinated.assign(NN ,100);
    state.compartment_when_changing_behavior.assign(NN, 100);
    return state;
}

std::vector<unsigned int> load_degrees(std::string inputname)
{
    std::vector<unsigned int> ids;
    unsigned int id;
    std::ifstream infile(inputname);
    while (infile >> id)
    {
        ids.push_back(id);
    }
    if (ids.size() == 0)
    {
        std::cout << "ERROR: length of the degrees vector=0. Probably the id file has a header, that must be removed. \n";
        exit(EXIT_FAILURE);
    }
    return ids;
}

std::vector<unsigned int> load_ages(std::string inputname)
{
    std::vector<unsigned int> ids;
    unsigned int id;
    std::ifstream infile(inputname);
    while (infile >> id)
    {
        ids.push_back(id);
    }
    if (ids.size() == 0)
    {
        std::cout << "ERROR: length of the ages vector=0. Probably the id file has a header, that must be removed. \n";
        exit(EXIT_FAILURE);
    }
    return ids;
}

std::vector<unsigned int> loadID(std::string inputname)
{
    std::vector<unsigned int> ids;
    unsigned int id;
    std::ifstream infile(inputname);
    while (infile >> id)
    {
        ids.push_back(id);
    }
    if (ids.size() == 0)
    {
        std::cout << "ERROR: length of the ids vector=0. Probably the id file has a header, that must be removed. \n";
        exit(EXIT_FAILURE);
    }
    if(ids.back()!=(ids.size()-1)){
        std::cout<<"ERROR: last element of ids does not coincide with ids.size()-1. Probably the cause is bad remapping \n";
        exit(EXIT_FAILURE);
    }
    return ids;
}

std::vector<std::vector<CONTACT>> loadData(std::string inputname, unsigned int T_data)
{
    // Define list of contact lists and list of nodes, etc.
    std::string line;
    unsigned int node_i, node_j;
    unsigned int day;
    float weight;
    CONTACT contact;
    std::vector<CONTACT> contact_list;
    std::vector<std::vector<CONTACT> > temp_network;

    // Read first line of inputfile as list of characters and get t,i,j:
    // getline skips first line of the text file
    input.open(inputname);
    if(input.fail()){
        std::cout<<"ERROR: FILE " << inputname << " DOES NOT EXIST \n";
        exit(EXIT_FAILURE);
    }
    getline(input, line);
    input >> node_i >> node_j >> day >> weight;

    // Loop over tt and create list of contact lists:
    for (unsigned int tt = 1; tt <= T_data; tt++)
    {
        while (day == tt && !input.eof())
        {
            contact.i = node_i;
            contact.j = node_j; // feed contact with the info just read
            contact.w = weight;
            contact_list.push_back(contact); // add contact to the contact list

            // Read line and get t,i,j:
            getline(input, line);
            input >> node_i >> node_j >> day >> weight;
        }
        temp_network.push_back(contact_list); // add the contact list (one newtork) to the list of temporal networks
        contact_list.clear();                 // re-initialize contact list
    }
    input.close();
    return temp_network;
}

//not the most genial function of the world: returns 0 if the guy is waiting vaccine or is not vaccined, VE otherwise
/*float Map_Vaccinal_Status_to_VES_or_VEI(int vacc_status, float VE){
    if(vacc_status < 0){
        return 0.0;
    }else if(vacc_status==0){
        return 0.0;
    }else if(vacc_status==1){
        return VE;
    }else{
        std::cout << "ERROR: Got invalid vacc_status \n";
        exit(EXIT_FAILURE);
    }
}*/

/*VACCINE_EFFICACY_VALUES Map_Vaccinal_status_to_struct(int vacc_status, float veffS, float veffI){
    VACCINE_EFFICACY_VALUES vaccine_efficay_values;
    //if vacc_status<0, it means that that node is waiting for being vaccinated. Therefore, now it has vacc_eff =0.
    if(vacc_status < 0){
        vaccine_efficay_values.VES = 0.0;
        vaccine_efficay_values.VEI = 0.0;
    }else{
        switch(vacc_status){
            case 0:
                vaccine_efficay_values.VES = 0.0;
                vaccine_efficay_values.VEI = 0.0;
                break;
            case 1:
                vaccine_efficay_values.VES = veffS;
                vaccine_efficay_values.VEI = veffI;
                break;
            default:
                std::cout << "ERROR: Got invalid vacc_status \n";
                exit(EXIT_FAILURE);
                break;
        }
    }
    return vaccine_efficay_values;    
}*/

/* EXPLANATION 
A negative vaccinal status means that that node has no vaccine because it is waiting either to receive a dose,
or it has received it already and it is waiting that it becomes effective.
Positive vaccine status means that the nodes is vaccined. I do not remember why, but I decided that
vacc_status=1 means 1 dose of a 3rd gen vaccine, vacc_status=2 means 2 doses of a 3rd gen vaccine (not used yet),
vacc_status=3 means 1st gen vaccine, and vacc_status=0 no vaccination.
*/

float Map_Vaccinal_Status(int vacc_status, std::vector<float> VE_vector){
    if(vacc_status < 0){
        return 0.0;
    }else{
        return VE_vector[vacc_status];
    }
}

void Save_Results_Matrix(std::string filename, std::vector<std::vector<int> > results, unsigned int run){
    if (run == 1)                               //create new file if run=1, otherwise append
    {
        output.open(filename, std::ios::out);
    }
    else
    {
        output.open(filename, std::ios::out | std::ios::app);
    }

    // first iteration (time 0) (ii = 0)
    for (int jj = 0; jj < (results[0].size() - 1); jj++)
    {
        output << results[0][jj] << ',';
    }
    output << results[0][results[0].size() - 1] << '\n'; // at the last row we want \n and not ,

    // and from now if t=0 it means that the simulation finished before
    for (int ii = 1; ii < results.size(); ii++)
    { // START FROM ii=1 HERE!
        if (results[ii][0] != 0)
        {
            for (int jj = 0; jj < (results[0].size() - 1); jj++)
            {
                output << results[ii][jj] << ',';
            }
            output << results[ii][results[0].size() - 1] << '\n'; // at the last row we want \n and not ,
        }
        else
            break;
    }
    output.close();
}

void Save_State_Matrix(std::string filename, STATE &state, unsigned int run)
{
    if (run == 1)
    {
        output.open(filename, std::ios::out);
    }
    else
    {
        output.open(filename, std::ios::out | std::ios::app);
    }
    for (int jj = 0; jj < state.compartment.size(); jj++)
    {
        output << state.ids[jj] << ',' << state.compartment[jj] << ',' << state.vacc_status[jj] << ',' << state.infecting_id[jj] << ','
               << state.generation[jj] << ',' << state.detected[jj] << ',' << state.inf_time[jj] << ',' << state.gen_time[jj] << ','
               << state.vacc_time[jj] << ',' << state.exposure_time[jj] << ',' << state.age[jj] << ',' << state.behavior[jj] << ','
               << state.change_time[jj] << ',' << state.compartment_when_vaccinated[jj] << ',' << state.compartment_when_changing_behavior[jj] << ','
               << state.run[jj] << '\n';
    }
    output.close();
}

void Save_Weights_Vec(std::string filename, std::vector<std::vector<CONTACT> > &temp_network_copy,  unsigned int T_simul,
                     unsigned int run)
{
    if (run == 1)
    {
        output.open(filename, std::ios::out);
    }
    else
    {
        output.open(filename, std::ios::out | std::ios::app);
    }
    for (int tt = 1; tt <= T_simul; tt++){
        for (auto link_iterator = temp_network_copy[tt - 1].begin(); link_iterator != temp_network_copy[tt - 1].end(); link_iterator++){
            output << run << ',' << tt << ',' << (*link_iterator).w << '\n';
        } 
    }
    
    output.close();
}
