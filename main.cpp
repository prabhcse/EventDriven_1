// Prabhmanmeet Singh

#include <iostream>
#include "rv.h"
#include "event.h"
#define A 1
#define U 2
#define D 3
using namespace std;
// Simulates an M/M/1 queueing system.  The simulation terminates
// once 100000 customers depart from the system.
int main()
{
  
  EventList Elist;                // Create event list
  enum {ARR,DEP};                 // Define the event types

  double lambda_1[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};            // ADMIN Arrival rate
  double mu = 0;                  // Service rate
  double lambda_2 = 0.0;		  // USER Arrival rate
  double clock = 0.0;             // System clock
  double prev = 0.0;			  // Previous Event time
  int N = 0;                      // Number of customers in system
  int M_CPU = 0;                  // Number of processors
  int Threshold_User = 0;         // Max limit of User jobs allowed
  int Total_Capacity = 0;		  // Max Limit of Admin jobs allowed
  int Ndep = 0;                   // Number of departures from system
  double EN[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};                // For calculating E[N], mean number of customers
  double ET[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};				// For calculating E[T], mean time spent in system
  double User_jobs[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};		    // Allowed User Jobs ( Enter )
  double Admin_jobs[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};		// Allowed Admin Jobs ( Enter )
  int go = 1;
  double Admin_Block[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};			  // Blocked Admin Jobs
  double User_Block[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};			  // Blocked User Jobs
  double Arrival_time[100000];
  double Depart_time[100000];
  int done = 0;                   // End condition satisfied?
// Edited
//while(go){
  // Initializing times
  for(int i = 0; i < 100000; i++){
  	Arrival_time[i] = 0.0;
  	Depart_time[i] = 0.0;
  }

  cout<<"Input the Lambda 2 rate: "<<endl;
  cin>>lambda_2;
  cout<<"Input the Service rate: "<<endl;
  cin>>mu;
  cout<<"Input the User Job Limit (L): "<<endl;
  cin>>Threshold_User;
  cout<<"Input the Total Capacity (K): "<<endl;
  cin>>Total_Capacity;
  cout<<"Input the Total Processors (m): "<<endl;   
  cin>>M_CPU;
  
  if(M_CPU > Total_Capacity){
    cout<<"Processors have to be less than or equal to Total Capacity.\nProgram will Exit"<<endl;
    cin.get();
    exit(0);
  }

for(int i = 0; i < 10; i++){

  lambda_1[i] = 0.1 * M_CPU * mu * (i+1);
  
  // Edited
  //double multiplier = 0.0;
  //cout<<"Multiplier for Lambda_1"<<endl;
  //cin>>multiplier;
  
  //lambda_1[i] = multiplier * M_CPU * mu * (i+1);
  // Edited
  Event* CurrentEvent;

  Elist.insert(exp_rv(lambda_1[i]),ARR,A); 		// Generate first Admin arrival event
  Admin_jobs[i]++;
  Elist.insert(exp_rv(lambda_2),ARR,U); 		// Generate first User arrival event
  User_jobs[i]++;
  int a = 0;									// iteration for arrival clock of jobs
  int d = 0;									// iteration for departure clock of jobs
  while (!done)
  {
    //CurrentEvent = Elist.get();               // Get next Event from list
    prev = clock;                        		// Store old clock value
	
	if(Elist.event_count != 0){
		if(N < Threshold_User){
			CurrentEvent = Elist.get();
			clock = CurrentEvent->time;
			//cout<< "1: "<< clock <<" Job: "<<CurrentEvent->job<<endl;
			//cin.get();
		}
		else if(N < Total_Capacity){
			CurrentEvent = Elist.get();
			
			if(CurrentEvent->job == U){			// User Job not allowed, need to skip
				User_Block[i]++;					// Create next arrival
				// Need to make another arrival
				Elist.insert(clock+exp_rv(lambda_2),ARR,U);
				continue; 
			}
			else{
				clock = CurrentEvent->time;
				//cout<< "2: "<< clock <<" Job: "<<CurrentEvent->job<<endl;
				//cin.get();
			}
		}
		else{
			CurrentEvent = Elist.get();				// N equals Total Capacity here
			if(CurrentEvent->job == U){				// User Job not allowed, need to skip
				User_Block[i]++;
				// Need to make another arrival
				Elist.insert(clock+exp_rv(lambda_2),ARR,U);
				continue; 
			}
			else if(CurrentEvent->job == A){    	// Admin job not allowed, need to skip
				Admin_Block[i]++;
				// Need to make another arrival
				Elist.insert(clock+exp_rv(lambda_1[i]),ARR,A);
				continue;
			}
			else if(CurrentEvent->job == D){
				clock = CurrentEvent->time;     	// Only Departures allowed
				//cout<< "3: "<< clock <<" Job: "<<CurrentEvent->job<<endl;
				//cin.get();
			}
		}
	}
	else{											// No events at all
		cout << "\nNo Departures and no Arrivals"<<endl;
		Elist.insert(clock+exp_rv(lambda_1[i]),ARR,A); 		// Generate first Admin arrival event
		Admin_jobs[i]++;
  		Elist.insert(clock+exp_rv(lambda_2),ARR,U); 		// Generate first User arrival event
  		User_jobs[i]++;
  		continue;
	}

    switch (CurrentEvent->type) {
    case ARR:                                 		//  If arrival 
      //cout<< "Time: "<<clock-prev<<endl;
      EN[i] += N*(clock-prev);                   //  update system statistics, need to calculate time
      N++;                                    		//  update system size
      if(CurrentEvent->job == A)
      { Elist.insert(clock+exp_rv(lambda_1[i]),ARR,A); Admin_jobs[i]++; }  //  generate Admin arrival
      if(CurrentEvent->job == U)
  	  { Elist.insert(clock+exp_rv(lambda_2),ARR,U); User_jobs[i]++;  }  //  generate User arrival
  	
      if (N <= M_CPU) {                             //  If Servers Available
        Elist.insert(clock+exp_rv(mu),DEP,D);   	//  generate Departure events
      }
      Arrival_time[a] = clock;					// Arrival time of ith customer
      a++;
      break;
    case DEP:                                 		//  If departure
      //cout<< "Time: "<<clock-prev<<endl;
      EN[i] += N*(clock-prev);                   //  update system statistics
      Ndep++;                                 		//  increment num. of departures
      if ((N-M_CPU) > 0) {                     		//  If customers in queue
        Elist.insert(clock+exp_rv(mu),DEP,D);   	//  generate next departure
      }
	  N--;                                    		//  decrement system size 
	  Depart_time[d] = clock;					// Departure time of ith customer
	  d++;
      break;
    }
    delete CurrentEvent;
    if (Ndep > 99999){
	done=1;        // End condition
	
/*    for(int j = 0; j < d; j++){
    	if(j == 0)
    	continue;
    	if(j > 0){
    		if( (Depart_time[j]-Depart_time[j-1]) <= 0){
    			cout<<"\nDEPART ERROR "<<i<<" AND "<<j<<endl;
    			continue;
    		}
    		if( (Arrival_time[j]-Arrival_time[j-1]) <= 0){
    			cout<<"\nARRIVAL ERROR "<<i<<" AND "<<j<<endl;
    			//continue;
    		}
    	}
    	
    	ET[i] += ( Depart_time[j]-Arrival_time[j] );

    }
*/    
    ET[i] = EN[i] / Ndep;							// Big Correction done for Expected time.
    EN[i] /= clock;
  }
}
  // Resetting for next iteration of lambda_1
  Ndep = 0;
  N = 0;
  prev = 0.0;
  done = 0;
  clock = 0.0;
  for(int i = 0; i < 100000; i++){
  	Arrival_time[i] = 0.0;
  	Depart_time[i] = 0.0;
  }
  Elist.clear();
}

for(int i = 0; i < 10; i++){
	
  double p0 = 0.0;
  double p1 = 0.0;
  double p2 = 0.0;
  double p3 = 0.0;
  double p4 = 0.0;
  if(Threshold_User == 1){
   p0 = 1/ ( 1 + (lambda_1[i]+lambda_2)/mu + lambda_1[i]*(lambda_1[i]+lambda_2)/(2*(mu*mu)) + 
  			(lambda_1[i]*lambda_1[i])*(lambda_1[i]+lambda_2)/(4*(mu*mu*mu)) + (lambda_1[i]*lambda_1[i]*lambda_1[i])*(lambda_1[i]+lambda_2)/(8*(mu*mu*mu*mu)) );
  
   p1 = ( (lambda_1[i]+lambda_2)/mu ) * p0;
   p2 = ( lambda_1[i] / (2*mu) ) * p1;
   p3 = ( lambda_1[i] / (2*mu) ) * p2;
   p4 = ( lambda_1[i] / (2*mu) ) * p3;
  
  double EstN = 0*p0 + 1*p1 + 2*p2 + 3*p3 + 4*p4;
  double EstT = EstN / ( mu*(p1 + 2*p2 + 2*p3 + 2*p4) );
  double U_Block = 1-p0;
  double A_Block = p4;
  cout<<"*****THEORETICAL*****"<<endl;
  cout<<"Expected number of customers: E[N]: "<<EstN<<endl;
  cout<<"Expected time spent by customer: E[T]: "<<EstT<<endl;
  cout<<"User Blocking Probability: "<<U_Block<<endl;
  cout<<"Admin Blocking Probability: "<<A_Block<<endl;
  
  }
  else if(Threshold_User == 3){
  	p0 = 1/ ( 1 + (lambda_1[i]+lambda_2)/mu + 0.5*((lambda_1[i]+lambda_2)/mu)*((lambda_1[i]+lambda_2)/mu) + 
	  		0.25*((lambda_1[i]+lambda_2)/mu)*((lambda_1[i]+lambda_2)/mu)*((lambda_1[i]+lambda_2)/mu) + 
			  0.125*(lambda_1[i]*(lambda_1[i]+lambda_2)*(lambda_1[i]+lambda_2)*(lambda_1[i]+lambda_2))/(mu*mu*mu*mu) );
	  		
	p1 = ( (lambda_1[i]+lambda_2)/mu ) * p0;
	p2 = ( (lambda_1[i]+lambda_2)/(2*mu) ) * p1;
	p3 = ( (lambda_1[i]+lambda_2)/(2*mu) ) * p2;
	p4 = ( (lambda_1[i]/(2*mu)) ) * p3;
	
  double EstN = 0*p0 + 1*p1 + 2*p2 + 3*p3 + 4*p4;
  double EstT = EstN / ( mu*(p1 + 2*p2 + 2*p3 + 2*p4) );
  double U_Block = p3+p4;
  double A_Block = p4;
  cout<<"*****THEORETICAL*****"<<endl;
  cout<<"Expected number of customers: E[N]: "<<EstN<<endl;
  cout<<"Expected time spent by customer: E[T]: "<<EstT<<endl;
  cout<<"User Blocking Probability: "<<U_Block<<endl;
  cout<<"Admin Blocking Probability: "<<A_Block<<endl;
  
  }
  cout<<"*****SIMULATED*****"<<endl;
  cout<<"Lambda: " << 0.1 * (i+1) * mu * M_CPU << endl;
  double User_Arrivals = 0.0;
  User_Arrivals = User_Block[i] + User_jobs[i];
  
  double Admin_Arrivals = 0.0;
  Admin_Arrivals = Admin_Block[i] + Admin_jobs[i];
  
  /*
  double Prob_Block_User =  User_Block[i]/(User_Block[i]+Admin_Block[i]);
  double Prob_User = User_Arrivals/(User_Arrivals + Admin_Arrivals);
  double Prob_Admin = Admin_Arrivals/(User_Arrivals + Admin_Arrivals);
  double Prob_Block_Admin = Admin_Block[i]/(User_Block[i]+Admin_Block[i]);
  
  double Prob_User_Block = ( Prob_Block_User * Prob_User ) / ( Prob_Block_User * Prob_User + Prob_Block_Admin * Prob_Admin );
  double Prob_Admin_Block = ( Prob_Block_Admin * Prob_Admin ) / ( Prob_Block_User * Prob_User + Prob_Block_Admin * Prob_Admin );
  */
  // output simulation results for N, E[N] 
  cout << "Number of User Blocked: " << User_Block[i] << endl;
  cout << "Number of Admin Blocked: " << Admin_Block[i] << endl;
  cout << "Expected number of customers (simulation): " << EN[i] << endl;
  cout << "Expected time spent by customer (simulation): " << ET[i] << endl;
  cout << "Blocked User Jobs Probability (simulation): " << User_Block[i] / (User_Block[i]+User_jobs[i]) << endl;
  cout << "Blocked Admin Jobs Probability (simulation): " << Admin_Block[i] / (Admin_Block[i]+Admin_jobs[i]) << endl;
//  cout << "Blocked User Jobs Probability: " << Prob_User_Block << endl;
//  cout << "Blocked Admin Jobs Probability: " << Prob_Admin_Block << endl;

  // output derived value for E[N]
  double rho_1 = lambda_1[i]/(M_CPU*mu); 
  cout << "RHO value: " << rho_1 << endl;
  double rho_2 = lambda_2/M_CPU*mu;
  double E_A = rho_1/(1-rho_1);
  double E_U = rho_2/(1-rho_2);
//  cout << "Expected number of Admin Jobs (analysis): " << E_A << endl;
//  cout << "Expected number of User Jobs (analysis): " << E_U << endl;
  cout<<"*****END*****\n"<<endl;
}
// Edited
//cout<<"To continue enter - 1, To exit enter - 0"<<endl;
//cin>>go;
//}
}

