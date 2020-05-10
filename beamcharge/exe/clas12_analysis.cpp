#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"


int main(int argc, char** argv) {
  // Need this to make sure root doesn't break
  ROOT::EnableThreadSafety();

  int NUM_THREADS = 2;
  if (getenv("NUM_THREADS") != NULL) NUM_THREADS = atoi(getenv("NUM_THREADS"));
  //if (NUM_THREADS > argc - NUM_THREADS) NUM_THREADS = 1;

  //std::cout<<"argc =="<<argc<<std::endl;


  // Make a vector of vectors of strings the size of the number of threads
//  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
  // Get the output file name
  std::string outfilename;

  if (argc >= 2) {
    // First argument is the output file
    outfilename = argv[1];
    // All other files are split evently by the under of threads
    double total_charge = 0;
    for (int i = 2; i < argc; i++) {
//	infilenames[i-2].push_back(argv[i]);
	auto chain = std::make_shared<TChain>("clas12");	
	chain->Add(argv[i]);
	auto data = std::make_shared<Branches12>(chain);
	int num_of_events = (int)chain->GetEntries();
	float minimum = 0;
	float maximum = 0;
	double difference = 0;
	double liveTime1 = 0;
	double liveTime2 = 0;
	double livetime = 0;
	for(int event = 0; event < num_of_events; event++) {
           chain->GetEntry(event);	
	   if(data->beamCharge()==0){continue;}
           else{
           	minimum = data->beamCharge();
		liveTime1 = data->liveTime();		
		//std::cout<<"livetime1  "<<liveTime1<<std::endl;
		break;
             }
         }
	
	for(int event = num_of_events-1; event > 0; event--) {
	   chain->GetEntry(event);
	   if(data->beamCharge()==0){continue;}
	   else{
           	maximum = data->beamCharge();
		liveTime2 = data->liveTime();
		//std::cout<<"livetime2  "<<liveTime2<<std::endl;
		break;
               }
         }
	//std::cout << "maximum   " << maximum << std::endl;
        //std::cout << "minimum   " << minimum << std::endl;
 
	if((liveTime1>0) && (liveTime2>0)){
	  	livetime = (liveTime1 + liveTime2)/2;
		maximum = maximum * livetime;
		minimum = minimum * livetime;
	   }
	 else if((liveTime1>0) && (liveTime2<0)){
	   	maximum = maximum * liveTime1;
                minimum = minimum * liveTime1;
	    } 
	 else if((liveTime1<0) && (liveTime2>0)){
                maximum = maximum * liveTime2;
                minimum = minimum * liveTime2;
            }
	else{
		maximum = maximum * 0.97;
                minimum = minimum * 0.97;
		std::cout<<"livetime error"<<std::endl;
	     }


	difference = maximum - minimum;

	//if(difference < 0 ) {
	//	std::cout << "i   " << i << std::endl;
	//}
    
	std::cout << "difference   " << difference << std::endl;
	
	total_charge +=difference; 
     }
	std::cout<<"total_charge = "<<total_charge<<std::endl;
} 
else{std::cout<<"error"<<std::endl;}
 return 0;     
}
  // Make a set of threads (Futures are special threads which return a value)
//  std::future<size_t> threads[NUM_THREADS];

