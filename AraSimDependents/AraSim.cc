#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "TTreeIndex.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h" 
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
//#include "TObject.h"

//#include <fftw3.h>

using namespace std;

#include "Tools.h"
#include "Constants.h"
#include "Vector.h"
#include "Position.h"
#include "EarthModel.h"
#include "IceModel.h"
#include "Efficiencies.h"
#include "Spectra.h"
#include "Event.h"
#include "Trigger.h"
#include "Detector.h"
#include "Settings.h"
#include "counting.hh"
#include "Primaries.h"
#include "signal.hh"
#include "secondaries.hh"

#include "Ray.h"
#include "RaySolver.h"
#include "Report.h"

//#include "UsefulIcrrStationEvent.h"

#ifdef ARA_UTIL_EXISTS
#include "UsefulIcrrStationEvent.h"
ClassImp(UsefulIcrrStationEvent);
#endif

class EarthModel; //class

void test();

string outputdir="outputs";

//--------------------------------------------------
// extern"C" {
//     void model_(int *ii);
// }
//-------------------------------------------------- 


//int main() {
int main(int argc, char **argv) {   // read setup.txt file



    // below is replace by Settings class Initialize() member.
/*    
  const int NNU=100;

  // NEED TO FIGURE OUT A GOOD WAY TO READ THIS IN AND STORE THEM.
  // INPUT FILE AGAIN?  SOMETHING ELSE?
  //These were moved here from IceModel under the new compilation scheme
  int ICE_MODEL=0; //Select ice model to be used.  0 = Crust 2.0 , 1 = BEDMAP.
  int NOFZ=1; // 1=depth dependent index of refraction,0=off
  int CONSTANTCRUST=0; // set crust density and thickness to constant values.
  int CONSTANTICETHICKNESS=0; // set ice thickness to constant value
  int FIXEDELEVATION=0; // fix the elevation to the thickness of ice.
  int MOOREBAY=0; //1=use Moore's Bay measured ice field attenuation length for the west land, otherwise use South Pole data
  double EXPONENT=19.; // 10^19 eV neutrinos only
*/  

  Settings *settings1 = new Settings();


  cout<<"\n\tDefault values!"<<endl;
  cout<<"NNU : "<<settings1->NNU<<endl;
  cout<<"ICE_MODEL : "<<settings1->ICE_MODEL<<endl;
  cout<<"NOFZ : "<<settings1->NOFZ<<endl;
  cout<<"CONSTANTICETHICKNESS : "<<settings1->CONSTANTICETHICKNESS<<endl;
  cout<<"FIXEDELEVATION : "<<settings1->FIXEDELEVATION<<endl;
  cout<<"MOOREBAY : "<<settings1->MOOREBAY<<endl;
  cout<<"EXPONENT : "<<settings1->EXPONENT<<endl;
  cout<<"DETECTOR : "<<settings1->DETECTOR<<endl;

  //string setupfile = "setup.txt";
  string setupfile;
  string run_no;
  string antennaVpolfile = "ARA_bicone6in_output.txt";

  if (argc<2) { // no setup file input, use default
      setupfile = "setup.txt";
      cout<<"setupfile : "<<setupfile<<endl;
  }
  else if (argc > 1) { // read file!!
      setupfile = string( argv[1] );
      cout<<"setupfile : "<<setupfile<<endl;
  }
  if (argc > 2) { // read file!!
     
      run_no = string( argv[2] );
      cout<<"run number : "<<run_no<<endl;
  }
  if (argc > 3) { // read file!!
      
      outputdir = string( argv[3] );
      if(outputdir[outputdir.size()-1]=='/') outputdir=outputdir.substr(0,outputdir.size()-1); // make sure outputdir doesn't have a / at the end
      cout<<"outputdir : "<<outputdir<<endl;
  }
  if (argc > 4) { 

      antennaVpolfile = string( argv[4] );
      cout<<"antenna file : "<<antennaVpolfile <<endl;
  }

//   else { // no mode for argc > 2!
//       cout<<"too many info! just use default setup.txt file!"<<endl;
//       setupfile = "setup.txt";
//   }

  settings1->ReadFile(setupfile);
  cout<<"Read "<<setupfile<<" file!"<<endl;


  cout<<"\n\tNew values!"<<endl;
  cout<<"NNU : "<<settings1->NNU<<endl;
  cout<<"ICE_MODEL : "<<settings1->ICE_MODEL<<endl;
  cout<<"NOFZ : "<<settings1->NOFZ<<endl;
  cout<<"CONSTANTICETHICKNESS : "<<settings1->CONSTANTICETHICKNESS<<endl;
  cout<<"FIXEDELEVATION : "<<settings1->FIXEDELEVATION<<endl;
  cout<<"MOOREBAY : "<<settings1->MOOREBAY<<endl;
  cout<<"EXPONENT : "<<settings1->EXPONENT<<endl;
  cout<<"DETECTOR : "<<settings1->DETECTOR<<endl;
  cout<<"POSNU_RADIUS : "<<settings1->POSNU_RADIUS<<endl;



  // set gRandom as TRandom3 when settings1->RANDOM_MODE = 1
  if (settings1->RANDOM_MODE == 1) {

    // test TRandom3
    TRandom3 *test_randm3 = new TRandom3 (0);
    gRandom = test_randm3;
  }
    //cout<<"first random from TRandom3 : "<<test_randm3->Rndm()<<"\n";
    cout<<"first random : "<<gRandom->Rndm()<<"\n";







//  IceModel *icemodel=new IceModel(ICE_MODEL + NOFZ*10,CONSTANTICETHICKNESS * 1000 + CONSTANTCRUST * 100 + FIXEDELEVATION * 10 + 0,MOOREBAY);// creates Antarctica ice model
  IceModel *icemodel=new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10,settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0,settings1->MOOREBAY);// creates Antarctica ice model
  //IceModel inherits from EarthModel  

  cout<<endl;
  cout<<"Surface at (log:0, lat:0) : "<<icemodel->Surface(0., 0.)<<endl;
  cout<<"SurfaceAboveGeoid at (log:0, lat:0) : "<<icemodel->SurfaceAboveGeoid(0., 0.)<<endl;
  
  //Detector *detector=new Detector(settings1, icemodel); // builds antenna array, 0 for testbed
  Detector *detector=new Detector(settings1, icemodel, setupfile, antennaVpolfile); // builds antenna array, 0 for testbed
  cout<<"end calling detector"<<endl;
//  Detector *detector=new Detector(settings1->DETECTOR); // builds antenna array, 0 for testbed

  Trigger *trigger=new Trigger(detector, settings1); // builds the trigger  
//  Efficiencies *efficiencies=new Efficiencies(detector->getnRx(),outputdir); // keeps track of efficiencies at each stage of the simulation
  Efficiencies *efficiencies=new Efficiencies(100,outputdir); // keeps track of efficiencies at each stage of the simulation
  cout<<"called Efficiencies"<<endl;
  
  Spectra *spectra=new Spectra(settings1->EXPONENT); // gets library (or whatever) of neutrino spectra
  cout<<"called Spectra"<<endl;

  Ray *ray = new Ray(); // construct Ray class
  cout<<"called Ray"<<endl;
  

    // 
    // test PickUnbiased in IceModel.
  Counting *count1 = new Counting();
  cout<<"called Counting"<<endl;

  Primaries *primary1 = new Primaries();
  cout<<"called Primaries"<<endl;

  int whichray = 0; // for test
//--------------------------------------------------
//   Interaction *interaction1=new Interaction("nu",primary1,settings1,whichray,count1);
//   cout<<"called Interaction1"<<endl;
//-------------------------------------------------- 

  Event *event = new Event();
  cout<<"called Event"<<endl;

  //Report *report = new Report(detector, settings1);
  Report *report = new Report();
  cout<<"called Evt"<<endl;


  TFile *AraFile;
   if (argc > 2) {
        AraFile=new TFile((outputdir+"/AraOut."+setupfile.substr(setupfile.find_last_of("/")+1)+".run"+run_no+".root").c_str(),"RECREATE","ara");
   }
   else {
        AraFile=new TFile((outputdir+"/AraOut.root").c_str(),"RECREATE","ara");
   }

  TTree *AraTree=new TTree("AraTree","AraTree");    // for single entry
  TTree *AraTree2=new TTree("AraTree2","AraTree2"); //for many entries
  cout<<"assign AraFile, AraTrees"<<endl;

  AraTree->Branch("detector",&detector);
  cout<<"branch detector"<<endl;
  AraTree->Branch("icemodel",&icemodel);
  cout<<"branch icemodel"<<endl;
  AraTree->Branch("trigger",&trigger);
  cout<<"branch trigger"<<endl;
  AraTree->Branch("settings",&settings1);
  cout<<"branch settings"<<endl;
  AraTree->Branch("spectra",&spectra);
  cout<<"branch spectra"<<endl;
  AraTree2->Branch("event",&event);
  cout<<"branch Evt"<<endl;
  AraTree2->Branch("report",&report);
  cout<<"branch report"<<endl;

  cout<<"finish tree assign"<<endl;


RaySolver *raysolver = new RaySolver();
cout<<"called RaySolver"<<endl;

    cout << "Make output file that is readable by AraRoot" << endl;

#ifdef ARA_UTIL_EXISTS
    UsefulIcrrStationEvent *theEvent = 0;
    
    TTree *eventTree;
    eventTree = new TTree("eventTree","Tree of ARA Events");
    //eventTree->Branch("event",&theEvent);
    eventTree->Branch("UsefulARAStationEvent",&theEvent);
#endif


cout<<"will call secondaries"<<endl;
Secondaries *sec1 = new Secondaries (settings1);
cout<<"will call signal"<<endl;
//Signal *signal = new Signal;
Signal *signal = new Signal (settings1);
signal->SetMedium(0);   // set medium as ice
cout<<"finish calling secondaries and signal"<<endl;





//--------------------------------------------------
// TH1F *hy=new TH1F("hy","hy",100,0.,1.); // histogram for inelasticity
//-------------------------------------------------- 

// before start looping events set noise values (this case, thermal)
trigger->SetMeanRmsDiode(settings1, detector, report);

if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0) {// noise waveforms will be generated for each evts
    trigger->ClearNoiseWaveforms();
}

// now in Trigger class, there will be meandiode, rmsdiode values for noise (we need this for trigger later)


double max_dt = 0.; // max arrival time difference

int Total_Global_Pass = 0;  // total global trigger passed number 
double Total_Weight = 0.;
double Total_Probability = 0.;


/*

 TCanvas *cFull_window = new TCanvas("cFull_window","A Simple Graph Example",200,10,10000,11200);
 cFull_window->Divide(1,16);

 TGraph *g_Full_window;

TGraph *G_V_threshold_diode;
G_V_threshold_diode = new TGraph(2, threshold_x, threshold_y);

 TCanvas *cFull_window_V = new TCanvas("cFull_window_V","A Simple Graph Example",200,10,3200,2400);
 cFull_window_V->Divide(4,4);

 TGraph *g_Full_window_V;

 */






 double x_V[settings1->NFOUR/2];
 double y_V[settings1->NFOUR/2];



 double xbin[settings1->DATA_BIN_SIZE];
 for (int i=0; i<settings1->DATA_BIN_SIZE; i++) {
     xbin[i] = i;
 }

double threshold_y[2];
double threshold_x[2];

threshold_x[0] = 0.;
threshold_x[1] = (double)settings1->DATA_BIN_SIZE-1.;
threshold_y[0] = (trigger->rmsdiode) * (trigger->powerthreshold);
threshold_y[1] = (trigger->rmsdiode) * (trigger->powerthreshold);





cout<<"powerthreshold : "<<trigger->powerthreshold<<endl;



int check_station_DC;

    ofstream TrigWind;
    TrigWind.open("outputs/TrigWindowStudy.txt");
    
    
//    for (int iTrigWind = 10; iTrigWind < 201; iTrigWind = iTrigWind+10){
    
//        double TRIG_WINDOW_Size = double (iTrigWind * 1.0E-9);
//        settings1->TRIG_WINDOW = TRIG_WINDOW_Size;
        
        Total_Global_Pass = 0;
cout<<"begin looping events!!"<<endl;


double pre_posnu_x;
double pre_posnu_y;
double pre_posnu_z;

double cur_posnu_x;
double cur_posnu_y;
double cur_posnu_z;


    cout << "Calpulser_on: " << settings1->CALPULSER_ON << endl;

    // test Detector set correctly
    cout<<"number of stations : "<<detector->params.number_of_stations << endl;
    cout<<"total number of antennas : "<<detector->params.number_of_antennas << endl;
    int ch_count = 0;
    for (int i=0; i<detector->params.number_of_stations; i++) {
        for (int j=0; j<detector->stations[i].strings.size(); j++) {
            for (int k=0; k<detector->stations[i].strings[j].antennas.size(); k++) {
                ch_count++;
                cout<<"station["<<i<<"].strings["<<j<<"].antennas["<<k<<"] no_ch:"<<ch_count<<endl;
            }
        }
    }



        
    // check if settings have to compatibility problems
    // if there's any, stop AraSim
    int settings_compatibility_error = settings1->CheckCompatibilities(detector);
    if (settings_compatibility_error > 0) {
        cerr<<"There are "<< settings_compatibility_error<<" errors from settings. Check error messages."<<endl;
        return -1;
    }
    
#ifndef ARA_UTIL_EXISTS
    if (settings1->DETECTOR == 3 && settings1->READGEOM == 1){
        cerr << "ERROR::InstalledStation geometry not available without AraRoot installation!" << endl;
        return -1;
    }
#endif


// test autoflush
//AraTree2->SetAutoFlush(0);


    // reset accumulative trig search bin info 
    settings1->ACCUM_TRIG_SEARCH_BINS_STATION0 = 0.;


    int nuLimit =0;
    if (settings1->ONLY_PASSED_EVENTS == 1){
      nuLimit = settings1->NNU_PASSED;
    } else {
      nuLimit = settings1->NNU;
    }
    
    int inu = 0;
    int Events_Thrown = 0;
    int Events_Passed = 0;
    //       for (int inu=0;inu<settings1->NNU;inu++) { // loop over neutrinos
    while (inu < nuLimit){

        check_station_DC = 0;
        check_station_DC = 0;
        
       if ( settings1->DEBUG_MODE_ON==0 ) {
           std::cerr<<"*";
           if ( Events_Thrown%100 == 0 )
               cout<<"Thrown "<<Events_Thrown<<endl;
       }



       //event = new Event ( settings1, spectra, primary1, icemodel, detector, signal, sec1 );
       event = new Event ( settings1, spectra, primary1, icemodel, detector, signal, sec1, Events_Thrown );
        event->inu_passed = -1;
        
       report = new Report(detector, settings1);

#ifdef ARA_UTIL_EXISTS
       theEvent = new UsefulIcrrStationEvent();
#endif


       // go further only if we picked up usable posnu
       if (event->Nu_Interaction[0].pickposnu>0) {


           /*
           if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0) {// noise waveforms will be generated for each evts
               trigger->GetNewNoiseWaveforms(settings1, detector, report);
           }
           */


    //--------------------------------------------------
    //        cout<<"inu : "<<inu<<endl;
    //        cout<<"event->pnu : "<<event->pnu<<endl;
    //        cout<<"posnu : ";
    //        event->Nu_Interaction[0].posnu.Print();
    //        cout<<"nnu : ";
    //        event->Nu_Interaction[0].nnu.Print();
    //        cout<<"event->n_interactions : "<<event->n_interactions<<endl;
    //        cout<<"nu_flavor : "<<event->nuflavor<<endl;
    //        cout<<"event->Nu_Interaction[0].vmmhz1m[0] : "<<event->Nu_Interaction[0].vmmhz1m[0]<<endl;
    //        cout<<"pickposnu : "<<event->Nu_Interaction[0].pickposnu<<endl;
    //-------------------------------------------------- 

           // connect Interaction class (nu interaction with ice) and Detector class (detector properties and layout)
           // save signal, noise at each antennas to Report class
           //report->Connect_Interaction_Detector (event, detector, raysolver, signal, icemodel, settings1, trigger);

           //report->Connect_Interaction_Detector (event, detector, raysolver, signal, icemodel, settings1, trigger, theEvent);
           report->Connect_Interaction_Detector (event, detector, raysolver, signal, icemodel, settings1, trigger, Events_Thrown);
           //report->Connect_Interaction_Detector (event, detector, raysolver, signal, icemodel, settings1, trigger, theEvent, Events_Thrown);
                      
#ifdef ARA_UTIL_EXISTS
           int stationID = 0;
           if (report->stations[stationID].Global_Pass) {
               report->MakeUsefulEvent(detector, settings1, trigger, stationID, theEvent);
           }
#endif
               
           report->ClearUselessfromConnect(detector, settings1, trigger);


           // test if posnu is exactly same in calpulser mode
           if (settings1->CALPULSER_ON == 1) {
               cur_posnu_x = event->Nu_Interaction[0].posnu.GetX();
               cur_posnu_y = event->Nu_Interaction[0].posnu.GetY();
               cur_posnu_z = event->Nu_Interaction[0].posnu.GetZ();
               //cout<<"posnu x:"<<cur_posnu_x<<" y:"<<cur_posnu_y<<" z:"<<cur_posnu_z<<endl;
               if (inu>0) {
                   if (pre_posnu_x==cur_posnu_x && pre_posnu_y==cur_posnu_y && pre_posnu_z==cur_posnu_z) {
                   }
                   else cout<<"posnu location changed!"<<endl;
               }
               pre_posnu_x = event->Nu_Interaction[0].posnu.GetX();
               pre_posnu_y = event->Nu_Interaction[0].posnu.GetY();
               pre_posnu_z = event->Nu_Interaction[0].posnu.GetZ();
           }







           for (int i=0; i<detector->params.number_of_stations; i++) {
#ifdef ARA_UTIL_EXISTS
               if (settings1->DETECTOR == 3 && i == 0){ theEvent->numRFChans = 14; }
               else { theEvent->numRFChans = 16; }
#endif
               
               if (max_dt < report->stations[i].max_arrival_time - report->stations[i].min_arrival_time) max_dt = report->stations[i].max_arrival_time - report->stations[i].min_arrival_time;
               // check the total global trigger passed
               if (report->stations[i].Global_Pass) {

                   cout<<"\nGlobal_Pass : "<<report->stations[i].Global_Pass<<" evt : "<<inu<<" added weight : "<<event->Nu_Interaction[0].weight<<endl;
                   
                   if ( check_station_DC == 0) { // count trigger pass only once per event
                       
                       Total_Global_Pass ++;
                       Total_Weight += event->Nu_Interaction[0].weight;
                       Total_Probability += event->Nu_Interaction[0].probability;
                       
                       // test increment weight
                       if (settings1->INTERACTION_MODE==1) {
                           count1->incrementEventsFound( event->Nu_Interaction[0].weight, event );
                       }
                       else if (settings1->INTERACTION_MODE==0) {
                           count1->incrementEventsFound( event->Nu_Interaction[0].probability, event );
                       }
                       else if (settings1->INTERACTION_MODE==3) {
                           count1->incrementEventsFound( event->Nu_Interaction[0].probability, event );
                       }
                       
                       
                   }
                   
                   check_station_DC++;
                                      
                   event->inu_passed = Events_Passed;
                   

               }

           }



           /*
       AraTree2->Fill();   //fill interaction every events

       // for 1, save all events whether passed trigger or not
       if (settings1->WRITE_ALL_EVENTS==1) {
           //theEvent = &report->theUsefulEvent;
           eventTree->Fill();
       }
       // for 0, save events which passed trigger
       else if (settings1->WRITE_ALL_EVENTS==0) {
           //if ( Global_Pass_Flag == 1 ) {
           if ( check_station_DC > 0 ) {
               //theEvent = &report->theUsefulEvent;
               eventTree->Fill();
           }
       }


       theEvent = NULL;

 delete theEvent;
 */

            settings1->ACCUM_TRIG_SEARCH_BINS_STATION0 += report->stations[0].total_trig_search_bin;


       } // if pickposnu > 0
       else {
           //cout<<"pickposnu : "<<event->Nu_Interaction[0].pickposnu<<endl;
           report->delete_all();
           event->delete_all();
       }



       // test FILL_TREE_MODE
       if (settings1->FILL_TREE_MODE==0) { // fill event event
           
           AraTree2->Fill();   //fill interaction every events

#ifdef ARA_UTIL_EXISTS

           // for 1, save all events whether passed trigger or not
           if (settings1->WRITE_ALL_EVENTS==1) {
               //theEvent = &report->theUsefulEvent;
               eventTree->Fill();
           }
           // for 0, save events which passed trigger
           else if (settings1->WRITE_ALL_EVENTS==0) {
               if ( check_station_DC > 0 ) {
                   //theEvent = &report->theUsefulEvent;
                   eventTree->Fill();
               }
           }
#endif
       }
       else if (settings1->FILL_TREE_MODE==1) { // fill only usable posnu event
           
           if (event->Nu_Interaction[0].pickposnu>0) {
               AraTree2->Fill();   //fill interaction every events
#ifdef ARA_UTIL_EXISTS
               // for 1, save all events whether passed trigger or not
               if (settings1->WRITE_ALL_EVENTS==1) {
                   //theEvent = &report->theUsefulEvent;
                   eventTree->Fill();
               }
               // for 0, save events which passed trigger
               else if (settings1->WRITE_ALL_EVENTS==0) {
                   if ( check_station_DC > 0 ) {
                       //theEvent = &report->theUsefulEvent;
                       eventTree->Fill();
                   }
               }
#endif
           }
       }
       else if (settings1->FILL_TREE_MODE==2) { // fill only triggered event
           
           if (check_station_DC>0) {
               AraTree2->Fill();   //fill interaction every events
#ifdef ARA_UTIL_EXISTS

               // for 1, save all events whether passed trigger or not
               if (settings1->WRITE_ALL_EVENTS==1) {
                   //theEvent = &report->theUsefulEvent;
                   eventTree->Fill();
               }
               // for 0, save events which passed trigger
               else if (settings1->WRITE_ALL_EVENTS==0) {
                   if ( check_station_DC > 0 ) {
                       //theEvent = &report->theUsefulEvent;
                       eventTree->Fill();
                   }
               }
#endif
           }
       }

       
       if (settings1->ONLY_PASSED_EVENTS == 1){
           if (check_station_DC > 0){
               inu++;
           }
       } else {
           inu++;
       }
        if (check_station_DC > 0){
            Events_Passed++;
        }
       Events_Thrown++;
       

       //theEvent = NULL;





       /*
       if (settings1->NOISE_WAVEFORM_GENERATE_MODE == 0) {// noise waveforms will be generated for each evts
           // remove noise waveforms for next evt
           trigger->ClearNoiseWaveforms();
       }
       */



 delete event;
 delete report;
#ifdef ARA_UTIL_EXISTS
 delete theEvent;
#endif


  } // end loop over neutrinos

    settings1->NNU = Events_Thrown;
    settings1->NNU_PASSED = Total_Global_Pass;

//        TrigWind << TRIG_WINDOW_Size << "\t" << Total_Global_Pass << endl;
//        cout << "TRIG_WINDOW_Size:Total_Global_Pass:: " << TRIG_WINDOW_Size << " : " << Total_Global_Pass << endl;

//    }// end trigger window loop
    TrigWind.close();

                       
//--------------------------------------------------
//    cFull_window_V->Print("test_V_mimic.pdf");
//-------------------------------------------------- 





   ofstream weight_file;
   //weight_file.open(("./weight_output/weight_"+setupfile).c_str());
   if (argc == 3) {
        weight_file.open(("./weight_output/weight_"+setupfile+".run"+run_no).c_str());
   }
   else if (argc >3 ){ // add the subdirectory for outputs
        weight_file.open((outputdir+"/weight_output/weight_"+setupfile+".run"+run_no).c_str());
   }
   else {
        weight_file.open(("./weight_output/weight_"+setupfile).c_str());
   }


   cout<<" end loop"<<endl;
   cout << "Total Events Thrown: " <<     settings1->NNU << endl;
   cout<<"Total_Global_Pass : "<<Total_Global_Pass<<endl;
   cout<<"Total_Weight : "<<Total_Weight<<endl;
   cout<<"Total_Probability : "<<Total_Probability<<endl;
               
   if (settings1->INTERACTION_MODE==1) {
       weight_file << "Total_Weight="<<Total_Weight<<endl;
   }
   else if (settings1->INTERACTION_MODE==0) {
       weight_file << "Total_Probability="<<Total_Probability<<endl;
   }
   else if (settings1->INTERACTION_MODE==3) {
       weight_file << "Total_Probability="<<Total_Probability<<endl;
   }

   cout<<"weight bin values : ";
   for (int i=0; i<count1->NBINS-1; i++) {
       cout<<count1->eventsfound_binned[i]<<", ";
       weight_file << count1->eventsfound_binned[i]<<" ";
   }
       cout<<count1->eventsfound_binned[count1->NBINS-1];
       weight_file << count1->eventsfound_binned[count1->NBINS-1]<<"\n";
       weight_file.close();
   cout<<"\n\n";


   // if using picknear method
   if (settings1->INTERACTION_MODE==1) {
       double IceVolume;
       IceVolume = PI * (settings1->POSNU_RADIUS) * (settings1->POSNU_RADIUS) * icemodel->IceThickness( detector->stations[0] );
       cout<<"IceVolume : "<<IceVolume<<endl;

       double Veff_test;

       // error bar for weight
       double error_plus = 0;
       double error_minus = 0;
       Counting::findErrorOnSumWeights( count1->eventsfound_binned, error_plus, error_minus );

       /*       
       Veff_test = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * Total_Weight / (double)(settings1->NNU);

       // account all factors to error
       error_plus = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * error_plus / (double)(settings1->NNU);
       error_minus = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * error_minus / (double)(settings1->NNU);
       */

       Veff_test = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * Total_Weight / (double)(settings1->NNU);
       error_plus = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * error_plus / (double)(settings1->NNU);
       error_minus = IceVolume * 4. * PI * signal->RHOICE / signal->RHOH20 * error_minus / (double)(settings1->NNU);


       cout<<"test Veff : "<<Veff_test<<" m3sr, "<<Veff_test*1.E-9<<" km3sr"<<endl;
       cout<<"And Veff error plus : "<<error_plus*1.E-9<<" and error minus : "<<error_minus*1.E-9<<endl;
   }


   // if using pickunbiased method
   //
   else if (settings1->INTERACTION_MODE==0) {
       double IceSurf;
       IceSurf = 2. * PI * icemodel->R_EARTH * icemodel->R_EARTH * (1. - cos(icemodel->GetCOASTLINE()*RADDEG) );
       cout<<"total IceSurf : "<<IceSurf<<" m2"<<endl;

       double Aeff_test;

       // error bar for weight
       double error_plus = 0;
       double error_minus = 0;
       Counting::findErrorOnSumWeights( count1->eventsfound_binned, error_plus, error_minus );

       /*
       Aeff_test = IceSurf * PI * Total_Probability / (double)(settings1->NNU);

       // account all factors to error
       error_plus = IceSurf * PI * error_plus / (double)(settings1->NNU);
       error_minus = IceSurf * PI * error_minus / (double)(settings1->NNU);
       */

       Aeff_test = IceSurf * PI * Total_Probability / (double)(settings1->NNU);

       // account all factors to error
       error_plus = IceSurf * PI * error_plus / (double)(settings1->NNU);
       error_minus = IceSurf * PI * error_minus / (double)(settings1->NNU);


       cout<<"test Aeff : "<<Aeff_test<<" m2sr, "<<Aeff_test*1.E-6<<" km2sr"<<endl;
       cout<<"And Aeff error plus : "<<error_plus*1.E-6<<" and error minus : "<<error_minus*1.E-6<<endl;
   }


   // if using picknearunbiased method
   //
   else if (settings1->INTERACTION_MODE==3) {
       double SphereSurf;
       SphereSurf = 4. * PI * settings1->PICKNEARUNBIASED_R * settings1->PICKNEARUNBIASED_R;
       cout<<"total SphereSurf : "<<SphereSurf<<" m2"<<endl;

       double Aeff_test;

       // error bar for weight
       double error_plus = 0;
       double error_minus = 0;
       Counting::findErrorOnSumWeights( count1->eventsfound_binned, error_plus, error_minus );

       /*
       Aeff_test = IceSurf * PI * Total_Probability / (double)(settings1->NNU);

       // account all factors to error
       error_plus = IceSurf * PI * error_plus / (double)(settings1->NNU);
       error_minus = IceSurf * PI * error_minus / (double)(settings1->NNU);
       */

       Aeff_test = SphereSurf * PI * Total_Probability / (double)(settings1->NNU);

       // account all factors to error
       error_plus = SphereSurf * PI * error_plus / (double)(settings1->NNU);
       error_minus = SphereSurf * PI * error_minus / (double)(settings1->NNU);


       cout<<"test Aeff : "<<Aeff_test<<" m2sr, "<<Aeff_test*1.E-6<<" km2sr"<<endl;
       cout<<"And Aeff error plus : "<<error_plus*1.E-6<<" and error minus : "<<error_minus*1.E-6<<endl;
   }



//--------------------------------------------------
//   cout<<"Total NNU : "<<settings1->NNU<<", PickUnbiased passed NNU : "<<nnu_pass<<endl;
//-------------------------------------------------- 
    


   // remove noisewaveform info if DATA_SAVE_MODE == 2
   // remove noisewaveform info if DATA_SAVE_MODE is not 0
   if (settings1->DATA_SAVE_MODE != 0) {
       trigger->v_noise_timedomain.clear();
       trigger->v_noise_timedomain_diode.clear();
   }
   if (settings1->DATA_SAVE_MODE == 2) {// in DATA_SAVE_MODE==2, remove noise spectrum before Rayleigh dist.
       trigger->Vfft_noise_before.clear();
   }

  
  AraTree->Fill();  // fill tree for one entry

  AraFile->Write();
  //AraFile->Close();





 efficiencies->summarize(); // summarize the results in an output file  


 double freq[detector->GetFreqBin()], Filter[detector->GetFreqBin()];
 double Filter_E[detector->GetFreqBin()];

 for (int i=0; i<detector->GetFreqBin(); i++) {
     freq[i] = detector->GetFreq(i);    // in Hz
     Filter[i] = detector->GetFilterGain(i);    // in dB
     Filter_E[i] = pow(10., (detector->GetFilterGain(i))/20.);
 }




 cout<<"max_dt : "<<max_dt<<endl;


//--------------------------------------------------
// int ii = 1;
// model_(&ii);
//-------------------------------------------------- 
 
 cout<<"rmsdiode= "<<trigger->rmsdiode<<endl;

 delete raysolver;
 delete icemodel;
 delete efficiencies;
 delete ray;
 delete detector;
 delete settings1;
 delete count1;
  delete primary1;
 delete trigger;
 delete spectra;
 delete sec1;
 delete signal;


 test();

 cout<<"outputdir= "<<outputdir<<endl;


 
 return 0;
  
} //end main



void test() {

  cout << "test is " << 0 << "\n";
}

