/*! \file MtasProcessor.cpp
 *
 * The MTAS processor handles detectors of type mtas  
 * cleanup for KDK marked with ** NTB
*/


//------------------------------INCLUDES-------------------------------//
#include "damm_plotids.h"
#include "param.h"
#include "MtasProcessor.h"
#include "DetectorDriver.h"
#include "RawEvent.h"

#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
//--------------------------------------------------------------------//


using std::cout;
using std::endl;
using std::vector;
using std::string;



MtasProcessor::MtasProcessor() :
    EventProcessor(),
    mtasSummary(NULL),
    siliSummary(NULL),
    logiSummary(NULL),
    refmodSummary(NULL) //Goetz added refmodSummary subtype
{
    firstTime = -1.;
    name = "mtas";
    associatedTypes.insert("mtas");
    associatedTypes.insert("sili");//silicons
}




void MtasProcessor::DeclarePlots(void) const
{
    using namespace dammIds::mtas;
    
    const int EnergyBins = SE; 
    
    DeclareHistogram1D(MTAS_POSITION_ENERGY+200, EnergyBins, "Total Mtas");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+210, EnergyBins, "Total Central");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+220, EnergyBins, "Total Inner");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+230, EnergyBins, "Total Middle");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+240, EnergyBins, "Total Outer");
    /* only filling totals	
    //beta gated TAS spectras
    DeclareHistogram1D(MTAS_POSITION_ENERGY+300, EnergyBins, "Total Mtas, Beta-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+310, EnergyBins, "Total Central, Beta-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+320, EnergyBins, "Total Inner, Beta-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+330, EnergyBins, "Total Middle, Beta-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+340, EnergyBins, "Total Outer, Beta-gated");
    
    //energy iff all central fire
    DeclareHistogram1D(MTAS_POSITION_ENERGY+201, EnergyBins, "Total Mtas, w/All Central PMTs");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+211, EnergyBins, "Total Central, w/All Central PMTs");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+301, EnergyBins, "Total Mtas, Beta-gated, w/All Central PMTs");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+311, EnergyBins, "Total Central, Beta-gated, w/All Central PMTs");
    
    //Sum spectras
    DeclareHistogram1D(MTAS_POSITION_ENERGY+225, EnergyBins, "Individual I");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+235, EnergyBins, "Individual M");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+245, EnergyBins, "Individual O");
    
    //Beta-gated sum spectras
    DeclareHistogram1D(MTAS_POSITION_ENERGY+305, EnergyBins, "Individual I,M,O, Beta-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+325, EnergyBins, "Individual I, Beta-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+335, EnergyBins, "Individual M, Beta-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+345, EnergyBins, "Individual O, Beta-gated");
    
    
    //2D spectras
    DeclareHistogram2D(MTAS_POSITION_ENERGY+350, SA, SA, "MTAS / 10  vs I, M, O / 10");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+351, SA, SA, "MTAS / 10  vs C / 10");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+352, SA, SD, "MTAS / 10  vs beta-gamma time (10 ns)");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+353, SA, SD, "MTAS / 10  vs beta-IMO time (10 ns)");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+354, S5, SD, "MTAS Module number vs beta-gamma time (10 ns)");
    
    DeclareHistogram2D(MTAS_POSITION_ENERGY+355, S9,SE, "Trigger E /10 v Total MTAS/10 ");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+356, S9,SE, "Trigger E /10 v Total MTAS/10 Region<=12");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+356, SB,SB, "Trigger TraceE /10 v Total MTAS/10 ");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+358, SA,SA, "Trigger TraceE /10 v Trigger E /10 ");
    
    DeclareHistogram2D(MTAS_POSITION_ENERGY+275, SA, S4, "Energy / 10 vs Multiplicity, Central");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+276, SA, S4, "Min Energy / 10 vs Multiplicity, Central");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+375, SA, S4, "Energy / 10 vs Multiplicity, Central Beta-gated");
    
    DeclareHistogram1D(MTAS_POSITION_ENERGY+274, S4, "Central Multiplicity");
    
    DeclareHistogram2D(MTAS_POSITION_ENERGY+360,SE,SC,"Si Raw Energy/10 vs. Time(min)");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+361,SB,SC,"Si Energy with MTAS coin. vs. Time");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+372, SB,SB, "Trigger E v TraceE/Height");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+362, SB,SB, "Trigger TraceE v TraceE/Height");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+363, S9,SB, "MnMTAS gated TraceE v TraceE/Height");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+364, S9,SB, "ZnMTAS gated TraceE v TraceE/Height");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+365, S9,SA, "Trigger TraceE v Height/TraceE");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+366, S9,SB, "Region 1 TraceE/10 v MTAS/10");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+367, S9,SB, "Region 2 TraceE/10 v MTAS/10");
    //DeclareHistogram2D(MTAS_POSITION_ENERGY+368, S9,SB, "Region 3 TraceE/10 v MTAS/10 ");
    
    DeclareHistogram1D(MTAS_POSITION_ENERGY+390,SA,"Si Saturation Length");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+391,SE,"Si Corrected Energy (w/Sat.)/10");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+392,SE,"Si Post - Pre Pulse Baseline");
    //clareHistogram2D(MTAS_POSITION_ENERGY+393,SD,SC,"Si E/10 vs. Post - Pre Baseline");
    
    //sum of the front and back photomultipier
    for(int i=0; i<24; i++)
    {
    DeclareHistogram1D(MTAS_POSITION_ENERGY+100+i, EnergyBins, "Sum F+B, Beta-gated"); 
    }
    
    //Neutron stuff
    DeclareHistogram1D(MTAS_POSITION_ENERGY+400, EnergyBins, "C B-gated N-E-T-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+401, EnergyBins, "Tot B-gated N-IMO-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+402, EnergyBins, "C B-gated N-IMO-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+403, EnergyBins, "Tot B-gated N-MO-gated");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+404, EnergyBins, "C B-gated N-MO-gated");
    
    //Silicons
    DeclareHistogram1D(MTAS_POSITION_ENERGY+500, S4, "Silicon Multiplicity");
    
    //KSI pair spectrum
    DeclareHistogram1D(MTAS_POSITION_ENERGY+501, SE, "SK1 if Mult. 2");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+502, SE, "SK2 if Mult. 2");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+503, SE, "SK1 if Mult. 1");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+504, SE, "SK2 if Mult. 1");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+505, SE, "Avg. SK1&2 if Mult. 2");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+506, SE, "Avg. SK1&2");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+507, SE, "Diff. SK1-SK2/SKTotal");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+508, SE, SC, "Diff. SK1-SK2/SKTotal vs SKTotal");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+509, SE, SC, "Diff. SK1-SK2/SKTotal vs SKTotal Anticoincidence"); // altered by Goetz 02/12
    DeclareHistogram2D(MTAS_POSITION_ENERGY+510, SE, SC, "Diff. SK1-SK2/SKTotal vs SKTotal Coincidence"); // altered by NTB 02/13
    /// comment out later
    DeclareHistogram1D(MTAS_POSITION_ENERGY+511, SE, "KSI if Region 1");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+512, SE, "KSI if Region 2");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+513, SE, "KSI if Region 3");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+514, SE, "KSI if Region 4");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+515, SE, "KSI if Region 1");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+516, SE, "KSI if Region 2");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+517, SE, "KSI if Region 3");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+518, SE, "KSI if Region 4");
    ///
    DeclareHistogram2D(MTAS_POSITION_ENERGY+511, SB , S7, "KSI E vs PMT1 Region");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+512, SB , S7, "KSI E vs PMT2 Region");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+513, SA , S7, "MTAS vs PMT1 Region");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+514, SA , S7, "MTAS vs PMT2 Region");   
    
    DeclareHistogram2D(MTAS_POSITION_ENERGY+520, SB , SE, "KSI Trace if Noise, 0 Mtas and <500 KSI");
    DeclareHistogram2D(MTAS_POSITION_ENERGY+521, SB , SE, "KSI Trace if Signal,0 Mtas and <1200 KSI");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+523, SE, "SK1 if Mult. 1 && 1300<MTAS<1600");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+524, SE, "SK2 if Mult. 1 && 1300<MTAS<1600");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+525, SE, "Avg. SK1&2/100 if Mult. 2");
    
    //Silicon Gated spectra
    DeclareHistogram1D(MTAS_POSITION_ENERGY+456, EnergyBins, "Max Si Energy");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+457, EnergyBins, "Si Gate 0");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+458, EnergyBins, "Si Gate 1");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+459, EnergyBins, "Si Gate 2");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+460, EnergyBins, "Si Gate 3");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+461, EnergyBins, "Si Gate 4");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+462, EnergyBins, "Si Gate 5");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+463, EnergyBins, "Si Gate 6");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+464, EnergyBins, "Si Gate 7");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+465, EnergyBins, "Si Gate 8");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+470, EnergyBins, "Si Gate 1");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+471, EnergyBins, "Si Gate 2");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+472, EnergyBins, "Si Gate 3");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+473, EnergyBins, "Si Gate 4");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+474, EnergyBins, "Si Gate 5");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+475, EnergyBins, "Si Gate 6");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+476, EnergyBins, "Si Gate 7");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+477, EnergyBins, "Si Gate 8");
    DeclareHistogram1D(MTAS_POSITION_ENERGY+478, EnergyBins, "Si Gate 9");
    */
}



bool MtasProcessor::Process(RawEvent &event)
{
    using namespace dammIds::mtas;
    
    if (!EventProcessor::Process(event)){
	return false;
    }

  
    //---------------------------------------------------------------------//
    // first time through, grab the according detector summaries
    if (mtasSummary == NULL){
	mtasSummary = event.GetSummary("mtas");
    }
    
    if (siliSummary == NULL){	
	siliSummary = event.GetSummary("sili");
    }
    
    vector<ChanEvent*> mtasList= mtasSummary->GetList();
    vector<ChanEvent*> siliList= siliSummary->GetList();
    
    map<string, struct MtasData> mtasMap;
    // dvm map<string, struct KSIData> siliMap;
    //---------------------------------------------------------------------//
    
    
    
    //sort data
    
    
    //------------------------------MTAS-SORT-VARIABLES--------------------//
    double maxLocation =0;
    //---------------------------------------------------------------------//

    
    //------------------------------MTAS-SORT------------------------------//
    for(vector<ChanEvent*>::const_iterator mtasListIt = mtasList.begin();
	mtasListIt != mtasList.end(); mtasListIt++){
	
	string subtype = (*mtasListIt)->GetChanID().GetSubtype();
	
	// Skip all events with multiple signals in a single event  ---dvm
	if (mtasMap.count(subtype)>0){
	    cout<<"Error: Detector "<<subtype<<" has "<<
		mtasMap.count(subtype)+1<<" signals in one event"<<endl;
      continue;
	}
	
	// Ignore all events with energy of 0 or more than 30000 ---dvm
	if( (*mtasListIt)->GetEnergy() == 0 || (*mtasListIt)->GetEnergy() > 30000 ){ 
	    continue;
	}
	
	// Pair subtype with the ?? ---dvm
	mtasMap.insert(make_pair(subtype,MtasData((*mtasListIt))));
	
	
	if( maxLocation < (*mtasListIt)->GetChanID().GetLocation() ){
	    maxLocation = (*mtasListIt)->GetChanID().GetLocation();
	}
    }
    //---------------------------------------------------------------------//



    

    
    //Matt** Pairing PMT's 
    
    // Sets first event time
    double actualTime = -1.0;   
    
    if(mtasSummary->GetMult() > 0){ //I have at least one element in mtasList
	vector<ChanEvent*>::const_iterator mtasListIt = mtasList.begin();
	actualTime = (*mtasListIt)->GetTime() * pixie::clockInSeconds;
	if (firstTime == -1.){
	    firstTime = actualTime;         //Check first is earliest. 
	}
    } else if (siliSummary->GetMult() > 0){ //I have at least one element in KSI
	vector<ChanEvent*>::const_iterator siliListIt = siliList.begin();
	actualTime = (*siliListIt)->GetTime() * pixie::clockInSeconds;
	if (firstTime == -1.){
	    firstTime = actualTime;         //Check first is earliest. 
	}
    }
    
    


  
    //------------------------------SILI-SORT-VARIABLES--------------------//

    
    //---------------------------------------------------------------------//

    //------------------------------SILI-SORT------------------------------//
    for(vector<ChanEvent*>::const_iterator siliListIt = siliList.begin();
	siliListIt != siliList.end(); siliListIt++){
	
	
	string subtype = (*siliListIt)->GetChanID().GetSubtype();
	
	int siliSaturated = (*siliListIt)->IsSaturated();
	
	if (siliSaturated) {
	    int satLen = (*siliListIt)->GetSatLen();
	    double calEnergy = (*siliListIt)->GetCalEnergy();
	    //cout << "calB " << calEnergy;
	    calEnergy=calEnergy*exp((satLen*satLen)/(8.*235.*235.));
	    //calEnergy+=calEnergy*((satLen*satLen)/(8*200*200));
	    (*siliListIt)->SetCalEnergy(calEnergy);
	    (*siliListIt)->SetEnergy(calEnergy);
	    //cout << "calA " << calEnergy << " len " << satLen << endl;
	}
	
	
	// dvm siliMap.insert(make_pair(subtype,KSIData((*siliListIt))));
    }
    //---------------------------------------------------------------------//


    double earliestIMOTime = std::numeric_limits<double>::max();
    const int numMod = 19;
    double earliestModuleTime[numMod];

    
    for (int i=0; i<numMod; i++) {
	earliestModuleTime[i] = std::numeric_limits<double>::max();
    } //moved NTB
    


    //     Matt** Loop over all APD. Sili means APD

  
    
    
    // KSI Info
    
    //------------------------------SILI-SORT-VARIABLES--------------------//
    bool isBetaSignal = false;
    double betaTime = -100.0;
    static double firstBetaTime = -1.0;
    double maxSiliconSignal = -1.0;
    //---------------------------------------------------------------------//
    
    /*double tqdcSiliconSignal = -1.0;
      double tAvgBaseline = -1.0;
      int    tTraceMax = -1.0;
      int    tRegion = -1.0;
      bool   siliSaturated =0;
      static vector <int> traceCounter(4,0);
      double siliAvgBaselinePost = -1;*/					//siliEv->GetAvgBaselinePost();
    
    // dvm int siliMult = siliMap.size();
    
    /* ONLY PLOTTING TOTALS ---DVM
       plot(MTAS_POSITION_ENERGY+500, siliMult);
    */
/* dvm    
    //From CLEAN_9 updated NTB
    if(siliMult>0){ //process TIME only of beta event.
	const ChanEvent *siliEv = siliSummary->GetMaxEvent();
	betaTime = (siliEv->GetTime()) * pixie::clockInSeconds;
	if (firstBetaTime == -1.0){
	    firstBetaTime = betaTime;         
	}
    }
*/
    double siliSumEnergy = -1.0;
    double sili1Energy = 0.0;
    double sili2Energy = 0.0;
    int sili1Region = 0;
    int sili2Region = 0;
 /* dvm    
    vector<int> trace1, trace2;
    for(map<string, struct KSIData>::const_iterator siliMapIt = siliMap.begin();
	siliMapIt != siliMap.end(); siliMapIt++){
	std::string subtype = (*siliMapIt).first;
	if (subtype == "SK1" && (*siliMapIt).second.energy >= 1. ){
	    sili1Energy = (*siliMapIt).second.calEnergy;
	    trace1 = (*siliMapIt).second.trace;
	    sili1Region = (*siliMapIt).second.region;
	    //cout << " cal 1 " << sili1Energy << endl;
	} else if (subtype == "SK2" && (*siliMapIt).second.energy >= 1. ){
	    sili2Energy = (*siliMapIt).second.calEnergy;
	    trace2 = (*siliMapIt).second.trace;
	    sili2Region = (*siliMapIt).second.region;            
	    //cout << " cal 2 " << sili2Energy << endl;
	}
	
	//plot(MTAS_POSITION_ENERGY, (*siliMapIt).second.location-76);//hardcoded numbering from map file! NTB
    }
    
    if (sili1Energy > 0.0 && sili2Energy > 0.0 ) {
	isBetaSignal = true;
	siliSumEnergy = (sili1Energy + sili2Energy)/2.;
	// cout << " Sili E " << siliSumEnergy << " sat "<< endl;
    } else if (sili1Energy > 0.0) {
	siliSumEnergy = sili1Energy;
    } else {
	siliSumEnergy = sili2Energy;
    }
 */
    /* only filling totals	
       if (siliMult ==1){
       plot(MTAS_POSITION_ENERGY+503, sili1Energy);
       plot(MTAS_POSITION_ENERGY+504, sili2Energy);
       }
    */
    
    
    double siliPDiff = 1000.0 * ( sili1Energy - sili2Energy ) / siliSumEnergy + 2000.0;
    
    //Cashew correction NTB , CG
    double corrFactor = 0.0001049*siliPDiff*siliPDiff - 0.3833 * siliPDiff + 1446;
    siliSumEnergy *= 1125/corrFactor; // corrFactor should be ~ between .5 and 3 
    
    // dvm if (siliMult > 1) {
	/* only filling totals
	   plot(MTAS_POSITION_ENERGY+501, sili1Energy);
	   plot(MTAS_POSITION_ENERGY+502, sili2Energy);
	   plot(MTAS_POSITION_ENERGY+505, siliSumEnergy);
	   plot(MTAS_POSITION_ENERGY+525, siliSumEnergy/100);
	   plot(MTAS_POSITION_ENERGY+507, siliPDiff );
	   plot(MTAS_POSITION_ENERGY+508, siliSumEnergy, siliPDiff);
	   
	   /// comment out later
	   if //(siliPDiff > 2000 && siliPDiff < 2500)
	   {
	   plot(MTAS_POSITION_ENERGY, siliSumEnergy/10 , (betaTime-firstBetaTime)/60.);  
	   }
	   ///
	   
	   plot(MTAS_POSITION_ENERGY+360, siliSumEnergy/10 , (betaTime-firstBetaTime)/60.); 
	   //CLEAN_17
	   */	    
	// dvm maxSiliconSignal = siliSumEnergy; //using avg as energy! NTB
    // dvm }
    
    
    /* only filling totals
       plot(MTAS_POSITION_ENERGY+506, siliSumEnergy);
       if (isBetaSignal){
       plot(MTAS_POSITION_ENERGY+511, siliSumEnergy/50 , sili1Region);
       plot(MTAS_POSITION_ENERGY+512, siliSumEnergy/50 , sili2Region);
       }
    */
    
    
    //spectras
    
    //0- all mtas, 1 - Central, 2 - Inner, 3 - Middle, 4 - Outer
    vector <double> totalMtasEnergy (5,-1);
    // 0-5 Central, 6-11 Inner, 12-17 Middle, 18-23 Outer
    vector <double> sumFrontBackEnergy(24,0);
    int nrOfCentralPMT = 0;
    double theSmallestCEnergy = 60000;
    
    //Matt** Building the MTAS rings: See up
    for(map<string, struct MtasData>::const_iterator mtasMapIt = mtasMap.begin();
	mtasMapIt != mtasMap.end(); mtasMapIt++){
	
	double signalEnergy = (*mtasMapIt).second.calEnergy;
	double location = (*mtasMapIt).second.location;
	double time = (*mtasMapIt).second.time * pixie::clockInSeconds;
	int moduleIndex = (location -1)/2;
	
	if((*mtasMapIt).first[0] == 'C'){
	    totalMtasEnergy.at(1) += signalEnergy/12;
	    totalMtasEnergy.at(0) += signalEnergy/12;	
	    nrOfCentralPMT++;
	    if (earliestModuleTime[0] > time ){
		earliestModuleTime[0]=time;
	    }
	    
	    if(theSmallestCEnergy > signalEnergy){
		theSmallestCEnergy = signalEnergy;
	    }
	} else if((*mtasMapIt).first[0] == 'I'){
	    totalMtasEnergy.at(2) += signalEnergy/2.;
	    totalMtasEnergy.at(0) += signalEnergy/2.;
	    if ( (moduleIndex > 5 && moduleIndex < 12) && earliestModuleTime[moduleIndex-5] > time ){
		earliestModuleTime[moduleIndex-5]=time;
	    } 
	    if (earliestIMOTime > 0 && time < earliestIMOTime){
		earliestIMOTime = time;
	    }
	} else if((*mtasMapIt).first[0] == 'M'){
	    totalMtasEnergy.at(3) += signalEnergy/2.;
	    totalMtasEnergy.at(0) += signalEnergy/2.;		
	    if ( (moduleIndex > 11 && moduleIndex < 18) && earliestModuleTime[moduleIndex-5] > time ){
		earliestModuleTime[moduleIndex-5]=time;
	    }
	    
	    if (earliestIMOTime > 0 && time < earliestIMOTime){
		earliestIMOTime = time;
	    }
	} else if((*mtasMapIt).first[0] == 'O'){
	    totalMtasEnergy.at(4) += signalEnergy/2.;
	    totalMtasEnergy.at(0) += signalEnergy/2.;
	    if ( (moduleIndex > 17 && moduleIndex < 24) && earliestModuleTime[moduleIndex-5] > time ){
		earliestModuleTime[moduleIndex-5]=time;
	    } 
	    if (earliestIMOTime > 0 && time < earliestIMOTime){
		earliestIMOTime = time;
	    }
	}
	
	
	//F+B
	
	if(moduleIndex > 24){
	    cout<<"Warning: detector "<<(*mtasMapIt).first<<" location > 48"<<endl;
	    continue;
	}
	
	if(sumFrontBackEnergy.at(moduleIndex) == 0){ //it's the first signal from this hexagon module
	    sumFrontBackEnergy.at(moduleIndex) = -1* signalEnergy/2.;
	}
	else if(sumFrontBackEnergy.at(moduleIndex) < 0){ //second signal from this hexagon module
	    sumFrontBackEnergy.at(moduleIndex) = -1*sumFrontBackEnergy.at(moduleIndex) + signalEnergy/2.;
	}
	else{ //sumFrontBackEnergy.at(moduleIndex) > 0 - 3 or more signals in one event
	    cout<<"Warning: detector "<<(*mtasMapIt).first<<" has 3 or more signals"<<endl;
	}
    }
    
    if(nrOfCentralPMT < 12 && nrOfCentralPMT > 0){
	if( totalMtasEnergy.at(1) > 0.0 ) {
	    const double oldMTASC = totalMtasEnergy.at(1);
	    totalMtasEnergy.at(1) = 12.0 * oldMTASC / ( (double ) nrOfCentralPMT );        
	    totalMtasEnergy.at(0) +=  totalMtasEnergy.at(1) - oldMTASC;
	}
    }

/* ONLY FILLING TOTAL SPECTRA
    if( nrOfCentralPMT == 12 || nrOfCentralPMT == 0){ //Total if all or nothing in Center PMTs. 
	plot(MTAS_POSITION_ENERGY+201, totalMtasEnergy.at(0) );
	if(isBetaSignal){
	    plot(MTAS_POSITION_ENERGY+301, totalMtasEnergy.at(0) );
	}
    }
    
    if( nrOfCentralPMT == 12) {
	plot(MTAS_POSITION_ENERGY+211, totalMtasEnergy.at(1) ); //Center only
	if(isBetaSignal){
	    plot(MTAS_POSITION_ENERGY+311, totalMtasEnergy.at(1) );
	}
    }
	
    if(nrOfCentralPMT < 13){ 
	plot(MTAS_POSITION_ENERGY+274, nrOfCentralPMT);
    }
        
    plot(MTAS_POSITION_ENERGY+275, totalMtasEnergy.at(1) / 10.0, nrOfCentralPMT);
    plot(MTAS_POSITION_ENERGY+276, theSmallestCEnergy / 10.0, nrOfCentralPMT);
	
    if (siliMult == 1 && totalMtasEnergy.at(0)>1300 && totalMtasEnergy.at(0) < 1600) {
	plot(MTAS_POSITION_ENERGY+523, sili1Energy);
	plot(MTAS_POSITION_ENERGY+524, sili2Energy);
    }
    //altered by Goetz 02/12
    if(isBetaSignal){
	plot(MTAS_POSITION_ENERGY+375, totalMtasEnergy.at(1) / 10.0, nrOfCentralPMT);
	if(  totalMtasEnergy.at(0) == -1 ){ 
	    plot(MTAS_POSITION_ENERGY+509, siliSumEnergy, siliPDiff ); // totalMtasEnergy.at(0) == 0 NTB
	}
	else{ //if(  totalMtasEnergy.at(0) > 0 )
	    plot(MTAS_POSITION_ENERGY+510, siliSumEnergy, siliPDiff ); // totalMtasEnergy.at(0) > 0 NTB siliSumEnergy >0 redundant?
	    if (totalMtasEnergy.at(0)>700 && totalMtasEnergy.at(0)<1200){
		plot(MTAS_POSITION_ENERGY+361, siliSumEnergy/10 , (betaTime-firstBetaTime)/60.); 
	    }
	}
    }
	
	
    for(unsigned int i=0; i<sumFrontBackEnergy.size(); i++){
	if(sumFrontBackEnergy.at(i) < 0){//it was only one signal in the hexagon module
	    sumFrontBackEnergy.at(i)=-1;  
	}
    }
*/	

    //3200 - 3240, no B-gated and 3300 - 3340, B-gated
    for(int i=0; i<5; i++){
	plot(MTAS_POSITION_ENERGY+200+i*10, totalMtasEnergy.at(i));
/* ONLY FILLING TOTALS
	if(isBetaSignal){
	    plot(MTAS_POSITION_ENERGY+300+i*10, totalMtasEnergy.at(i));
	} 
*/   
    }
	
/* ONLY FILLING TOTALS        
    //3225, 3235, 3245 - sum spectra, no B-gated and 3305, 3325, 3335, 3345, B-gated 
    for(int i=0; i<6; i++){
	plot(MTAS_POSITION_ENERGY+225, sumFrontBackEnergy.at(i+6));//Sum I
	plot(MTAS_POSITION_ENERGY+235, sumFrontBackEnergy.at(i+12));//Sum M
	plot(MTAS_POSITION_ENERGY+245, sumFrontBackEnergy.at(i+18));//Sum O
      	    
	if(isBetaSignal) {
	    plot(MTAS_POSITION_ENERGY+305, sumFrontBackEnergy.at(i+18));
	    plot(MTAS_POSITION_ENERGY+305, sumFrontBackEnergy.at(i+12));
	    plot(MTAS_POSITION_ENERGY+305, sumFrontBackEnergy.at(i+6));
	    plot(MTAS_POSITION_ENERGY+325, sumFrontBackEnergy.at(i+6));//Sum I
	    plot(MTAS_POSITION_ENERGY+335, sumFrontBackEnergy.at(i+12));//Sum M
	    plot(MTAS_POSITION_ENERGY+345, sumFrontBackEnergy.at(i+18));//Sum O
	}
    }

    if(isBetaSignal) {
	//3100 - 3146 even, B-gated
	for(int i=0; i<24; i++){
	    plot(MTAS_POSITION_ENERGY+100+2*i, sumFrontBackEnergy.at(i));
	}
	    
	//3350 mtas tot vs I, M, O, B-gated
	plot(MTAS_POSITION_ENERGY+351, totalMtasEnergy.at(0) / 10.0,totalMtasEnergy.at(1) / 10.0);
	    
	for(int i=6; i<24; i++){					
	    plot(MTAS_POSITION_ENERGY+350, totalMtasEnergy.at(0) / 10.0,sumFrontBackEnergy.at(i) / 10.0);
	}
	plot(MTAS_POSITION_ENERGY+355, totalMtasEnergy.at(0)/10 , siliSumEnergy /10 );
	if (sili1Region <12 && sili2Region <12 ) {
	    plot(MTAS_POSITION_ENERGY+356, totalMtasEnergy.at(0)/10 , siliSumEnergy /10 );
	}
	if (siliMult ==2){		
	    plot(MTAS_POSITION_ENERGY+513, totalMtasEnergy.at(0)/10 , sili1Region);
	    plot(MTAS_POSITION_ENERGY+514, totalMtasEnergy.at(0)/10 , sili2Region);        
	}
	static int rowNoise =0;
	static int rowSignal=0;            
	if (totalMtasEnergy.at(0) < 1 && siliSumEnergy < 2000){
	    if (siliSumEnergy<500) {
		for(unsigned int i = 0; i < trace1.size(); i++){
		    plot(MTAS_POSITION_ENERGY+520,i,rowNoise,trace1[i]);
		    plot(MTAS_POSITION_ENERGY+520,i,rowNoise+1,trace2[i]);
		}
		rowNoise+=2;
	    } else {
		for(unsigned int i = 0; i < trace1.size(); i++){
		    plot(MTAS_POSITION_ENERGY+ 521,i,rowSignal,trace1[i]);
		    plot(MTAS_POSITION_ENERGY+ 521,i,rowSignal+1,trace2[i]);
		}
		rowSignal+=2;
	    }
	}     
	//plot(MTAS_POSITION_ENERGY+358, maxSiliconSignal /10 , tqdcSiliconSignal /10 );
	
	double dt_beta_gamma = (actualTime - betaTime) * 1.0e8;  //convert to pixie (100MHz) ticks.
	double dt_beta_imo = (earliestIMOTime - betaTime) * 1.0e8;
	double dt_shift = 4000.0;
	plot(MTAS_POSITION_ENERGY+352, totalMtasEnergy.at(0) / 10.0, dt_beta_gamma + dt_shift);
	plot(MTAS_POSITION_ENERGY+353, totalMtasEnergy.at(0) / 10.0, dt_beta_imo + dt_shift);
	double deltaModTime;
            
	for (int i=0; i<numMod; i++) {
	    if (earliestModuleTime[i]< std::numeric_limits<double>::max() ){
		deltaModTime = ( earliestModuleTime[i] - betaTime ) * 1.0e8;
		plot(MTAS_POSITION_ENERGY+354, i, deltaModTime+dt_shift );               
	    }
	}
	    
	//These values are read from plots above - they consitite a neutron gate
	double lowNeutronE = 6500.0;
	double highNeutronE = 8000.0;
	//Notice time gate is read from dt_beta_imo plot (with shift!) and is very restrictive
	double lowNeutronT = 40.0;
	double highNeutronT = 100.0;
	    
	if (	totalMtasEnergy.at(0) > lowNeutronE && 
		totalMtasEnergy.at(0) < highNeutronE &&
		dt_beta_imo + dt_shift > lowNeutronT &&
		dt_beta_imo + dt_shift < highNeutronT
		){
	    plot(MTAS_POSITION_ENERGY+400, totalMtasEnergy.at(1));
	}            
	    
	//These neutron gates are based on neutron capture in I M O rings
	double sumMO = totalMtasEnergy.at(3) + totalMtasEnergy.at(4);
	double sumIMO = sumMO + totalMtasEnergy.at(2);
	
	if (sumIMO > lowNeutronE && sumIMO < highNeutronE) {
	    plot(MTAS_POSITION_ENERGY+401, totalMtasEnergy.at(0));
	    plot(MTAS_POSITION_ENERGY+402, totalMtasEnergy.at(1));
	}
	if (sumMO > lowNeutronE && sumMO < highNeutronE) {
	    plot(MTAS_POSITION_ENERGY+403, totalMtasEnergy.at(0));
	    plot(MTAS_POSITION_ENERGY+404, totalMtasEnergy.at(1));
	}
	    
	plot(MTAS_POSITION_ENERGY+456, siliSumEnergy/50 );
	if( maxSiliconSignal > 2000.0 && maxSiliconSignal < 4000.0 ){// greater than 2505 level feeding
	    
	    plot(MTAS_POSITION_ENERGY+457, totalMtasEnergy.at(0) );//will have comptons
	    plot(MTAS_POSITION_ENERGY+458, totalMtasEnergy.at(1) );
	}
	if(maxSiliconSignal > 0. && maxSiliconSignal < 2000.){
	    plot(MTAS_POSITION_ENERGY+459, totalMtasEnergy.at(0) );
	}
	    
	if(maxSiliconSignal > 0. && maxSiliconSignal < 2000.){
	    plot(MTAS_POSITION_ENERGY+460, totalMtasEnergy.at(1) );
	}
	    
	if(maxSiliconSignal > 4000. && maxSiliconSignal < 10000.) {          
	    plot(MTAS_POSITION_ENERGY+461, totalMtasEnergy.at(0) );//
	    plot(MTAS_POSITION_ENERGY+462, totalMtasEnergy.at(1) );//
	}
    }
	
	
	
*/
	EndProcess(); // update the processing time
	
	return true;
}


MtasProcessor::MtasData::MtasData(ChanEvent *chan)
{
  detSubtype	= chan->GetChanID().GetSubtype();
  energy	= chan->GetEnergy();
  calEnergy = chan->GetCalEnergy();
  time = chan->GetTime();
  location = chan->GetChanID().GetLocation();
}
/*
MtasProcessor::KSIData::KSIData(string type)
{
  detSubtype	= type;
  energy	= -1.;
  calEnergy = -1.;
  time = -1.;
  location = -1.;
  trace.clear();
  
}

MtasProcessor::KSIData::KSIData(ChanEvent *chan)
{
  detSubtype	= chan->GetChanID().GetSubtype();
  energy	= chan->GetEnergy();
  calEnergy = chan->GetCalEnergy();
  time = chan->GetTime();
  location = chan->GetChanID().GetLocation();
  trace = chan->GetTraceRef();
  region = chan->GetRegion();
}
*/
