/*! \file MtasProcessor.cpp
 *
 * The MTAS processor handles detectors of type mtas  */

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
 
using std::cout;
using std::endl;
using std::vector;
using std::string;

static double measureOnTime = -1.;
static double firstTime = 0.;
double previousTimeany = 0.0; // Variable to store the previous event time
double time1stBetaevent = 0.0; // Variable to store the previous bothbeta event time
double time2ndBetaevent = 0.0; // Variable to store the previous bothbeta event time
double previousTimeBeta1st = 0.0; // Variable to store the previous first beta event time
double previousTime2nd = 0.0; // Variable to store the previous any event second to beta event time
double previousTime1stBeta = 0.0; // Variable to store the previous 1st beta event time
double previousTime2nd_notBeta = 0.0; // Variable to store the not beta event second to beta event time
bool is1stBetaevent = false ; // flag for the first beta event 
bool MtasProcessor::isTapeMoveOn = false;
bool MtasProcessor::isMeasureOn = true;
bool MtasProcessor::isBkgOn = false;
bool MtasProcessor::isLightPulserOn = false;
bool MtasProcessor::isIrradOn = false;
double MtasProcessor::measureOnTime = -1; 
unsigned MtasProcessor::cycleNumber = 0;

MtasProcessor::MtasProcessor():EventProcessor(), mtasSummary(NULL), siliSummary(NULL), geSummary(NULL), sipmSummary(NULL), logiSummary(NULL), refmodSummary(NULL){ //Goetz added refmodSummary subtype
	firstTime = -1.;	//IS THIS THE SAME AS THE global static double firstTime ? DOUBLE CHECK CPP RULES
	name = "mtas";
	associatedTypes.insert("mtas");
	associatedTypes.insert("sili");//silicons
	associatedTypes.insert("ge");
	associatedTypes.insert("mtaspspmt");
	associatedTypes.insert("logi");//logic sygnals
	//associatedTypes.insert("ionc");//ionization chamber
	// Goetz: added reference module type generic for March 2015 experiment 
	associatedTypes.insert("refmod"); 
}

void MtasProcessor::DeclarePlots(void) const{
	using namespace dammIds::mtas;
    
	const int EnergyBins = SE; 
	//MTAS spectras
	for(int i=0; i<4; i++){    	
		string titlePart;
		if(i == 1)
			titlePart = ", Background";
		if(i == 2)	
			titlePart = ", Light Pulser";
		if(i == 3)
			titlePart = ", Irradiation";       
		DeclareHistogram1D(MTAS_POSITION_ENERGY+200+i, EnergyBins, ("Total Mtas"+titlePart).c_str());
		DeclareHistogram1D(MTAS_POSITION_ENERGY+210+i, EnergyBins, ("Total Central"+titlePart).c_str());
		DeclareHistogram1D(MTAS_POSITION_ENERGY+220+i, EnergyBins, ("Total Inner"+titlePart).c_str());
		DeclareHistogram1D(MTAS_POSITION_ENERGY+230+i, EnergyBins, ("Total Middle"+titlePart).c_str());
		DeclareHistogram1D(MTAS_POSITION_ENERGY+240+i, EnergyBins, ("Total Outer"+titlePart).c_str());
	}
	
	//beta gated TAS spectras
	for(int i=0; i<5; i++){    	
		string titlePart;
		if(i == 1)
			titlePart = ", Background";
		if(i == 2)	
			continue;
		if(i == 3)
			titlePart = ", Irradiation"; 
		if(i == 4)
			titlePart = ", Irr Bkg";	      
		DeclareHistogram1D(MTAS_POSITION_ENERGY+300+i, EnergyBins, ("Total Mtas, Beta-gated"+titlePart).c_str());
		DeclareHistogram1D(MTAS_POSITION_ENERGY+310+i, EnergyBins, ("Total Central, Beta-gated"+titlePart).c_str());
		DeclareHistogram1D(MTAS_POSITION_ENERGY+320+i, EnergyBins, ("Total Inner, Beta-gated"+titlePart).c_str());
		DeclareHistogram1D(MTAS_POSITION_ENERGY+330+i, EnergyBins, ("Total Middle, Beta-gated"+titlePart).c_str());
		DeclareHistogram1D(MTAS_POSITION_ENERGY+340+i, EnergyBins, ("Total Outer, Beta-gated"+titlePart).c_str());
	}
	
	//Sum spectras
	DeclareHistogram1D(MTAS_POSITION_ENERGY+200+25, EnergyBins, "Sum I");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+200+35, EnergyBins, "Sum M");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+200+45, EnergyBins, "Sum O");
	
	//Full Story spectras
	vector<string> titlePart = {"Total MTAS","Total Central","Total Inner","Total Middle","Total Outer"};
	for (int i=0; i<5; i++){
		DeclareHistogram2D(MTAS_EVO_NOLOGIC+i,EnergyBins,SB,(titlePart.at(i)+"Full Story, No logic").c_str());
		DeclareHistogram2D(MTAS_EVO_NOLOGIC+i+10,EnergyBins,SB,(titlePart.at(i)+"Full Story, Beta-gated").c_str());
		DeclareHistogram2D(MTAS_EVO_NOLOGIC+i+20,EnergyBins,S9,(titlePart.at(i)+"Full Story, Beta-gated, cycle time").c_str());
	}
	//THOMAS RULAND PLOTS
	DeclareHistogram1D(MTAS_POSITION_ENERGY+600, EnergyBins, "250-370 m spectra, Beta-gated");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+601, EnergyBins, "1363-1483 min spectra, Beta-gated");

	//Beta-gated sum spectras
	DeclareHistogram1D(MTAS_POSITION_ENERGY+300+5, EnergyBins, "Sum I+M+O, Beta-gated");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+300+25, EnergyBins, "Sum I, Beta-gated");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+300+35, EnergyBins, "Sum M, Beta-gated");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+300+45, EnergyBins, "Sum O, Beta-gated");

	//Akhil's plots for time difference
	DeclareHistogram1D(MTAS_POSITION_ENERGY+700, SE, "TDiff Both Beta");
	//DeclareHistogram1D(MTAS_POSITION_ENERGY+701, SE, "TDiff 1st Beta only");
	//DeclareHistogram1D(MTAS_POSITION_ENERGY+702, SE, "TDiff any two events");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+703, SE, "TDiff 1st Beta and 2nd !beta");

	//2d plots with time difference for "TDiff 1st Beta only"
	DeclareHistogram2D(MTAS_POSITION_ENERGY+704, SA, SA, "MTAS / 10  vs I, M, O / 10");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+705, SA, SA, "MTAS / 10  vs C / 10");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+706, SA, S9, "MTAS / 10  vs cycle time");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+707, SA, S9, "MTAS C / 10  vs cycle time");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+708, SA, SE, "MTAS / 10  vs time diff (1st beta)");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+709, SA, SE, "MTAS C / 10  vs time diff (1st beta)");
	// DeclareHistogram2D(MTAS_POSITION_ENERGY+361, SA, S9, "MTAS / 10  vs cycle time 100ms/bin, B-gated");
	// DeclareHistogram2D(MTAS_POSITION_ENERGY+362, SA, S9, "MTAS C/ 10  vs cycle time , B-gated");

	//Plots with time difference "TDiff 1st Beta and 2nd !beta"
	//string titlePart1stBnB = "TDiff BnB: Ist Beta";
	DeclareHistogram1D(MTAS_POSITION_ENERGY+250, EnergyBins, "Total Mtas,TDiff BnB: Ist Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+251, EnergyBins, "Total Central,TDiff BnB: Ist Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+252, EnergyBins, "Total Inner,TDiff BnB: Ist Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+253, EnergyBins, "Total Middle,TDiff BnB: Ist Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+254, EnergyBins, "Total Outer,TDiff BnB: Ist Beta");

	//string titlePart2ndBnB = ", TDiff BnB: 2nd !Beta";	
	DeclareHistogram1D(MTAS_POSITION_ENERGY+255, EnergyBins, "Total Mtas,TDiff BnB: 2nd !Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+256, EnergyBins, "Total Central,TDiff BnB: 2nd !Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+257, EnergyBins, "Total Inner,TDiff BnB: 2nd !Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+258, EnergyBins, "Total Middle,TDiff BnB: 2nd !Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+259, EnergyBins, "Total Outer,TDiff BnB: 2nd !Beta");
	// DeclareHistogram2D(MTAS_POSITION_ENERGY+710, SA, SA, "MTAS / 10  vs I, M, O / 10 b/!b");
	// DeclareHistogram2D(MTAS_POSITION_ENERGY+711, SA, SA, "MTAS / 10  vs C / 10 b/!b");
	// DeclareHistogram2D(MTAS_POSITION_ENERGY+712, SA, S9, "MTAS / 10  vs cycle time b/!b");
	// DeclareHistogram2D(MTAS_POSITION_ENERGY+713, SA, S9, "MTAS C / 10  vs cycle time b/!b");
	// DeclareHistogram2D(MTAS_POSITION_ENERGY+714, SA, SE, "MTAS / 10  vs time diff b/!b ");
	// DeclareHistogram2D(MTAS_POSITION_ENERGY+715, SA, SE, "MTAS C / 10  vs time diff b/!b");

	//Plots with time difference "TDiff 1st Beta and 2nd beta"
	//string titlePart1stBnB = "TDiff BnB: Ist Beta";
	DeclareHistogram1D(MTAS_POSITION_ENERGY+280, EnergyBins, "Total Mtas,TDiff BB: Ist Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+281, EnergyBins, "Total Central,TDiff BB: Ist Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+282, EnergyBins, "Total Inner,TDiff BB: Ist Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+283, EnergyBins, "Total Middle,TDiff BB: Ist Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+284, EnergyBins, "Total Outer,TDiff BB: Ist Beta");

	//string titlePart2ndBnB = ", TDiff BnB: 2nd !Beta";	
	DeclareHistogram1D(MTAS_POSITION_ENERGY+285, EnergyBins, "Total Mtas,TDiff BB: 2nd Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+286, EnergyBins, "Total Central,TDiff BB: 2nd Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+287, EnergyBins, "Total Inner,TDiff BB: 2nd Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+288, EnergyBins, "Total Middle,TDiff BB: 2nd Beta");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+289, EnergyBins, "Total Outer,TDiff BB: 2nd Beta");


	//2D spectras
	DeclareHistogram2D(MTAS_POSITION_ENERGY+350, SA, SA, "MTAS / 10  vs I, M, O / 10");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+351, SA, SA, "MTAS / 10  vs C / 10");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+260, SA, S9, "MTAS / 10  vs cycle time");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+262, SA, S9, "MTAS C / 10  vs cycle time");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+360, SA, S9, "MTAS / 10  vs cycle time, MSU, BGD, B-gated");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+361, SA, S9, "MTAS / 10  vs cycle time 100ms/bin, B-gated");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+362, SA, S9, "MTAS C/ 10  vs cycle time , B-gated");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+363, SA, S9, "MTAS C / 10  vs cycle time 100ms/bin, MSU, BGD, B-gated");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+364, SA, S9, "MTAS / 10  vs cycle time 1 min/bin, B-gated");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+365, SA, S8, "MTAS / 10  vs beta-gamma time (10 ns)");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+366, SA, S8, "MTAS / 10  vs beta-IMO time (10 ns)");

	DeclareHistogram2D(MTAS_POSITION_ENERGY+355, SA, SA, "Gamma-Gamma Matrix / 10");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+275, SA, S4, "Energy / 10 vs Multiplicity, Central");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+276, SA, S4, "Min Energy / 10 vs Multiplicity, Central");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+277, SA, S4, "Max Energy / 10 vs Multiplicity, Central");
	//cycle number
	DeclareHistogram2D(MTAS_POSITION_ENERGY+270, SA, SA, "MTAS / 10  vs cycle nymber");	
	DeclareHistogram2D(MTAS_POSITION_ENERGY+271, S4, SA, "Strip Number  vs cycle nymber");
	DeclareHistogram2D(MTAS_POSITION_ENERGY+272, S4, SA, "Logic signals vs cycle number");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+273, EnergyBins, "MTAS Total scalar rate (s)");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+274, S4, "Central Multiplicity");
	//sum of the front and back photomultipier
	for(int i=0; i<24; i++){
		DeclareHistogram1D(MTAS_POSITION_ENERGY+100+3*i, EnergyBins, "Sum F+B"); 
		DeclareHistogram1D(MTAS_POSITION_ENERGY+100+3*i+1, EnergyBins, "Sum F+B, Beta-gated"); 
		DeclareHistogram1D(MTAS_POSITION_ENERGY+100+3*i+2, EnergyBins, "Sum F+B, Light Pulser");
	}

	//Neutron stuff
	DeclareHistogram1D(MTAS_POSITION_ENERGY+400, EnergyBins, "C B-gated N-E-T-gated");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+401, EnergyBins, "Tot B-gated N-IMO-gated");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+402, EnergyBins, "C B-gated N-IMO-gated");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+403, EnergyBins, "Tot B-gated N-MO-gated");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+404, EnergyBins, "C B-gated N-MO-gated");

	//Ge stuff
	DeclareHistogram2D(MTAS_POSITION_ENERGY+416, EnergyBins, SA, "Ge energy vs cycle time");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+417, EnergyBins, "Implant gated Ge");
	
	//Silicons
	DeclareHistogram1D(MTAS_POSITION_ENERGY+500, S4, "Silicon Multiplicity");
	DeclareHistogram1D(MTAS_POSITION_ENERGY+510, S4, "Silicon Position, left chan - top detector");
	//DeclareHistogram2D(MTAS_POSITION_ENERGY+520, S4, SC, "Energy vs Silicon Position");	
	
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
}

bool MtasProcessor::Process(RawEvent &event)
{
	using namespace dammIds::mtas;
	if (!EventProcessor::Process(event))
		return false;

	// first time through, grab the appropriate detector summaries
	if (mtasSummary == NULL)
		mtasSummary = event.GetSummary("mtas");//czy musze to robic tak? //Do I have to do this?
	if (siliSummary == NULL)	
		siliSummary = event.GetSummary("sili");
	if (geSummary == NULL)
		geSummary = event.GetSummary("ge");
	if (logiSummary == NULL)
		logiSummary = event.GetSummary("logi");
	if (sipmSummary == NULL)	
		sipmSummary = event.GetSummary("mtaspspmt");
	//if (refmodSummary == NULL)
	//	refmodSummary = event.GetSummary("refmod"); //added by Goetz

	mtasList = mtasSummary->GetList();
	siliList = siliSummary->GetList();
	geList = geSummary->GetList();
	sipmList = sipmSummary->GetList();
	logiList = logiSummary->GetList();
	refmodList = event.GetSummary("refmod")->GetList(); //added by Goetz
	
	//Map structures (class MtasData) are init'd in the header and emptied here at the beginning of the fill stage like they should be
	FillMtasMap();
	FillSiliMap();
	FillGeMap();
	FillSipmMap();
 	FillRefModMapAndEnergy();
	FillLogicMap(); //Maps are filled with data

	//Set up time for mtas and cycle logic
	double cycleTime = -1.0;
	double actualTime = -1.0;
	double earliestIMOTime = std::numeric_limits<double>::max();
	double actualLogiTime = -1.0;
	double cycleLogiTime = -1.0;
	

	if(mtasSummary->GetMult() > 0)//I have at leat one element in mtasList
	{
		vector<ChanEvent*>::const_iterator mtasListIt = mtasList.begin();
		actualTime = (*mtasListIt)->GetTime() * pixie::clockInSeconds;
		if (firstTime == -1.) //COMPARING EQUALITY FOR FLOATS DOESN'T ALWAYS WORK RIGHT
			firstTime = actualTime;
			cycleTime = actualTime - measureOnTime;
	}

	if(logiSummary->GetMult() > 0){//I have at leat one element in mtasList
		vector<ChanEvent*>::const_iterator logiListIt = logiList.begin();
		actualLogiTime = (*logiListIt)->GetTime() * pixie::clockInSeconds;
		cycleLogiTime = actualLogiTime - measureOnTime;
	}

	//ADDED BY THOMAS RULAND, FOR THE RU105 MEASUREMENT
/*	isBkgOn = false;
	isLightPulserOn = false;
	isTapeMoveOn = false;
	if( (actualTime - firstTime)/60 > 102 ){
		isMeasureOn = true;
		isIrradOn = false;
	}else{
		isMeasureOn = false;
		isIrradOn = true;
	}
	cycleTime = actualTime - firstTime;
*/	vector<double> starttimes = {250*60,1363*60};
	vector<double> endtimes = {370*60,1483*60};


	SetCycleState();   //sets booleans to let us know where we are in the cycle
	//options are isTapeMoveOn, isMeasureOn, isBkgOn, isLightPulserOn, isIrradOn, and cycleNumber

        //Spectrum number convention
        //0- all mtas, 1 - Central, 2 - Inner, 3 - Middle, 4 - Outer
        vector <double> totalMtasEnergy (5,-1);
        // 0-5 Central, 6-11 Inner, 12-17 Middle, 18-23 Outer
	vector <double> sumFrontBackEnergy(24,0);
	//int nrOfCentralPMT = 0; //this is easily confused with nrOfCentralPMTs set in FillMtasMap;
	double theSmallestCEnergy = 60000;


	for(map<string, struct MtasData>::const_iterator mtasMapIt = mtasMap.begin(); mtasMapIt != mtasMap.end(); mtasMapIt++){
		double signalEnergy = (*mtasMapIt).second.energy;
		double location = (*mtasMapIt).second.location;
		double time = (*mtasMapIt).second.time * pixie::clockInSeconds;

		if((*mtasMapIt).first[0] == 'C'){
			totalMtasEnergy.at(1) += signalEnergy/12.0;
			totalMtasEnergy.at(0) += signalEnergy/12.0;
			isCenter = true;	
         		//nrOfCentralPMT ++;//redundant?
			//if(theSmallestCEnergy > signalEnergy) //ORIGINAL, THIS LOGIC SEEMS BACKWARDS CHANGED TO MORE STRAIGHT FORWARD WAY -TR 1/27/2020
			if( signalEnergy < theSmallestCEnergy )
				theSmallestCEnergy = signalEnergy;
		}else if((*mtasMapIt).first[0] == 'I'){
			totalMtasEnergy.at(2) += signalEnergy/2.0;
			totalMtasEnergy.at(0) += signalEnergy/2.0;		
			isInner = true;
          		if (earliestIMOTime > 0 && time < earliestIMOTime)
            			earliestIMOTime = time;
		}else if((*mtasMapIt).first[0] == 'M'){
			totalMtasEnergy.at(3) += signalEnergy/2.0;
			totalMtasEnergy.at(0) += signalEnergy/2.0;		
			isMiddle = true;
			if (earliestIMOTime > 0 && time < earliestIMOTime)
				earliestIMOTime = time;
		}else if((*mtasMapIt).first[0] == 'O'){
			totalMtasEnergy.at(4) += signalEnergy/2.0;
			totalMtasEnergy.at(0) += signalEnergy/2.0;
			isOuter = true;
			if (earliestIMOTime > 0 && time < earliestIMOTime)
				earliestIMOTime = time;
		}
					 
		//F+B
		int moduleIndex = (location -1)/2;
		if(moduleIndex > 24){
			cout<<"Warning: detector "<<(*mtasMapIt).first<<" location > 48"<<endl;
			continue;
		}
		
		if(sumFrontBackEnergy.at(moduleIndex) == 0){		//it's the first signal from this hexagon module
			sumFrontBackEnergy.at(moduleIndex) = -1* signalEnergy/2.;
		}else if(sumFrontBackEnergy.at(moduleIndex) < 0){	//second signal from this hexagon module
			sumFrontBackEnergy.at(moduleIndex) = -1*sumFrontBackEnergy.at(moduleIndex) + signalEnergy/2.;
		}else{							//sumFrontBackEnergy.at(moduleIndex) > 0 - 3 or more signals in one event
			cout<<"Warning: detector "<<(*mtasMapIt).first<<" has 3 or more signals"<<endl;
		}
		
	}

	SetIfOnlyRingBool(); //Sets booleans for ring gating. Options are 
			     //isCenter, isInner, isMiddle, isOuter, also 
			     //isCenterOnly, isInnerOnly, isMiddleOnly, isOuterOnly, isAll
			     //Usefull options should be no more than two of these.

	//Added by Goetz for March 2015 experiment
	//this is set up to work with any generic type detector, simply add in correct subtype 
	// reference module is ref subtype
	
	/*vector <double> refModule(1,-1);
	for(map<string, struct MtasData>::const_iterator refmodMapIt = refmodMap.begin(); refmodMapIt != refmodMap.end(); refmodMapIt++)
	  {
	    double signalEnergy = (*refmodMapIt).second.energy;
	    // double location = (*refmodMapIt).second.location;

	   if((*refmodMapIt).first == "ref")
	      {
		refModule.at(0) += signalEnergy;
	      }
	    
	  }*/
	//This is likely over kill. Maybe it doesn't matter. Energy set in first loop for MapFill.

	if(isMeasureOn && !isLightPulserOn && !isTapeMoveOn){
		//silicon spectras
		plot(MTAS_POSITION_ENERGY+500, siliMap.size());
		int siliconNumber;
		//double siliconEnergy;
		for(map<string, struct MtasData>::const_iterator siliMapIt = siliMap.begin(); siliMapIt != siliMap.end(); siliMapIt++){
			siliconNumber = 0;
       			if((*siliMapIt).first[2] == 'B')//bottom - channel from 9 to 15
				siliconNumber = 8;
			siliconNumber += (*siliMapIt).first[1]-48;//48 - position of '0' character in Ascii Table
			plot(MTAS_POSITION_ENERGY+510, siliconNumber);
			plot(MTAS_POSITION_ENERGY+271, siliconNumber, cycleNumber);// Jan 03 2011
		}
	}
  
	//logic signals
	if (logicSignalsValue > 0){
	    incplot(MTAS_POSITION_ENERGY+272, cycleLogiTime, cycleNumber, logicSignalsValue);
	}

	// K. C. Goetz for March 2015 Experiment: reference crystal vs time in 1 minute intervals
	// Histogram # 5100
	//plot(MTAS_REFERENCE_EVO, refmodEnergy, (actualTime - firstTime)/60); //not using refModule.at(0) NTB
  
	if(nrOfCentralPMTs != 12 ) //redundant?|| nrOfCentralPMTs == 0)
        	totalMtasEnergy.at(1) =-1;

	//Plotting begins here. 
	if(nrOfCentralPMTs < 15)
        	plot(MTAS_POSITION_ENERGY+274, nrOfCentralPMTs);
        plot(MTAS_POSITION_ENERGY+275, totalMtasEnergy.at(1) / 10.0, nrOfCentralPMTs);
	plot(MTAS_POSITION_ENERGY+276, theSmallestCEnergy / 10.0, nrOfCentralPMTs);
	for(unsigned int i=0; i<sumFrontBackEnergy.size(); i++){
		if(sumFrontBackEnergy.at(i) < 0){//it was only one signal in the hexagon module
			sumFrontBackEnergy.at(i)=-1;  
		}
	}

	//Background  
	if(isMeasureOn && isBkgOn && !isLightPulserOn && !isTapeMoveOn){
		//3201 - 3241, no B-gated and 3301 - 3341, B-gated
		for(int i=0; i<5; i++){
			//cout<<"Total Mtas Energy: "<<totalMtasEnergy.at(i)<<endl;
			plot(MTAS_POSITION_ENERGY+201+i*10, totalMtasEnergy.at(i));
			if(isBetaSignal)
				plot(MTAS_POSITION_ENERGY+301+i*10, totalMtasEnergy.at(i));
		}
	}	

	//Timediff related plots 
	//Time difference between two beta events
	//if(isMeasureOn && isBetaSignal && !isBkgOn && !isLightPulserOn && !isTapeMoveOn){
	//bool was1stBetaevent = false;
	//bool is2ndBetaevent = false;
	if(isMeasureOn && isBetaSignal && !isLightPulserOn && !isTapeMoveOn){
		for(int j=0; j<totalMtasEnergy.size(); j++){
			plot(MTAS_POSITION_ENERGY+280+j, totalMtasEnergy.at(j));
		}
		if (time1stBetaevent != 0.0 && time2ndBetaevent == 0.0) {
			time2ndBetaevent = actualTime;
		}else if (time1stBetaevent == 0.0 && time2ndBetaevent == 0.0){
			time1stBetaevent = actualTime;
		}else if(time1stBetaevent != 0.0 && time2ndBetaevent != 0.0 ){
			double timeDiffBB = (time2ndBetaevent - time1stBetaevent)*100000000;
			plot(MTAS_POSITION_ENERGY+700, timeDiffBB);
			if (timeDiffBB > 300.0 && timeDiffBB < 750.0)
			{
				for(int k=0; k<totalMtasEnergy.size(); k++){
					plot(MTAS_POSITION_ENERGY+285+k, totalMtasEnergy.at(k));
				}
			}
			time1stBetaevent = time2ndBetaevent;
			time2ndBetaevent = actualTime;
		}
 
	}

	//Time difference between two events with first beta and second not beta
	//bool is1stBetaevent = false;
	bool is2ndNotBetaevent = false;
	if(isMeasureOn && !isLightPulserOn && !isTapeMoveOn){
		if(isBetaSignal && actualTime != -1.0){
			previousTime1stBeta = actualTime;
			//is1stBetaevent = true;
			for(int j=0; j<totalMtasEnergy.size(); j++){
				plot(MTAS_POSITION_ENERGY+250+j, totalMtasEnergy.at(j));
			}
		}else if (!isBetaSignal && previousTime1stBeta != 0.0 && actualTime != -1.0){
			previousTime2nd_notBeta = actualTime;// Update previousTime with the current actualTime previousTimeBetaBoth
			is2ndNotBetaevent = true;
		}
		if (previousTime1stBeta != 0.0 && previousTime2nd_notBeta != 0.0 /*&& previousTimeBeta1st != -1.0 && previousTime2nd != -1.0*/) {
			double timeDiffBnotB = (previousTime2nd_notBeta - previousTime1stBeta)*100000000;
			// std::cout<<" time diff:"<<timeDiffBnotB<<std::endl;
			plot(MTAS_POSITION_ENERGY+703, timeDiffBnotB);	
			if (is2ndNotBetaevent && timeDiffBnotB > 300.0 && timeDiffBnotB < 550.0)
			{
				for(int k=0; k<totalMtasEnergy.size(); k++){
					plot(MTAS_POSITION_ENERGY+255+k, totalMtasEnergy.at(k));
				}
			}		
		}	
	}
	//older version of the above
	// if(isMeasureOn && !isLightPulserOn && !isTapeMoveOn){
	// 	if(isBetaSignal && actualTime != -1.0){
	// 		previousTime1stBeta = actualTime;
	// 	}else if (previousTime1stBeta != 0.0 && !isBetaSignal && actualTime != -1.0){
	// 		previousTime2nd_notBeta = actualTime;// Update previousTime with the current actualTime previousTimeBetaBoth
	// 	}
	// 	if (previousTime1stBeta != 0.0 && previousTime2nd_notBeta != 0.0) {
	// 		double timeDiffBnotB = (previousTime2nd_notBeta - previousTime1stBeta)*100000000;
	// 		plot(MTAS_POSITION_ENERGY+703, timeDiffBnotB);
	// 		plot(MTAS_POSITION_ENERGY+714, totalMtasEnergy.at(0) / 10.0, timeDiffBnotB); //"MTAS / 10  vs time diff (1st beta)");
	// 		plot(MTAS_POSITION_ENERGY+715, totalMtasEnergy.at(1) / 10.0, timeDiffBnotB); //"MTAS C / 10  vs time diff (1st beta)");	
	// 		if (isBetaSignal && timeDiffBnotB > 300.0 && timeDiffBnotB < 600.0) {
	// 			plot(MTAS_POSITION_ENERGY+711, totalMtasEnergy.at(0) / 10.0, totalMtasEnergy.at(1) / 10.0);
	// 			for(int i=6; i<24; i++){					
	// 				plot(MTAS_POSITION_ENERGY+710, totalMtasEnergy.at(0) / 10.0, sumFrontBackEnergy.at(i) / 10.0);
	// 			}
	// 			plot(MTAS_POSITION_ENERGY+712, totalMtasEnergy.at(0) / 10.0, cycleTime);
	// 			plot(MTAS_POSITION_ENERGY+713, totalMtasEnergy.at(1) / 10.0, cycleTime); //Central vs. cycle time
	// 		}		
	// 	}
	// 	// if(isBetaSignal){
	// 	// 	previousTime1stBeta = actualTime;
	// 	// }else if (previousTime1stBeta != 0.0 && !isBetaSignal){
	// 	// 	previousTime2nd_notBeta = actualTime;// Update previousTime with the current actualTime previousTimeBetaBoth
	// 	// }	
	// }
	// //Time difference between two events with first beta only
	// if(isMeasureOn && !isLightPulserOn && !isTapeMoveOn){
	// 	if (previousTimeBeta1st != 0.0 && previousTime2nd != 0.0 && previousTimeBeta1st != -1.0 && previousTime2nd != -1.0) {
	// 		double timeDiffBA = (previousTime2nd - previousTimeBeta1st)*100000000;
	// 		//std::cout<<"  previousTime2nd  "<<previousTime2nd<<"  previousTimeBeta1st  "<<previousTimeBeta1st<<"  2-1 *10^8  "<<timeDiffBA<<std::endl;
	// 		plot(MTAS_POSITION_ENERGY+701, timeDiffBA);
	// 		plot(MTAS_POSITION_ENERGY+708, totalMtasEnergy.at(0) / 10.0, timeDiffBA); //"MTAS / 10  vs time diff (1st beta)");
	// 		plot(MTAS_POSITION_ENERGY+709, totalMtasEnergy.at(1) / 10.0, timeDiffBA); //"MTAS C / 10  vs time diff (1st beta)");	
	// 		if (timeDiffBA > 310.0 && timeDiffBA < 530.0) {
	// 			plot(MTAS_POSITION_ENERGY+705, totalMtasEnergy.at(0) / 10.0, totalMtasEnergy.at(1) / 10.0);
	// 			for(int i=6; i<24; i++){					
	// 				plot(MTAS_POSITION_ENERGY+704, totalMtasEnergy.at(0) / 10.0, sumFrontBackEnergy.at(i) / 10.0);
	// 			}
	// 			plot(MTAS_POSITION_ENERGY+706, totalMtasEnergy.at(0) / 10.0, cycleTime);
	// 			plot(MTAS_POSITION_ENERGY+707, totalMtasEnergy.at(1) / 10.0, cycleTime); //Central vs. cycle time
	// 		}		
	// 	}
	// 	if(isBetaSignal){
	// 		previousTimeBeta1st = actualTime;
	// 	}else if (previousTimeBeta1st != 0.0){
	// 		previousTime2nd = actualTime;// Update previousTime with the current actualTime previousTimeBetaBoth
	// 	}	
	// }	

	// //Time difference between any two events
	// if(isMeasureOn && !isBkgOn && !isLightPulserOn && !isTapeMoveOn){
	// 	if (previousTimeany != 0.0) {
	// 		double timeDiffAA = (actualTime - previousTimeany)*100000000;
	// 		plot(MTAS_POSITION_ENERGY+702, timeDiffAA);
	// 	}
	// 	previousTimeany = actualTime; // Update previousTime with the current actualTime previousTimeBetaBoth
	// }		

	//Light Pulser  
	if(isLightPulserOn){
		//3202 - 3242, no B-gated
		for(int i=0; i<5; i++)
			plot(MTAS_POSITION_ENERGY+202+i*10, totalMtasEnergy.at(i));
		
		//3101 - 3147 odd, no B-gated
		for(unsigned int i=0; i<sumFrontBackEnergy.size(); i++)
			plot(MTAS_POSITION_ENERGY+100+3*i+2, sumFrontBackEnergy.at(i));

		//for (int i=0; i<5; i++)
		//  	plot(MTAS_LIGHTPULSER_EVO+i, totalMtasEnergy.at(i), cycleTime);
	}  
	
	for(map<string, struct MtasData>::const_iterator geMapIt = geMap.begin(); geMapIt != geMap.end(); geMapIt++){
		    plot(MTAS_POSITION_ENERGY+416, (*geMapIt).second.calEnergy, cycleTime);
		    if (implantMult == 2)
     		    	plot(MTAS_POSITION_ENERGY+417, (*geMapIt).second.calEnergy);
	}
		
	//Irradiation  
	if(isIrradOn && !isBkgOn){
		//3203 - 3243, no B-gated and 3303 - 3343, B-gated
		for(int i=0; i<5; i++){
			plot(MTAS_POSITION_ENERGY+203+i*10, totalMtasEnergy.at(i));
			if(isBetaSignal)
				plot(MTAS_POSITION_ENERGY+303+i*10, totalMtasEnergy.at(i));
		}
	} 
  
	//Irradiation and Bkg 
	if(isIrradOn && isBkgOn && isBetaSignal){
		//3304 - 3344, B-gated
		for(int i=0; i<5; i++)
			plot(MTAS_POSITION_ENERGY+304+i*10, totalMtasEnergy.at(i));
	}   
  
	//3260 mtas vs cycle time, no logic conditions
	plot(MTAS_POSITION_ENERGY+260, totalMtasEnergy.at(0) / 10.0, cycleTime);
	plot(MTAS_POSITION_ENERGY+262, totalMtasEnergy.at(1) / 10.0, cycleTime); //Central vs. cycle time
	plot(MTAS_POSITION_ENERGY+270, totalMtasEnergy.at(0) / 10.0, cycleNumber);// Jan 03 2011
		
	if(isBetaSignal){
        	plot(MTAS_POSITION_ENERGY+273, actualTime - firstTime);
    	}
        
	for(int i=0; i<5; i++){
        	plot(MTAS_EVO_NOLOGIC+i, totalMtasEnergy.at(i), (actualTime - firstTime)/60);
		if( isBetaSignal )
        		plot(MTAS_EVO_NOLOGIC+i+10, totalMtasEnergy.at(i), (actualTime - firstTime)/60);
	}

	//"Regular" measurement 
	if(isMeasureOn && !isBkgOn && !isLightPulserOn && !isTapeMoveOn){
		//3200 - 3240, no B-gated and 3300 - 3340, B-gated
		for(int i=0; i<5; i++){
			plot(MTAS_POSITION_ENERGY+200+i*10, totalMtasEnergy.at(i));
			if(isBetaSignal) {
				plot(MTAS_POSITION_ENERGY+300+i*10, totalMtasEnergy.at(i));
            //if (i==1) {
            //see 3276 plot(MTAS_POSITION_ENERGY+306, (actualTime - firstTime));//NTB 12/15/15
            //				}
			}    
		}
		
        
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
		//3100 - 3146 even, B-gated
		for(int i=0; i<24; i++)
			plot(MTAS_POSITION_ENERGY+100+3*i, sumFrontBackEnergy.at(i));

		if(isBetaSignal) {
			//ADDED BY THOMAS RULAND 3600 AND 3601
			for( auto ii = 0; ii < starttimes.size(); ii++){
				if( cycleTime >= starttimes.at(ii) and cycleTime <= endtimes.at(ii) )
					plot(MTAS_POSITION_ENERGY+600+ii, totalMtasEnergy.at(0));
			}
			for(int i=0; i<5; i++)
        			plot(MTAS_EVO_NOLOGIC+i+20, totalMtasEnergy.at(i), cycleTime);

			//3100 - 3146 even, B-gated
			for(int i=0; i<24; i++)
				plot(MTAS_POSITION_ENERGY+100+3*i+1, sumFrontBackEnergy.at(i));
			
			//3350 mtas tot vs I, M, O, B-gated
			plot(MTAS_POSITION_ENERGY+351, totalMtasEnergy.at(0) / 10.0, totalMtasEnergy.at(1) / 10.0);
			for(int i=6; i<24; i++)					
				plot(MTAS_POSITION_ENERGY+350, totalMtasEnergy.at(0) / 10.0, sumFrontBackEnergy.at(i) / 10.0);
			//3355 Gamma-Gamma matrix, B-gated I,M,O
			for(unsigned int i=6; i<sumFrontBackEnergy.size(); i++){
				for(unsigned int j=6; j<sumFrontBackEnergy.size(); j++)
					if(i != j){
						 plot(MTAS_POSITION_ENERGY+355, sumFrontBackEnergy.at(i) / 10.0, sumFrontBackEnergy.at(j) / 10.0);
						 plot(MTAS_POSITION_ENERGY+355, sumFrontBackEnergy.at(j) / 10.0, sumFrontBackEnergy.at(i) / 10.0);
					}
			}
			plot(MTAS_POSITION_ENERGY+360, totalMtasEnergy.at(0) / 10.0, cycleTime);
			plot(MTAS_POSITION_ENERGY+361, totalMtasEnergy.at(0) / 10.0, cycleTime * 10.0);
			plot(MTAS_POSITION_ENERGY+362, totalMtasEnergy.at(1) / 10.0, cycleTime);//C vs Time (s
			plot(MTAS_POSITION_ENERGY+363, totalMtasEnergy.at(1) / 10.0, cycleTime * 10.0);//C vs. Time (100ms)
			plot(MTAS_POSITION_ENERGY+364, totalMtasEnergy.at(1) / 10.0, cycleTime / 60.0 );//C vs Time (min)
			double dt_beta_gamma = (actualTime - betaTime) * 1.0e8;
			double dt_beta_imo = (earliestIMOTime - betaTime) * 1.0e8;
			double dt_shift = 100.0;
			plot(MTAS_POSITION_ENERGY+365, totalMtasEnergy.at(0) / 10.0, dt_beta_gamma + dt_shift);
			plot(MTAS_POSITION_ENERGY+366, totalMtasEnergy.at(0) / 10.0, dt_beta_imo + dt_shift);
			// These values are read from plots above - they consitite a neutron gate
			double lowNeutronE = 6500.0;
			double highNeutronE = 8000.0;
			// Notice time gate is read from dt_beta_imo plot (with shift!) and is very restrictive
			double lowNeutronT = 40.0;
			double highNeutronT = 100.0;
			
			if (totalMtasEnergy.at(0) > lowNeutronE && totalMtasEnergy.at(0) < highNeutronE && dt_beta_imo + dt_shift > lowNeutronT && dt_beta_imo + dt_shift < highNeutronT) 
				plot(MTAS_POSITION_ENERGY+400, totalMtasEnergy.at(1));
            		// These neutron gates are based on neutron capture in I M O rings
			double sumMO = totalMtasEnergy.at(3) + totalMtasEnergy.at(4);
			double sumIMO = sumMO + totalMtasEnergy.at(2);

			if (sumIMO > lowNeutronE && sumIMO < highNeutronE){
				plot(MTAS_POSITION_ENERGY+401, totalMtasEnergy.at(0));
				plot(MTAS_POSITION_ENERGY+402, totalMtasEnergy.at(1));
			}
			if (sumMO > lowNeutronE && sumMO < highNeutronE){
				plot(MTAS_POSITION_ENERGY+403, totalMtasEnergy.at(0));
				plot(MTAS_POSITION_ENERGY+404, totalMtasEnergy.at(1));
			}

			plot(MTAS_POSITION_ENERGY+456, maxSiliconSignal );
			if( maxSiliconSignal > 700.0 && maxSiliconSignal < 3200.0 )// greater than 2505 level feeding
			//if( dt_beta_gamma < -3.0 && maxSiliconSignal > 1000.0)// good for low energy beta cuts.
			{//
				plot(MTAS_POSITION_ENERGY+457, totalMtasEnergy.at(0) );//will have comptons
				const double E_IMO = totalMtasEnergy.at(2)+totalMtasEnergy.at(3)+totalMtasEnergy.at(4);
				const bool noCenter = totalMtasEnergy.at(1) < 1.0;
				// This has few comptons due to silicon timing walk for low energy (< 320 keV) betas
				// But also has few high energy betas from 1332 level feeding
				if( dt_beta_gamma < -10.0 ) plot(MTAS_POSITION_ENERGY+458, totalMtasEnergy.at(0) );//

				if( maxSiliconSignal > 2400.0 ) plot(MTAS_POSITION_ENERGY+459, totalMtasEnergy.at(0) );//
				if( maxSiliconSignal > 2500.0 ) plot(MTAS_POSITION_ENERGY+460, totalMtasEnergy.at(0) );//
				if( maxSiliconSignal > 2600.0 ) plot(MTAS_POSITION_ENERGY+461, totalMtasEnergy.at(0) );//
				if( maxSiliconSignal > 2700.0 ) plot(MTAS_POSITION_ENERGY+462, totalMtasEnergy.at(0) );//
				if( maxSiliconSignal > 2800.0 ) plot(MTAS_POSITION_ENERGY+463, totalMtasEnergy.at(0) );//
				if( maxSiliconSignal > 2900.0 ) plot(MTAS_POSITION_ENERGY+464, totalMtasEnergy.at(0) );//
				if( maxSiliconSignal > 3000.0 ) plot(MTAS_POSITION_ENERGY+465, totalMtasEnergy.at(0) );//

				if( maxSiliconSignal > 2100.0 && noCenter ) plot(MTAS_POSITION_ENERGY+470, E_IMO );//
				if( maxSiliconSignal > 2200.0 && noCenter ) plot(MTAS_POSITION_ENERGY+471, E_IMO );//
				if( maxSiliconSignal > 2300.0 && noCenter ) plot(MTAS_POSITION_ENERGY+472, E_IMO );//
				if( maxSiliconSignal > 2400.0 && noCenter ) plot(MTAS_POSITION_ENERGY+473, E_IMO );//
				if( maxSiliconSignal > 2500.0 && noCenter ) plot(MTAS_POSITION_ENERGY+474, E_IMO );//
				if( maxSiliconSignal > 2600.0 && noCenter ) plot(MTAS_POSITION_ENERGY+475, E_IMO );//
				if( maxSiliconSignal > 2700.0 && noCenter ) plot(MTAS_POSITION_ENERGY+476, E_IMO );//

				if( maxSiliconSignal > 2000.0 && noCenter && dt_beta_gamma < 10.0 ) plot(MTAS_POSITION_ENERGY+477, E_IMO );//
				if( maxSiliconSignal > 2000.0 && noCenter && dt_beta_gamma < 0.0 ) plot(MTAS_POSITION_ENERGY+478, E_IMO );//
			}//
		}
	}  
	EndProcess(); // update the processing time
	return true;
}

MtasProcessor::MtasData::MtasData(string type){
	detSubtype     = type;
	energy         = -1.;
	calEnergy      = -1.;
	time           = -1.;
	location       = -1.;
}

MtasProcessor::MtasData::MtasData(ChanEvent *chan){
	detSubtype = chan->GetChanID().GetSubtype();
	energy = chan->GetEnergy();
	calEnergy = chan->GetCalEnergy();
	time = chan->GetTime();
	location = chan->GetChanID().GetLocation();
}

void MtasProcessor::SetIfOnlyRingBool(){
	if(isCenter & isInner & isMiddle & isOuter) isAll = true;
	if(isCenter & !isInner & !isMiddle & !isOuter) isCenterOnly = true;
	if(!isCenter & isInner & !isMiddle & !isOuter) isInnerOnly = true;
	if(!isCenter & !isInner & isMiddle & !isOuter) isMiddleOnly = true;
	if(!isCenter & !isInner & !isMiddle & isOuter) isOuterOnly = true;
}

void MtasProcessor::SetCycleState(){
	//Check flags and set main (static) flags

	//tape		
	if(isTapeMoveOffSignal && isTapeMoveOnSignal)
		cout<<"Error: tape movement signal and end of tape movement signal in the same event"<<endl;
	if(isTapeMoveOnSignal && isTapeMoveOn) 
		cout<<"Error: No end of tape movement signal in the last tape cicle"<<endl;	
	if(isTapeMoveOnSignal)
		isTapeMoveOn = true;
	if(isTapeMoveOffSignal)
		isTapeMoveOn = false;

	//measurement  
	if(isMeasureOffSignal && isMeasureOnSignal)
		cout<<"Error: measurement signal and no measurement signal in the same event"<<endl;
	if(isMeasureOnSignal && isMeasureOn) 
		cout<<"Error: No end of measurement signal in the last tape cicle"<<endl;
	if(isMeasureOnSignal)
		isMeasureOn = true;
	if(isMeasureOffSignal)
		isMeasureOn = false; 

	//background		
	if(isBkgOffSignal && isBkgOnSignal)
		cout<<"Error: background signal and no background signal in the same event"<<endl;
	if(isBkgOnSignal && isBkgOn) 
		cout<<"Error: No end of background signal in the last tape cicle"<<endl;
	if(isBkgOnSignal)
		isBkgOn = true;
	if(isBkgOffSignal)
		isBkgOn = false; 
		
	//light pulser	
	if(isLightPulserOffSignal && isLightPulserOnSignal)
		cout<<"Error: light pulser signal and no light pulser signal in the same event"<<endl;
	if(isLightPulserOnSignal && isLightPulserOn) 
		cout<<"Error: No end of light pulser signal in the last tape cicle"<<endl;
	if(isLightPulserOnSignal)
		isLightPulserOn = true;
	if(isLightPulserOffSignal)
		isLightPulserOn = false;		
		
	//irradiation		
	if(isIrradOffSignal && isIrradOnSignal)
		cout<<"Error: irradiation signal and no irradiation signal in the same event"<<endl;	
	if(isIrradOnSignal && isIrradOn) 
		cout<<"Error: No end of irradiation signal in the last tape cicle"<<endl;
	if(isIrradOnSignal)
		isIrradOn = true;
	if(isIrradOffSignal)
		isIrradOn = false;
}

void MtasProcessor::FillMtasMap(){
	mtasMap.clear();
	maxLocation =0; 
	nrOfCentralPMTs = 0;   
	for(vector<ChanEvent*>::const_iterator mtasListIt = mtasList.begin(); mtasListIt != mtasList.end(); mtasListIt++){
		string subtype = (*mtasListIt)->GetChanID().GetSubtype();
		if(subtype[0] =='C')
		    nrOfCentralPMTs ++;
		if (mtasMap.count(subtype)>0){
			cout<<"Error: Detector "<<subtype<<" has "<< mtasMap.count(subtype)+1<<" signals in one event"<<endl;
			continue;//should I skip such events?
		}
			
		if ((*mtasListIt)->GetEnergy() == 0 || (*mtasListIt)->GetEnergy() > 30000)
				continue;
				
		mtasMap.insert(make_pair(subtype,MtasData((*mtasListIt))));
		
		if(maxLocation < (*mtasListIt)->GetChanID().GetLocation())
			maxLocation = (*mtasListIt)->GetChanID().GetLocation();
	}
	//init bools to false
	isCenter = false;
	isInner = false;
	isMiddle = false;
	isOuter = false;
	isCenterOnly = false;
	isInnerOnly = false;
	isMiddleOnly = false;
	isOuterOnly = false;
	isAll = false;
}

void MtasProcessor::FillSiliMap(){
	siliMap.clear();
	isBetaSignal = false;
	betaTime = -100.0;
	maxSiliconSignal = -1.0;
	extern DetectorDriver driver;
	for(vector<ChanEvent*>::const_iterator siliListIt = siliList.begin(); siliListIt != siliList.end(); siliListIt++){
		string subtype = (*siliListIt)->GetChanID().GetSubtype();
		if (mtasMap.count(subtype)>0)
			cout<<"Error: Detector "<<subtype<<" has "<< siliMap.count(subtype)+1<<" signals in one event"<<endl;
		
		Calibration cal = driver.cal.at((*siliListIt)->GetID());
		if ((*siliListIt)->GetEnergy() < 200 || (*siliListIt)->GetEnergy() > 30000) {
//		if ((*siliListIt)-> GetEnergy() < cal.GetMinThreshold() || (*siliListIt)->GetEnergy() > 30000) Oct '15 use hard coded threshold for online.
			continue;
        	}		
		siliMap.insert(make_pair(subtype,MtasData((*siliListIt))));
	}
	
	if(siliMap.size()>0){
            maxSiliconSignal = siliSummary->GetMaxEvent()->GetCalEnergy();//Jan 03 2011 Ola K
            if(maxSiliconSignal > 2.0) {
                isBetaSignal = true;
                betaTime = siliSummary->GetMaxEvent()->GetTime() * pixie::clockInSeconds;
            }
	}
}

void MtasProcessor::FillRefModMapAndEnergy(){
	// added by Goetz for March 2015 experiment allows for any detector type generic to be imported into mtas processor (subtype is independent)
	refmodMap.clear(); 
        refmodEnergy = 0.;
	for(vector<ChanEvent*>::const_iterator refmodListIt = refmodList.begin(); refmodListIt != refmodList.end(); refmodListIt++){
		string subtype = (*refmodListIt)->GetChanID().GetSubtype();
		if ((*refmodListIt)->GetEnergy() == 0 || (*refmodListIt)->GetEnergy() > 30000)
			continue;
		refmodMap.insert(make_pair(subtype,MtasData((*refmodListIt))));
		refmodEnergy += (*refmodListIt)->GetEnergy();
	}
}

void MtasProcessor::FillSipmMap(){
	// added by Goetz for March 2015 experiment allows for any detector type generic to be imported into mtas processor (subtype is independent)
 	sipmMap.clear();
	implantEnergy = 0.;
	implantMult = 0.;
	for(vector<ChanEvent*>::const_iterator sipmListIt = sipmList.begin(); sipmListIt != sipmList.end(); sipmListIt++){
		string subtype = (*sipmListIt)->GetChanID().GetSubtype();
		if ((*sipmListIt)->GetEnergy() == 0 || (*sipmListIt)->GetEnergy() > 30000)
			continue;
		sipmMap.insert(make_pair(subtype,MtasData((*sipmListIt))));
		if (subtype == "implant_xa" || subtype == "implant_xb")
    			implantEnergy += (*sipmListIt)->GetEnergy();
		implantMult +=1;
	}
}

void MtasProcessor::FillGeMap(){
	// added for Nov 2019 experiment at ANL with implant detector
	//geEnergy = 0.;
	geMap.clear();
	for(vector<ChanEvent*>::const_iterator geListIt = geList.begin(); geListIt != geList.end(); geListIt++){
		string subtype = (*geListIt)->GetChanID().GetSubtype();
		if ((*geListIt)->GetEnergy() == 0 || (*geListIt)->GetEnergy() > 30000)
			continue;
		geMap.insert(make_pair(subtype,MtasData((*geListIt))));
		// geEnergy += (*geListIt)->GetEnergy();
	}
}

void MtasProcessor::FillLogicMap(){
	logiMap.clear();
	isTriggerOnSignal = false;
	isTapeMoveOnSignal = false;
	isTapeMoveOffSignal = false;
	isMeasureOnSignal = false;
	isMeasureOffSignal = false;	
 	isBkgOnSignal = false;
	isBkgOffSignal = false;
 	isLightPulserOnSignal = false;
	isLightPulserOffSignal = false;	
 	isIrradOnSignal = false;
	isIrradOffSignal = false;
	
	double logicTreshold = 1;	//logic threshold !!!!!!! (value?????)	
	logicSignalsValue = 0;

	for(vector<ChanEvent*>::const_iterator logiListIt = logiList.begin(); logiListIt != logiList.end(); logiListIt++){
		string subtype = (*logiListIt)->GetChanID().GetSubtype();
		if (logiMap.count(subtype)>0)
			cout<<"Error: Detector "<<subtype<<" has "<< logiMap.count(subtype)+1<<" signals in one event"<<endl;
		logiMap.insert(make_pair(subtype,MtasData(((*logiListIt)))));
		
		//set logic flags
		if((*logiListIt)->GetEnergy() > logicTreshold){
			if(subtype == "TRU") {
				isTriggerOnSignal = true;
				measureOnTime = (*logiListIt)->GetTime() * pixie::clockInSeconds;
				cycleNumber ++;
				logicSignalsValue +=1;
			}
			
			if(subtype == "IRU"){
				isIrradOnSignal = true;
				logicSignalsValue +=2;
			}
			if(subtype == "IRD"){
				isIrradOffSignal = true;
				logicSignalsValue +=4;
			}	
			if(subtype == "LPU"){
				isLightPulserOnSignal = true;
				logicSignalsValue +=8;
			}
			if(subtype == "LPD"){
				isLightPulserOffSignal = true;
				logicSignalsValue +=16;
			}
			if(subtype == "TMU"){
				isTapeMoveOnSignal = true;
				logicSignalsValue +=32;
			}
			if(subtype == "TMD"){
				isTapeMoveOffSignal = true;
				logicSignalsValue +=64;
			}
			
			if(subtype == "BGU"){
				isBkgOnSignal = true;
				logicSignalsValue +=128;
			}
			if(subtype == "BGD"){
				isBkgOffSignal = true;
				logicSignalsValue +=256;
			}
			if(subtype == "MSU"){
				isMeasureOnSignal = true;
				logicSignalsValue +=512;
			}
			if(subtype == "MSD"){
				isMeasureOffSignal = true;
				logicSignalsValue +=1024;
			}		

		}
		
	}    
}
