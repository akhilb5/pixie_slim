///@file PspmtProcessor.cpp
///@Processes information from a Position Sensitive PMT.  No Pixel work yet. 
///@author A. Keeler, S. Go, S. V. Paulauskas, NT Brewer
///@date Mar 12, 2019

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <signal.h>
#include <limits.h>
#include <cmath>

#include "damm_plotids.h"
#include "DetectorDriver.h"
#include "MtasPspmtProcessor.h"
#include "MtasProcessor.h"

using namespace std;
using namespace dammIds::mtaspspmt;

namespace dammIds {
    namespace mtaspspmt {
        //const int DD_MULTI = 4;
        const int DD_POS_DIAG = OFFSET+1;
        const int D_POS_IMPL = OFFSET+2;
        const int D_SUME_IMPL = OFFSET+3;
        const int D_SUME_DIAG = OFFSET+4; 
        const int D_SUME_IMPL10 = OFFSET+5;
        const int D_SUME_DIAG10 = OFFSET+6;        
    }
}


void MtasPspmtProcessor::DeclarePlots(void) const
{

    //cout << " all clear" << endl;
    DeclareHistogram2D(DD_POS_DIAG, SB, SB, "Diagnostic Detector Positions");
    DeclareHistogram1D(D_POS_IMPL, SD, "Implant Detector Positions");
    DeclareHistogram1D(D_SUME_IMPL, SD, "Sum Energy of Implant");
    DeclareHistogram1D(D_SUME_DIAG, SD, "Sum Energy of Diagnostic");
    DeclareHistogram1D(D_SUME_IMPL10, SD, "Sum Energy/10 of Implant");
    DeclareHistogram1D(D_SUME_DIAG10, SD, "Sum Energy/10 of Diagnostic");
    //DeclareHistogram2D(DD_MULTI,S3,S3, "Multiplicity");
}

const string MtasPspmtProcessor::defaultConfigFile="mtasPspmtConfig.txt";

MtasPspmtProcessor::MtasPspmtProcessor() : EventProcessor()
{
    name = "mtaspspmt";
    associatedTypes.insert("mtaspspmt");
}

//bool MtasPspmtProcessor::Init(DetectorDriver &driver, const string &ConfigFile)
bool MtasPspmtProcessor::Init(DetectorDriver &driver)
{
    

   //DetectorDriver driver;
   const std::string ConfigFile = defaultConfigFile;
   //Pspmt specific function to read in a config file and initialize params
    //  cout << "being init: reading pspmt config file" << endl;
    if (!EventProcessor::Init(driver))
	return false;

    ifstream in(ConfigFile.c_str());
    if (!in) {
      cout << "failed to open the filter parameter file" << endl;
      cout << "  using default values" << endl;
      return false;
    }

    while (!in.eof()) {
      if ( isdigit(in.peek()) ) {
	in >> VDtypeStr;
	in >> implantScale_  >> implantOffset_ >> implantThreshold_;
	in >> diagnosticScale_ >> diagnosticOffset_ >> diagnosticThreshold_;
	cout << "Pspmt configuration parameters are: "
	 << " Voltage Divider " << VDtypeStr <<"\n"
	 << " implant: " << implantScale_  << " " <<  implantOffset_ << " " <<  implantThreshold_ <<"\n"
	 << " diagnostic: " << diagnosticScale_ << " " << diagnosticOffset_ << " " << diagnosticThreshold_ 
         << endl;
	break;
      } else {
	// assume this is a comment
	in.ignore(1000, '\n');
      }
    }
    if (in.fail()) {
      cout << "problem reading pspmt config file" << endl;
      return false;
    }
    cout << "successful pspmt init" << endl;
    if(VDtypeStr == "1018" || VDtypeStr == "1730")
        dttype_ = implant;
    else if(VDtypeStr == "0926")
        dttype_ = diagnostic;
    else
        dttype_ = UNKNOWN;

    ThreshStr = implantThreshold_;

    return true;

}

bool MtasPspmtProcessor::Process(RawEvent &event) {
    //  cout << "processing" << endl;
    if (!EventProcessor::Process(event))
        return false;

    /*if (DetectorDriver::get()->GetSysRootOutput()) {
        PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
    }*/
    mtasPspmtMap.clear();

    //Lists declared in header
    //read in anode & dynode signals
    pspmtList = event.GetSummary("mtaspspmt")->GetList();
    //cout << " size: " << pspmtList.size() << endl; 
    //cout << " energy: " << (*pspmtList.begin())->GetEnergy() << " " << (*pspmtList.begin())->GetChanID().GetLocation() << endl; 
    FillMtasPspmtMap();

    //set up position calculation for low / high gain yso signals and ion scint
    position_mtas.first = 0., position_mtas.second = 0.;

    //initalize all the things
    double xa = 0., ya = 0., xb = 0., yb = 0.;


    double sumImpl = 0;
    double sumDiag = 0;

    for(map<string, struct PSData>::const_iterator pspmtIt = mtasPspmtMap.begin();
            pspmtIt != mtasPspmtMap.end(); pspmtIt ++)
    {

        std::string dttype = (*pspmtIt).first;
        if( dttype.substr(0,7) == "implant")
        {

           double energy = (*pspmtIt).second.energy;
           std::string tag = (*pspmtIt).second.tag;
           dttype_ = implant;
	   if (energy < implantThreshold_)
               continue;
            //parcel out position signals by tag. here tag is designated by loc in map.txt
           if (tag == "xa" && xa == 0) {
               xa = energy;
               sumImpl += energy;
           }
           if (tag == "xb" && xb == 0) {
               xb = energy;
               sumImpl += energy;
           }
	   //cout <<"here " << xa << " , " << xb  << endl;
           if (xa > 0 && xb > 0){
           
              position_mtas.first = CalculatePosition(xa, xb, ya, yb, dttype_).first;
              position_mtas.second = CalculatePosition(xa, xb, ya, yb, dttype_).second;
              plot(D_POS_IMPL, position_mtas.first * implantScale_ + implantOffset_);
              plot(D_SUME_IMPL, sumImpl);
              plot(D_SUME_IMPL10, sumImpl/10.);
           }
	}
        //cout << "sub " << dttype.substr(0,10) << endl;
        if( dttype.substr(0,10) == "diagnostic")
        {
           double energy = (*pspmtIt).second.energy;
           std::string tag = (*pspmtIt).second.tag;
           dttype_ = diagnostic;
           
	   if (energy < diagnosticThreshold_) {
               continue;
	   }
            //parcel out position signals by tag. here tag is designated by loc in map.txt
           if (tag == "xa" && xa == 0) {
               xa = energy;
               sumDiag += energy;
           }
           if (tag == "xb" && xb == 0) {
               xb = energy;
               sumDiag += energy;
           }
           if (tag == "ya" && ya == 0) {
               ya = energy;
               sumDiag += energy;
           }
           if (tag == "yb" && yb == 0) {
               yb = energy;
               sumDiag += energy;
           }
           //cout << "xa,b,ya,b" << xa << " , " << xb << " , " << ya 
	   //	<< " , " << yb << endl;
           if (xa > 0 && xb > 0 && ya > 0 && yb > 0){
	//	cout << "xa,b,ya,b" << xa << " , " << xb << " , " << ya 
	//   	<< " , " << yb << endl;
              position_mtas.first = CalculatePosition(xa, xb, ya, yb, dttype_).first;
              position_mtas.second = CalculatePosition(xa, xb, ya, yb, dttype_).second;
              plot(DD_POS_DIAG, position_mtas.first * diagnosticScale_ + diagnosticOffset_,
                   position_mtas.second * diagnosticScale_ + diagnosticOffset_);
              plot(D_SUME_DIAG, sumDiag);
              plot(D_SUME_DIAG10, sumDiag/10.);
           }
	}
    }


    //Begin plotting

    EndProcess();
    return true;

}

pair<double, double> MtasPspmtProcessor::CalculatePosition(
        double &xa, double &xb, double &ya, double &yb,
        const MtasPspmtProcessor::DTTYPES &dttype) {
    double x = 0, y = 0;
    switch(dttype){
        case implant:
            x = (xa - xb ) / (xa+xb);
            y = 0.5;
            break;
        case diagnostic:
            x = ((xa + ya)-(xb + yb))/(xa+xb+ya+yb);
            y = ((xa + xb)-(ya + yb))/(xa+xb+ya+yb);
            break;
        case UNKNOWN:
        default:
            cout << "We received a Detector Type we didn't recognize: " << dttype << endl;
    }
    return make_pair(x, y);
}


void MtasPspmtProcessor::FillMtasPspmtMap()
{
	//different from utkscan no plotting done here and all subtypes in event filled at once.

    for (vector<ChanEvent*>::const_iterator it = pspmtList.begin(); it != pspmtList.end(); it++) {
        //cout << "inserted " <<  (*it)->GetEnergy() << " tag " 
        //     << (*it)->GetChanID().GetLocation() <<  endl;
        string subtype = (*it)->GetChanID().GetSubtype();
        mtasPspmtMap.insert(make_pair(subtype,PSData((*it))));
    }
    //cout << "mapsize" << mtasPspmtMap.size() << endl;
    
    /*map<char,int> mymap;
    for (int i=0; i<5; i++)
    {
	mymap.insert ( std::pair<char,int>('a',i*2) );
    }*/
}


MtasPspmtProcessor::PSData::PSData(string type)
{
    subtype        = type;
    energy         = -1.;
    calEnergy      = -1.;
    time           = -1.;
    loc            = -1.;
    tag            = "na";
}

MtasPspmtProcessor::PSData::PSData(ChanEvent *chan)
{
    subtype	= chan->GetChanID().GetSubtype();
    energy	= chan->GetEnergy();
    calEnergy = chan->GetCalEnergy();
    time = chan->GetTime();
    loc  = chan->GetChanID().GetLocation();
    tag = subtype.substr(subtype.find("_")+1,2);
    
}
