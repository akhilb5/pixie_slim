///@file PspmtProcessor.cpp
///@Processes information from a Position Sensitive PMT.  No Pixel work yet. 
///@author A. Keeler, S. Go, S. V. Paulauskas 
///@date July 8, 2018

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
#include "PspmtProcessor.h"

using namespace std;
using namespace dammIds::pspmt;

namespace dammIds {
    namespace pspmt {
        const int DD_DYNODE_QDC = OFFSET+0;
        const int DD_POS_LOW = OFFSET+1;
        const int DD_POS_HIGH = OFFSET+2;
        const int DD_PLASTIC_EN = OFFSET+3;
        const int DD_MULTI = OFFSET+4;
        const int DD_DY_SUM_LG =OFFSET+5;
        const int DD_DY_SUM_HG =OFFSET+6;

        const int D_TRANS_EFF_YSO = OFFSET+10;
        const int DD_SEPAR_GATED_LOW = OFFSET+11;
        const int DD_DESI_GATED_LOW = OFFSET+12;
        const int D_DESI_ENERGY = OFFSET+15;
        const int D_DESI_YSO_GATED = OFFSET+16;
        const int DD_SEPAR_ENERGY = OFFSET+17;
        const int DD_SEPAR_YSO_GATED = OFFSET+18;

        const int DD_POS_ION = OFFSET+20;
        const int DD_SEPAR_GATED_ION = OFFSET+21;
        const int DD_DESI_GATED_ION = OFFSET+22;
    }
}


void PspmtProcessor::DeclarePlots(void) {
    DeclareHistogram2D(DD_DYNODE_QDC, SD, S2, "Dynode QDC- Low gain 0, High gain 1");
    DeclareHistogram2D(DD_POS_LOW, SB, SB, "Low-gain Positions");
    DeclareHistogram2D(DD_POS_HIGH, SB, SB, "High-gain Positions");
    DeclareHistogram2D(DD_PLASTIC_EN,SD,S4, "Plastic Energy, 0-3 = VETO, 5-8 = Ion Trigger");
    DeclareHistogram2D(DD_MULTI,S3,S3, "Dynode:Anode(+2) Multi Low gain 0, High gain 1");
    DeclareHistogram2D(DD_DY_SUM_LG,SA,SA,"Low Gain Dynode vs Anode Sum");
    DeclareHistogram2D(DD_DY_SUM_HG,SA,SA,"High Gain Dynode vs Anode Sum");
    DeclareHistogram1D(D_TRANS_EFF_YSO, S3, "Separator events (0) in ion scint (1), YSO (2), and veto (3)");
    DeclareHistogram2D(DD_SEPAR_GATED_LOW, SB, SB, "Separator-gated low-gain positions");
    DeclareHistogram2D(DD_DESI_GATED_LOW, SB, SB, "Silicon dE-gated low-gain positions");
    DeclareHistogram1D(D_DESI_ENERGY, SC, "Separator-gated dE-silicon events");
    DeclareHistogram1D(D_DESI_YSO_GATED, SC, "YSO-gated dE-silicon");
    DeclareHistogram2D(DD_SEPAR_ENERGY, S2, SC, "dE-silicon-gated separator events");
    DeclareHistogram2D(DD_SEPAR_YSO_GATED, S2, SC, "YSO-gated separator events");

    DeclareHistogram2D(DD_POS_ION, SB, SB, "Ion-scint positions - ungated");
    DeclareHistogram2D(DD_SEPAR_GATED_ION, SB, SB, "Ion-scint positions - separator-gated");
    DeclareHistogram2D(DD_DESI_GATED_ION, SB, SB, "Ion-scint positions - silicon dE-gated");
}

const string PspmtProcessor::defaultConfigFile="pspmtConfig.txt";

PspmtProcessor::PspmtProcessor() : EventProcessor()
{
    name = "pspmt";
    associatedTypes.insert("pspmt");
}

bool PspmtProcessor::Init(DetectorDriver &driver, const string &ConfigFile)
{
    //Pspmt specific function to read in a config file and initialize params

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
	in >> positionScale_  >> positionOffset_ >> threshold_;
	in >> front_positionScale_ >> front_positionOffset_ >> front_threshold_;
	in >> rotation;
	cout << "Pspmt configuration parameters are: "
	 << " Voltage Divider " << VDtypeStr
	 << " YSO: " << positionScale_ << " " << positionOffset_ << " " << threshold_
	 << " Front: " << front_positionScale_ << " " << front_positionOffset_ << " " << front_threshold_
	 << "Rotation " <<  rotation
         << endl;
	// scale thresholds by the length of integration (i.e. rise time)
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

    if(VDtypeStr == "1018" || VDtypeStr == "1730")
        vdtype_ = corners;
    else if(VDtypeStr == "0926")
        vdtype_ = sides;
    else
        vdtype_ = UNKNOWN;

    ThreshStr = threshold_;
    rotation_ = rotation* 3.1415926 / 180.; 

    return true;

}

bool PspmtProcessor::Process(RawEvent &event) {

    if (!EventProcessor::Process(event))
        return false;

    /*if (DetectorDriver::get()->GetSysRootOutput()) {
        PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
    }*/
    pspmtMap.clear();

    //Lists declared in header
    //read in anode & dynode signals
    hiDynode = event.GetSummary("pspmt:dynode_high")->GetList();
    lowDynode = event.GetSummary("pspmt:dynode_low")->GetList();
    hiAnode = event.GetSummary("pspmt:anode_high")->GetList();
    lowAnode = event.GetSummary("pspmt:anode_low")->GetList();
    //read in other signals
    veto = event.GetSummary("pspmt:veto")->GetList();
    ionTrig = event.GetSummary("pspmt:ion")->GetList();
    desi = event.GetSummary("pspmt:desi")->GetList();
    separatorScint = event.GetSummary("pspmt:f11")->GetList();

    FillPspmtMap();

    //set up position calculation for low / high gain yso signals and ion scint
    position_low.first = 0, position_low.second = 0;
    position_high.first = 0, position_high.second = 0;
//initalized all the things
    double energy = 0;
    double xa_l = 0, ya_l = 0, xb_l = 0, yb_l = 0;
    double xa_h = 0, ya_h = 0, xb_h = 0, yb_h = 0;
    double top_l = 0, top_r = 0, bottom_l = 0, bottom_r = 0;
    bool hasPosition_low = false, hasPosition_high = false, hasPosition_ion = false, hasUpstream = false,
    hasDeSi = false, hasVeto = false;

    plot(DD_MULTI, lowDynode.size(), 0);
    plot(DD_MULTI, hiDynode.size(), 1);

    plot(DD_MULTI, lowAnode.size(), 2);
    plot(DD_MULTI, hiAnode.size(), 3);


    double lowAnodeSum = 0;
    double highAnodeSum = 0;

    int numOfVetoChans = veto.size();

    for(map<string, struct PSData>::const_iterator pspmtIt = pspmtMap.begin();
            pspmtIt != pspmtMap.end(); pspmtIt ++)
    {
        if( (*pspmtIt).first == "anode_low")
        {
           double energy = (*pspmtIt).second.calEnergy;
           double tag = (*pspmtIt).second.tag;

	   if (energy < threshold_)
               continue;
            //parcel out position signals by tag. here tag is designated by loc in map.txt
            if (tag == 0 && xa_l == 0) {
                xa_l = energy;
                lowAnodeSum += energy;
            }
            if (tag == 1 && xb_l == 0) {
                xb_l = energy;
                lowAnodeSum += energy;
            }
            if (tag == 2 && ya_l == 0) {
                ya_l = energy;
                lowAnodeSum += energy;
            }
            if (tag == 3 && yb_l == 0) {
                yb_l = energy;
                lowAnodeSum += energy;
            }
	}

        if( (*pspmtIt).first == "anode_high")
        {
           double energy = (*pspmtIt).second.calEnergy;
           double tag = (*pspmtIt).second.tag;

	   if (energy < threshold_ || energy > 63000)
               continue;
            //parcel out position signals by tag. here tag is designated by loc in map.txt
            if (tag == 0 && xa_h == 0) {
                xa_h = energy;
                highAnodeSum += energy;
            }
            if (tag == 1 && xb_h == 0) {
                xb_h = energy;
                highAnodeSum += energy;
            }
            if (tag == 2 && ya_h == 0) {
                ya_h = energy;
                highAnodeSum += energy;
            }
            if (tag == 3 && yb_h == 0) {
                yb_h = energy;
                highAnodeSum += energy;
            }
	}

        if( (*pspmtIt).first == "veto")
        {
           double energy = (*pspmtIt).second.calEnergy;
           double tag = (*pspmtIt).second.tag;

           plot(DD_PLASTIC_EN, energy, tag);

           if (energy > 1 && energy < 10000){
               hasVeto = true;
           }


	}

        if( (*pspmtIt).first == "ion")
        {
            double energy = (*pspmtIt).second.calEnergy;
            if (energy < 10)
                continue;

            double tag = (*pspmtIt).second.tag;
           //plot(DD_PLASTIC_EN, energy, loc + numOfVetoChans + 1);
            //parcel out position signals by tag. here tag is designated by loc in map.txt
            if (tag == 0 && top_l == 0) {
                top_l = energy;
            }
            if (tag == 1 && top_r == 0) {
                top_r = energy;
            }
            if (tag == 2 && bottom_l == 0) {
                bottom_l = energy;
            }
            if (tag == 3 && bottom_r == 0) {
                bottom_r = energy;
            }
	}

        if( (*pspmtIt).first == "f11")
        {
           double energy = (*pspmtIt).second.calEnergy;
           if (energy > 1 && energy < 10000){
               hasUpstream = true;
           }
	}
        if( (*pspmtIt).first == "desi")
        {
           double energy = (*pspmtIt).second.calEnergy;
           if (energy > 1 && energy < 10000){
               hasDeSi = true;
           }
        }
    }

    //Begin plotting

    if (xa_l > 0 && xb_l > 0 && ya_l > 0 && yb_l > 0){
        hasPosition_low = true;
        position_low.first = CalculatePosition(xa_l, xb_l, ya_l, yb_l, vdtype_, rotation_).first;
        position_low.second  = CalculatePosition(xa_l, xb_l, ya_l, yb_l, vdtype_, rotation_).second;
        plot(DD_POS_LOW, position_low.first * positionScale_ + positionOffset_,
             position_low.second * positionScale_ + positionOffset_);
    }

    if (xa_h > 0 && xb_h > 0 && ya_h > 0 && yb_h > 0){
        hasPosition_high = true;
        position_high.first = CalculatePosition(xa_h, xb_h, ya_h, yb_h, vdtype_, rotation_).first;
        position_high.second = CalculatePosition(xa_h, xb_h, ya_h, yb_h, vdtype_, rotation_).second;
        plot(DD_POS_HIGH, position_high.first * positionScale_ + positionOffset_,
             position_high.second * positionScale_ + positionOffset_ );
    }

    if (top_l > 0 && top_r > 0 && bottom_l > 0 && bottom_r > 0){
        hasPosition_ion = true;
        position_ion.first = (top_l + bottom_l - top_r - bottom_r) / (top_l + top_r + bottom_l + bottom_r);
        position_ion.second  = (top_l + top_r - bottom_l - bottom_r) / (top_l + top_r + bottom_l + bottom_r);
        plot(DD_POS_ION, position_ion.first * front_positionScale_ + front_positionOffset_,
             position_ion.second * front_positionScale_ + front_positionOffset_);
    }
      /*   plot(DD_DYNODE_QDC, (*it)->GetTrace().GetQdc(), 0);
           plot(DD_DYNODE_QDC, (*it)->GetTrace().GetQdc(), 1);
      */
        //plot valid YSO positions and dE silicon events gated on upstream events
        //plot upstream events gated on dE silicon
        //plot transmission efficiency from upstream to YSO and veto
    
    if(hasUpstream) {
        for(vector<ChanEvent*>::const_iterator de_it = desi.begin(); de_it != desi.end(); de_it++) {
            plot(D_DESI_ENERGY,(*de_it)->GetCalEnergy());
        }
        if(hasPosition_low){
            plot(DD_SEPAR_GATED_LOW, position_low.first * positionScale_ + positionOffset_,
                    position_low.second * positionScale_ + positionOffset_);
        }
        if(hasPosition_ion){
            plot(DD_SEPAR_GATED_ION, position_ion.first * positionScale_ + positionOffset_,
                    position_ion.second * positionScale_ + positionOffset_);
        }
    }

    if(hasDeSi) {
        for (vector<ChanEvent*>::const_iterator it_sep = separatorScint.begin(); it_sep != separatorScint.end(); it_sep++) {
            if ((*it_sep)->GetChanID().GetLocation() == 0) {
                plot(DD_SEPAR_ENERGY, (*it_sep)->GetCalEnergy(), 0);
            } else if ((*it_sep)->GetChanID().GetLocation() == 1) {
                plot(DD_SEPAR_ENERGY, (*it_sep)->GetCalEnergy(), 1);
            }
        }

        if (hasPosition_low) {
            plot(DD_DESI_GATED_LOW, position_low.first * positionScale_ + positionOffset_,
                 position_low.second * positionScale_ + positionOffset_);
        }
        if (hasPosition_ion) {
            plot(DD_DESI_GATED_ION, position_ion.first * positionScale_ + positionOffset_,
                 position_ion.second * positionScale_ + positionOffset_);
        }
    }

    if(hasPosition_low) {
        for (vector<ChanEvent*>::const_iterator de_it = desi.begin(); de_it != desi.end(); de_it++) {
            plot(D_DESI_YSO_GATED, (*de_it)->GetCalEnergy());
        }
        for (vector<ChanEvent*>::const_iterator it_sep = separatorScint.begin(); it_sep != separatorScint.end(); it_sep++) {
            if ((*it_sep)->GetChanID().GetLocation() == 0) {
                plot(DD_SEPAR_YSO_GATED, (*it_sep)->GetCalEnergy(), 0);
            } else if ((*it_sep)->GetChanID().GetLocation() == 1) {
                plot(DD_SEPAR_YSO_GATED, (*it_sep)->GetCalEnergy(), 1);
            }
        }
    }

    if(hasUpstream)
        plot(D_TRANS_EFF_YSO, 0);
    if(hasUpstream && hasPosition_ion)
        plot(D_TRANS_EFF_YSO, 1);
    if(hasUpstream && hasPosition_low)
        plot(D_TRANS_EFF_YSO, 2);
    if (hasUpstream && hasVeto)
        plot(D_TRANS_EFF_YSO, 3);



/*    if (!lowDynode.empty())
        plot(DD_DY_SUM_LG,lowDynode.front()->GetCalibratedEnergy(),lowAnodeSum);

    if (!hiDynode.empty())
        plot(DD_DY_SUM_HG,hiDynode.front()->GetCalibratedEnergy(),highAnodeSum);
*/
    EndProcess();
    return true;

}

pair<double, double> PspmtProcessor::CalculatePosition(double &xa, double &xb, double &ya, double &yb, const VDTYPES &vdtype, double &rot){

    double x = 0, y = 0, x_tmp = 0, y_tmp = 0, center = 0;

    switch(vdtype){
        case corners:
            x_tmp = (0.5 * (yb + xa)) / (xa + xb + ya + yb);
            y_tmp = (0.5 * (xa + xb)) / (xa + xb + ya + yb);
            center = 0.2;
            break;
        case sides:
            x_tmp = (xa - xb) / (xa + xb);
            y_tmp = (ya - yb) / (ya + yb);
            center = 0;
            break;
        case UNKNOWN:
        default:
            cerr<<"We recieved a VD_TYPE we didn't recognize " << vdtype << endl;

    }
    x = (x_tmp - center) * cos(rot) - (y_tmp - center) * sin(rot) + center;  //rotate positions about center of image by angle rot    
    y = (x_tmp - center) * sin(rot) + (y_tmp - center) * cos(rot) + center;
    return make_pair(x, y);
}


void PspmtProcessor::FillPspmtMap()
{
	//different from utkscan no plotting done here and all subtypes in event filled at once.

    for (vector<ChanEvent*>::const_iterator it = lowDynode.begin(); it != lowDynode.end(); it++) {
        string subtype = (*it)->GetChanID().GetSubtype();
        pspmtMap.insert(make_pair(subtype,PSData((*it))));
    }
    for (vector<ChanEvent*>::const_iterator it = hiDynode.begin(); it != hiDynode.end(); it++) {
        string subtype = (*it)->GetChanID().GetSubtype();
        pspmtMap.insert(make_pair(subtype,PSData((*it))));
    }
    for (vector<ChanEvent*>::const_iterator it = lowAnode.begin(); it != lowAnode.end(); it++) {
        string subtype = (*it)->GetChanID().GetSubtype();
        pspmtMap.insert(make_pair(subtype,PSData((*it))));
    }
    for (vector<ChanEvent*>::const_iterator it = hiAnode.begin(); it != hiAnode.end(); it++) {
        string subtype = (*it)->GetChanID().GetSubtype();
        pspmtMap.insert(make_pair(subtype,PSData((*it))));
    }
    for(vector<ChanEvent*>::const_iterator it = veto.begin(); it != veto.end(); it++){
        string subtype = (*it)->GetChanID().GetSubtype();
        pspmtMap.insert(make_pair(subtype,PSData((*it))));
    }
    for(vector<ChanEvent*>::const_iterator it = ionTrig.begin(); it != ionTrig.end(); it++){
        string subtype = (*it)->GetChanID().GetSubtype();
        pspmtMap.insert(make_pair(subtype,PSData((*it))));
    }
    for(vector<ChanEvent*>::const_iterator it = separatorScint.begin(); it != separatorScint.end(); it++){
        string subtype = (*it)->GetChanID().GetSubtype();
        pspmtMap.insert(make_pair(subtype,PSData((*it))));
    }
    for (vector<ChanEvent*>::const_iterator it = desi.begin(); it != desi.end(); it++){
        string subtype = (*it)->GetChanID().GetSubtype();
        pspmtMap.insert(make_pair(subtype,PSData((*it))));
    }
}


PspmtProcessor::PSData::PSData(string type)
{
    subtype        = type;
    energy         = -1.;
    calEnergy      = -1.;
    time           = -1.;
    tag            = -1.;
}

PspmtProcessor::PSData::PSData(ChanEvent *chan)
{
    subtype	= chan->GetChanID().GetSubtype();
    energy	= chan->GetEnergy();
    calEnergy = chan->GetCalEnergy();
    time = chan->GetTime();
    tag  = chan->GetChanID().GetLocation();
}

/*PspmtProcessor::PspmtProcessor(const std::string &vd, const double &yso_scale, const unsigned int &yso_offset,
			       const double &yso_threshold, const double &front_scale,
			       const unsigned int &front_offset, const double &front_threshold, const double &rotation)
				:EventProcessor(OFFSET, RANGE, "PspmtProcessor"){

//copied from utkscan and replaced by Init. left here as a reference only.
// DO NOT UNCOMMENT.
    if(vd == "SIB064_1018" || vd == "SIB064_1730")
        vdtype_ = corners;
    else if(vd == "SIB064_0926")
        vdtype_ = sides;
    else
        vdtype_ = UNKNOWN;

    VDtypeStr = vd;
    positionScale_ = yso_scale;
    positionOffset_ = yso_offset;
    threshold_ = yso_threshold;
    front_positionScale_ = front_scale;
    front_positionOffset_ = front_offset;
    front_threshold_ = front_threshold;
    ThreshStr = yso_threshold;
    rotation_ = rotation * 3.1415926 / 180.;     // convert from degrees to radians
    associatedTypes.insert("pspmt");

}*/

    //Plot Dynode QDCs
    /*for (auto it = lowDynode.begin(); it != lowDynode.end(); it++) {
        if (DetectorDriver::get()->GetSysRootOutput()) {
            PSstruct.energy = (*it)->GetCalibratedEnergy();
            PSstruct.time = (*it)->GetTimeSansCfd();
            PSstruct.subtype = (*it)->GetChanID().GetSubtype();
            PSstruct.tag = (*it)->GetChanID().GetGroup();
            pixie_tree_event_->pspmt_vec_.emplace_back(PSstruct);
            PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
        }
        plot(DD_DYNODE_QDC, (*it)->GetTrace().GetQdc(), 0);
    }
    for (auto it = hiDynode.begin(); it != hiDynode.end(); it++) {
        if (DetectorDriver::get()->GetSysRootOutput()) {
            PSstruct.energy = (*it)->GetCalibratedEnergy();
            PSstruct.time = (*it)->GetTimeSansCfd();
            PSstruct.subtype = (*it)->GetChanID().GetSubtype();
            PSstruct.tag = (*it)->GetChanID().GetGroup();
            pixie_tree_event_->pspmt_vec_.emplace_back(PSstruct);
            PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
        }
        plot(DD_DYNODE_QDC, (*it)->GetTrace().GetQdc(), 1);
    }*/

    /*//set up position calculation for low / high gain yso signals and ion scint
    position_low.first = 0, position_low.second = 0;
    position_high.first = 0, position_high.second = 0;
//initalized all the things
    double energy = 0;
    double xa_l = 0, ya_l = 0, xb_l = 0, yb_l = 0;
    double xa_h = 0, ya_h = 0, xb_h = 0, yb_h = 0;
    double top_l = 0, top_r = 0, bottom_l = 0, bottom_r = 0;
    bool hasPosition_low = false, hasPosition_high = false, hasPosition_ion = false, hasUpstream = false,
    hasDeSi = false, hasVeto = false;

    plot(DD_MULTI, lowDynode.size(), 0);
    plot(DD_MULTI, hiDynode.size(), 1);

    plot(DD_MULTI, lowAnode.size(), 2);
    plot(DD_MULTI, hiAnode.size(), 3);


    double lowAnodeSum = 0;
    for (auto it = lowAnode.begin(); it != lowAnode.end(); it++) {
        //check signals energy vs threshold
        energy = (*it)->GetCalibratedEnergy();

        if (DetectorDriver::get()->GetSysRootOutput()) {
            PSstruct.energy = (*it)->GetCalibratedEnergy();
            PSstruct.time = (*it)->GetTimeSansCfd();
            PSstruct.subtype = (*it)->GetChanID().GetSubtype();
            PSstruct.tag = (*it)->GetChanID().GetGroup();
            pixie_tree_event_->pspmt_vec_.emplace_back(PSstruct);
            PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
        }

        if (energy < threshold_)
            continue;
        //parcel out position signals by tag
        if ((*it)->GetChanID().GetGroup() == "xa" && xa_l == 0) {
            xa_l = energy;
            lowAnodeSum += energy;
        }
        if ((*it)->GetChanID().GetGroup() == "xb" && xb_l == 0) {
            xb_l = energy;
            lowAnodeSum += energy;
        }
        if ((*it)->GetChanID().GetGroup() == "ya" && ya_l == 0) {
            ya_l = energy;
            lowAnodeSum += energy;
        }
        if ((*it)->GetChanID().GetGroup() == "yb" && yb_l == 0) {
            yb_l = energy;
            lowAnodeSum += energy;
        }
    }


    double highAnodeSum = 0;
    for (auto it = hiAnode.begin(); it != hiAnode.end(); it++) {

        if (DetectorDriver::get()->GetSysRootOutput()) {
            PSstruct.energy = (*it)->GetCalibratedEnergy();
            PSstruct.time = (*it)->GetTimeSansCfd();
            PSstruct.subtype = (*it)->GetChanID().GetSubtype();
            PSstruct.tag = (*it)->GetChanID().GetGroup();
            pixie_tree_event_->pspmt_vec_.emplace_back(PSstruct);
            PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
        }

        //check signals energy vs threshold
        energy = (*it)->GetCalibratedEnergy();
        if (energy < threshold_ || energy > 63000)
            continue;
        //parcel out position signals by tag
        if ((*it)->GetChanID().GetGroup() == "xa" && xa_h == 0){
            xa_h = energy;
            highAnodeSum += energy;
        }
        if ((*it)->GetChanID().GetGroup() == "xb" && xb_h == 0){
            xb_h = energy;
            highAnodeSum += energy;
        }
        if ((*it)->GetChanID().GetGroup() == "ya" && ya_h == 0){
            ya_h = energy;
            highAnodeSum += energy;
        }
        if ((*it)->GetChanID().GetGroup() == "yb" && yb_h == 0){
            yb_h = energy;
            highAnodeSum += energy;
        }
    }
    //compute position only if all 4 signals are present
    if (xa_l > 0 && xb_l > 0 && ya_l > 0 && yb_l > 0){
        hasPosition_low = true;
        position_low.first = CalculatePosition(xa_l, xb_l, ya_l, yb_l, vdtype_, rotation_).first;
        position_low.second  = CalculatePosition(xa_l, xb_l, ya_l, yb_l, vdtype_, rotation_).second;
        plot(DD_POS_LOW, position_low.first * positionScale_ + positionOffset_,
             position_low.second * positionScale_ + positionOffset_);
    }

    if (xa_h > 0 && xb_h > 0 && ya_h > 0 && yb_h > 0){
        hasPosition_high = true;
        position_high.first = CalculatePosition(xa_h, xb_h, ya_h, yb_h, vdtype_, rotation_).first;
        position_high.second = CalculatePosition(xa_h, xb_h, ya_h, yb_h, vdtype_, rotation_).second;
        plot(DD_POS_HIGH, position_high.first * positionScale_ + positionOffset_,
             position_high.second * positionScale_ + positionOffset_ );
    }

    //---------------VETO LOOP------------------------------------------------
    int numOfVetoChans = (int) (DetectorLibrary::get()->GetLocations("pspmt", "veto")).size();

    for(auto it = veto.begin(); it != veto.end(); it++){
        int loc = (*it)->GetChanID().GetLocation();
        plot(DD_PLASTIC_EN, (*it)->GetCalibratedEnergy(), loc);
        if ((*it)->GetCalibratedEnergy() > 1 && (*it)->GetCalibratedEnergy() < 10000){
            hasVeto = true;
        }

        if (DetectorDriver::get()->GetSysRootOutput()) {
            PSstruct.energy = (*it)->GetCalibratedEnergy();
            PSstruct.time = (*it)->GetTimeSansCfd();
            PSstruct.subtype = (*it)->GetChanID().GetSubtype();
            PSstruct.tag = (*it)->GetChanID().GetGroup();
            pixie_tree_event_->pspmt_vec_.emplace_back(PSstruct);
            PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
        }
    }

    //------------Positions from ion scintillator---------------------------------
        //using top - bottom and left - right computation scheme

    for(auto it = ionTrig.begin(); it != ionTrig.end(); it++){
        //check signals energy vs threshold

        if (DetectorDriver::get()->GetSysRootOutput()) {
            PSstruct.energy = (*it)->GetCalibratedEnergy();
            PSstruct.time = (*it)->GetTimeSansCfd();
            PSstruct.subtype = (*it)->GetChanID().GetSubtype();
            PSstruct.tag = (*it)->GetChanID().GetGroup();
            pixie_tree_event_->pspmt_vec_.emplace_back(PSstruct);
            PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
        }
        energy = (*it)->GetCalibratedEnergy();
        if (energy < 10)
            continue;

        // damm plotting of energies
        int loc = (*it)->GetChanID().GetLocation();
        plot(DD_PLASTIC_EN, (*it)->GetCalibratedEnergy(), loc + numOfVetoChans + 1); //max veto chan +1 for readablility

        //parcel out position signals by tag
        if ((*it)->GetChanID().GetGroup() == "black" && top_l == 0 )
            top_l = energy;
        if ((*it)->GetChanID().GetGroup() == "blue" && top_r == 0 )
            top_r = energy;
        if ((*it)->GetChanID().GetGroup() == "white" && bottom_l == 0 )
            bottom_l = energy;
        if ((*it)->GetChanID().GetGroup() == "green" && bottom_r == 0 )
            bottom_r = energy;


    }

    if (top_l > 0 && top_r > 0 && bottom_l > 0 && bottom_r > 0){
        hasPosition_ion = true;
        position_ion.first = (top_l + bottom_l - top_r - bottom_r) / (top_l + top_r + bottom_l + bottom_r);
        position_ion.second  = (top_l + top_r - bottom_l - bottom_r) / (top_l + top_r + bottom_l + bottom_r);
        plot(DD_POS_ION, position_ion.first * front_positionScale_ + front_positionOffset_,
             position_ion.second * front_positionScale_ + front_positionOffset_);
    }

    //----------------------------------------------------------------------------
    //------------Check Transmission efficiencies---------------------------------

        //check for valid upstream events, dE silicon events, and vetos for gating
    for(auto it = separatorScint.begin(); it != separatorScint.end(); it++){
        if (DetectorDriver::get()->GetSysRootOutput()) {
            PSstruct.energy = (*it)->GetCalibratedEnergy();
            PSstruct.time = (*it)->GetTimeSansCfd();
            PSstruct.subtype = (*it)->GetChanID().GetSubtype();
            PSstruct.tag = (*it)->GetChanID().GetGroup();
            pixie_tree_event_->pspmt_vec_.emplace_back(PSstruct);
            PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
        }
        if ((*it)->GetCalibratedEnergy() > 1 && (*it)->GetCalibratedEnergy() < 10000){
            hasUpstream = true;
        }
    }

    for (auto it = desi.begin(); it != desi.end(); it++){
        if (DetectorDriver::get()->GetSysRootOutput()) {
            PSstruct.energy = (*it)->GetCalibratedEnergy();
            PSstruct.time = (*it)->GetTimeSansCfd();
            PSstruct.subtype = (*it)->GetChanID().GetSubtype();
            PSstruct.tag = (*it)->GetChanID().GetGroup();
            pixie_tree_event_->pspmt_vec_.emplace_back(PSstruct);
            PSstruct = processor_struct::PSPMT_DEFAULT_STRUCT;
        }
        if((*it)->GetCalibratedEnergy() > 1 && (*it)->GetCalibratedEnergy() < 10000){
            hasDeSi = true;
        }
    }*/


    //Begin plotting
      /*   plot(DD_DYNODE_QDC, (*it)->GetTrace().GetQdc(), 0);
           plot(DD_DYNODE_QDC, (*it)->GetTrace().GetQdc(), 1);

        //plot valid YSO positions and dE silicon events gated on upstream events
        //plot upstream events gated on dE silicon
        //plot transmission efficiency from upstream to YSO and veto

    if(hasUpstream) {
        for(auto de_it = desi.begin(); de_it != desi.end(); de_it++) {
            plot(D_DESI_ENERGY,(*de_it)->GetCalibratedEnergy());
        }
        if(hasPosition_low){
            plot(DD_SEPAR_GATED_LOW, position_low.first * positionScale_ + positionOffset_,
                    position_low.second * positionScale_ + positionOffset_);
        }
        if(hasPosition_ion){
            plot(DD_SEPAR_GATED_ION, position_ion.first * positionScale_ + positionOffset_,
                    position_ion.second * positionScale_ + positionOffset_);
        }
    }

    if(hasDeSi) {
        for (auto it_sep = separatorScint.begin(); it_sep != separatorScint.end(); it_sep++) {
            if ((*it_sep)->GetChanID().GetGroup() =="left") {
                plot(DD_SEPAR_ENERGY, (*it_sep)->GetCalibratedEnergy(), 0);
            } else if ((*it_sep)->GetChanID().GetGroup() == "right") {
                plot(DD_SEPAR_ENERGY, (*it_sep)->GetCalibratedEnergy(), 1);
            }
        }

        if (hasPosition_low) {
            plot(DD_DESI_GATED_LOW, position_low.first * positionScale_ + positionOffset_,
                 position_low.second * positionScale_ + positionOffset_);
        }
        if (hasPosition_ion) {
            plot(DD_DESI_GATED_ION, position_ion.first * positionScale_ + positionOffset_,
                 position_ion.second * positionScale_ + positionOffset_);
        }
    }

    if(hasPosition_low) {
        for (auto de_it = desi.begin(); de_it != desi.end(); de_it++) {
            plot(D_DESI_YSO_GATED, (*de_it)->GetCalibratedEnergy());
        }
        for (auto it_sep = separatorScint.begin(); it_sep != separatorScint.end(); it_sep++) {
            if ((*it_sep)->GetChanID().GetGroup() =="left") {
                plot(DD_SEPAR_YSO_GATED, (*it_sep)->GetCalibratedEnergy(), 0);
            } else if ((*it_sep)->GetChanID().GetGroup() =="right") {
                plot(DD_SEPAR_YSO_GATED, (*it_sep)->GetCalibratedEnergy(), 1);
            }
        }
    }

    if(hasUpstream)
        plot(D_TRANS_EFF_YSO, 0);
    if(hasUpstream && hasPosition_ion)
        plot(D_TRANS_EFF_YSO, 1);
    if(hasUpstream && hasPosition_low)
        plot(D_TRANS_EFF_YSO, 2);
    if (hasUpstream && hasVeto)
        plot(D_TRANS_EFF_YSO, 3);



    if (!lowDynode.empty())
        plot(DD_DY_SUM_LG,lowDynode.front()->GetCalibratedEnergy(),lowAnodeSum);

    if (!hiDynode.empty())
        plot(DD_DY_SUM_HG,hiDynode.front()->GetCalibratedEnergy(),highAnodeSum);

    EndProcess();
    return true;

}
*/
