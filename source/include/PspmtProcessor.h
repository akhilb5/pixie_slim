///@file PspmtProcessor.cpp
///@brief Processes information from a Position Sensitive PMT.
///@author A. Keeler, S. Go, S. V. Paulauskas
///@date July 8, 2018
#ifndef __PSPMTPROCESSOR_HPP__
#define __PSPMTPROCESSOR_HPP__

#include <string>
#include <map>
#include <vector>

#include "EventProcessor.h"
#include "RawEvent.h"

class DetectorSummary;
class RawEvent;
class ChanEvent;
//#include "ProcessorRootStruc.hpp"

///Class to handle processing of position sensitive pmts
class PspmtProcessor : public EventProcessor {
public:
    ///Default Constructor */
    PspmtProcessor();

    ///@brief Constructor that sets the configuration parameters
    virtual bool Init(DetectorDriver &driver, const std::string &ConfigFile=defaultConfigFile);
    /*PspmtProcessor(const std::string &vd, const double &yso_scale,
                   const unsigned int &yso_offset, const double &yso_threshold,
                   const double &front_scale,
                   const unsigned int &front_offset, const double &front_threshold, const double &rotation);
    */
    ///Default Destructor
    ~PspmtProcessor() {};

    ///Declare the plots used in the analysis
    void DeclarePlots(void);

    ///Preprocess the PSPMT data
    ///@param [in] event : the event to preprocess
    ///@return true if successful */
    bool Process(RawEvent &event);

    ///Enumeration that describes the different voltage dividers that have
    /// been used. Each board may have a different method to calculate the
    /// position of the interaction.
    enum VDTYPES {
        corners, sides, UNKNOWN
    };

    ///@return The x,y position of the interaction as calculated from the QDC
    /// found from the waveform in the WaveformAnalyzer, the on-board
    /// Pixie-16 energy filter, or the TraceFilterAnalyzer applying a
    /// trapezoidal filter to the trace.
    std::pair<double, double> GetPosition(const std::string &type) {
        if (type == "low")
            return position_low;
        else if (type == "high")
            return position_high;
        else if (type == "ion_scint")
            return position_ion;
        else
            return std::pair<double, double>(0., 0.);
    }

    ///@return The PSPMT Processor's TNAMED header. The Order is VD type then the software-based anode threshold.
    std::pair<std::string,std::string> GetPSPMTHeader(){
        return (make_pair(VDtypeStr,ThreshStr));
    }

private:

    static const std::string defaultConfigFile;

    double rotation;

    std::pair<double, double> position_low;
    std::pair<double, double> position_high;
    std::pair<double, double> position_ion;

    std::vector<ChanEvent*> hiDynode;
    std::vector<ChanEvent*> lowDynode;
    std::vector<ChanEvent*> hiAnode;
    std::vector<ChanEvent*> lowAnode;

    std::vector<ChanEvent*> veto;
    std::vector<ChanEvent*> ionTrig;
    std::vector<ChanEvent*> desi;
    std::vector<ChanEvent*> separatorScint;

    ///@brief A method to calculate the x position of the interaction with
    /// the scintillator
    ///@param[in] map : The map containing the set of x and y points that we
    /// will use to calculate the position.
    ///@return The x,y position of the interaction
    std::pair<double, double> CalculatePosition(double &xa, double &xb, double &ya,
                                                double &yb, const VDTYPES &vdtype, double &rot);


    VDTYPES vdtype_; ///< Local variable to store the type of voltage divider
    ///< we're using.
    double positionScale_; ///< The scale that we need for the DAMM output
    double front_positionScale_; ///< The scale that we need for the DAMM output
    unsigned int positionOffset_; ///< The offset that we need for the the
    unsigned int front_positionOffset_; ///< The offset that we need for the the
    ///< DAMM output
    double threshold_; ///< The threshold that the energy calculated by
    double front_threshold_; ///< The threshold that the energy calculated by
    ///< the Pixie-16 trapezoidal filter needs to reach
    ///< before we can analyze the signals.
    double rotation_; ///< rotation angle for Pspmt positions

    //processor_struct::PSPMT PSstruct; //!< PSPMT root Struct

    std::string VDtypeStr; //!< VD Type as a string
    std::string ThreshStr; //!< Threshold as a string

    struct PSData 
    {
        PSData(std::string type="");
	PSData(ChanEvent *chan);
        
        std::string subtype;
        double energy;
        double calEnergy;
        double time;
        double tag;

        //std::vector<int> trace;
    };

    void FillPspmtMap();
    std::map<std::string, struct PSData>  pspmtMap;
};

#endif // __PSPMTPROCESSOR_HPP__
