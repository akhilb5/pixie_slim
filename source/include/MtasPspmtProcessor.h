///@file MtasPspmtProcessor.cpp
///@brief Processes information from a Position Sensitive PMT.
///@author A. Keeler, S. Go, S. V. Paulauskas
///@date July 8, 2018
#ifndef __MtasPspmtProcessor_HPP__
#define __MtasPspmtProcessor_HPP__

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
class MtasPspmtProcessor : public EventProcessor {
public:
    ///Default Constructor */
    MtasPspmtProcessor();

    ///@brief Constructor that sets the configuration parameters
    //virtual bool Init(DetectorDriver &driver, const std::string &ConfigFile=defaultConfigFile) {return true;};
    //virtual bool Init();
    virtual bool Init(DetectorDriver &driver);
    /*MtasPspmtProcessor(const std::string &vd, const double &yso_scale,
                   const unsigned int &yso_offset, const double &yso_threshold,
                   const double &front_scale,
                   const unsigned int &front_offset, const double &front_threshold, const double &rotation);
    */
    ///Default Destructor
    ~MtasPspmtProcessor() {};

    ///Declare the plots used in the analysis
    void DeclarePlots(void) const;

    ///Preprocess the PSPMT data
    ///@param [in] event : the event to preprocess
    ///@return true if successful */
    bool Process(RawEvent &event);

    ///Enumeration that describes the different voltage dividers that have
    /// been used. Each board may have a different method to calculate the
    /// position of the interaction.
    enum DTTYPES {
        implant, diagnostic, UNKNOWN
    };
    ///@return The x,y position of the interaction as calculated from the QDC
    /// found from the waveform in the WaveformAnalyzer, the on-board
    /// Pixie-16 energy filter, or the TraceFilterAnalyzer applying a
    /// trapezoidal filter to the trace.
    /*std::pair<double, double> GetPosition(const std::string &type) {
        if (type == "low")
            return position_low;
        else if (type == "high")
            return position_high;
        else if (type == "ion_scint")
            return position_ion;
        else
            return std::pair<double, double>(0., 0.);
    }*/

    ///@return The PSPMT Processor's TNAMED header. The Order is VD type then the software-based anode threshold.
    std::pair<std::string,std::string> GetPSPMTHeader(){
        return (make_pair(VDtypeStr,ThreshStr));
    }

private:

    static const std::string defaultConfigFile;

    double rotation;

    std::pair<double, double> position_mtas;

    std::vector<ChanEvent*> pspmtList;

    ///@brief A method to calculate the x position of the interaction with
    /// the scintillator
    ///@param[in] map : The map containing the set of x and y points that we
    /// will use to calculate the position.
    ///@return The x,y position of the interaction
    std::pair<double, double> CalculatePosition(double &xa, double &xb, double &ya,
                                                double &yb, const DTTYPES &dttype);



    DTTYPES dttype_; ///< Local variable to store the type of voltage divider
    ///< we're using.
    double implantScale_; ///< The scale that we need for the DAMM output
    double diagnosticScale_; ///< The scale that we need for the DAMM output
    unsigned int implantOffset_; ///< The offset that we need for the the
    unsigned int diagnosticOffset_; ///< The offset that we need for the the
    ///< DAMM output
    double implantThreshold_; ///< The threshold that the energy calculated by
    double diagnosticThreshold_; ///< The threshold that the energy calculated by
    ///< the Pixie-16 trapezoidal filter needs to reach
    ///< before we can analyze the signals.
    
    
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
        double loc;
        std::string tag;
        //std::vector<int> trace; //needed or not?
    };

    void FillMtasPspmtMap();
    std::map<std::string, struct PSData>  mtasPspmtMap;
};

#endif // __MtasPspmtProcessor_HPP__
