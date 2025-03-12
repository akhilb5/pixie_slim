/** \file MtasProcessor.h
 *
 * Header file for MTAS analysis
 */

#ifndef __MTAS_PROCESSOR_H_
#define __MTAS_PROCESSOR_H_

#include "EventProcessor.h"
#include <vector>

class DetectorSummary;
class RawEvent;
class ChanEvent;
/**
 * \brief Handles detectors of type Mtas
 */
class MtasProcessor : public EventProcessor 
{
    private:
        DetectorSummary *mtasSummary; ///< all detectors of type Mtas
        DetectorSummary *siliSummary; 
        DetectorSummary *geSummary;
        DetectorSummary *sipmSummary; 
        DetectorSummary *logiSummary;
        DetectorSummary *refmodSummary; //added by Goetz

	double refmodEnergy;
	double implantEnergy;
	int implantMult;

        static bool isTapeMoveOn;
        static bool isMeasureOn;
        static bool isBkgOn;
        static bool isLightPulserOn;
        static bool isIrradOn;
        static unsigned cycleNumber;
    
        static double measureOnTime;
        double firstTime;

	bool isTriggerOnSignal;
	bool isTapeMoveOnSignal;
	bool isTapeMoveOffSignal;
	bool isMeasureOnSignal;
	bool isMeasureOffSignal;	
 	bool isBkgOnSignal;
	bool isBkgOffSignal;
 	bool isLightPulserOnSignal;
	bool isLightPulserOffSignal;	
 	bool isIrradOnSignal;
	bool isIrradOffSignal;
        int logicSignalsValue;

	//booleans for rings
	bool isCenter;
	bool isInner;
	bool isMiddle;
	bool isOuter;
	bool isCenterOnly;
	bool isInnerOnly;
	bool isMiddleOnly;
	bool isOuterOnly;
	bool isAll;

    public:
        MtasProcessor(); // no virtual c'tors
        virtual void DeclarePlots(void) const;
        virtual bool Process(RawEvent &event);
        
	bool GetIsTapeMove(void) const {return isTapeMoveOn;}
 	bool GetIsMeasure(void) const {return isMeasureOn;}
 	bool GetIsBackground(void) const {return isBkgOn;}
	unsigned GetCycleNumber(void) const {return cycleNumber;} 		    
    	
    private:
        void FillMtasMap();
        void FillSiliMap();
        void FillGeMap();
        void FillSipmMap();
        void FillRefModMapAndEnergy();
        void FillLogicMap();

        double maxLocation; 
	int nrOfCentralPMTs; 

	void SetCycleState();
        /*
        void FillMtasEnergyVectors();*/
        void SetIfOnlyRingBool();

        struct MtasData //to trzeba przerobic // Needs rewritten?
        {
            MtasData(std::string type="");
	    MtasData(ChanEvent *chan);
            
            std::string detSubtype;
            double energy;
            double calEnergy;
            double time;
            double location;

            //std::vector<int> trace; not necessary for MTAS
        };

        std::vector<ChanEvent*> logiList;
	std::vector<ChanEvent*> mtasList;
	std::vector<ChanEvent*> siliList;
	std::vector<ChanEvent*> geList;
	std::vector<ChanEvent*> sipmList;
	std::vector<ChanEvent*> refmodList;

	std::map<std::string, struct MtasData>  mtasMap;
	std::map<std::string, struct MtasData>  siliMap;
	std::map<std::string, struct MtasData>  geMap;
	std::map<std::string, struct MtasData>  sipmMap;
	std::map<std::string, struct MtasData>  refmodMap;
        std::map<std::string, struct MtasData>  logiMap;

        bool isBetaSignal;
        double betaTime;
    	double maxSiliconSignal;
};

#endif // __MTAS_PROCESSOR_H_
