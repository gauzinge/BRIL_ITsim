// -*- C++ -*-
//
// Package:    BRIL_ITsim/BIBAnalyzerOccupancy
// Class:      BIBAnalyzerOccupancy
//
/**\class BIBAnalyzerOccupancy BIBAnalyzerOccupancy.cc BRIL_ITsim/BIBAnalyzerOccupancy/plugins/BIBAnalyzerOccupancy.cc

Description: 
time of flight studies of BIB data at TEPX
plots produced: 
phi distribution 1D histo
X-Y distribution 2D histo
E spectrum 1D histo
tof-Q plots 2D histo

Implementation:
time of flight info is extracted from PSimHits that are stored in two branches (low + high)
*/
//
// Original Author:  Peter Kicsiny
//         Created:  Wed, 05 Dec 2020 15:10:16 GMT
//
//

// system include files
#include <algorithm>
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimTracker/SiPhase2Digitizer/plugins/Phase2TrackerDigitizerAlgorithm.h"

//matching recohits and psimhits
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include <TH2F.h>
#include <TTree.h>
#include <TStyle.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class BIBAnalyzerOccupancy : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit BIBAnalyzerOccupancy(const edm::ParameterSet&);
    ~BIBAnalyzerOccupancy();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    uint32_t getModuleID(bool, unsigned int, unsigned int, unsigned int);
    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis;
   
    // the pointers to geometry, topology and clusters
    // these are members so all functions can access them without passing as argument
    const TrackerTopology* tTopo = NULL;
    const TrackerGeometry* tkGeom = NULL;
    const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19

    const double occupancy_rd53b = 1.722332451499118e-04;  //100/(672*864) [%] 2x2 or 1x4 module size from RD53B manual
    const double occupancy = 1.706924480540207e-4;  //100/(1353*433) [%] size of actual TEPX module in CMSSW geo
    
    //1D histo for occupancy per disk per ring
    TH1F* m_histoOccAll[8][5];
    TH2F* m_histoRowCol[8][5];

    //event counter
    uint32_t m_nevents;

    //max channel counter
    int maxrow;
    int maxcol;

    //dict for storing occupied channels
    //keys are strings: "diskid_ringid_moduleid"
    //values are extendible vectors containing channel ids
    std::map<std::string, std::vector< int >> occupancy_dict;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BIBAnalyzerOccupancy::BIBAnalyzerOccupancy(const edm::ParameterSet& iConfig)
        : m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis")))
       {

    std::cout << "BIBAnalyzerOccupancy::BIBAnalyzerOccupancy" << std::endl; 
 
    //now do what ever initialization is needed
    m_nevents = 0;
    maxrow = 0;
    maxcol = 0;

    //std::cout << "constructor finished" << std::endl;
} // end of constructor


BIBAnalyzerOccupancy::~BIBAnalyzerOccupancy() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}

//
// member functions
//


// ------------ method called once each job just before starting event loop  ------------
void BIBAnalyzerOccupancy::beginJob() {

    edm::Service<TFileService> fs;
    
    //std::cout << "BIBAnalyzerOccupancy::beginJob()" << std::endl;
    /*
    * TEPX histograms for hits
    */    
    fs->file().cd("/");
    TFileDirectory oc = fs->mkdir("TEPX_occupancy");
    oc = fs->mkdir("TEPX_occupancy/Hits_occupancy");

    //loop over disks
    for (unsigned int i = 0; i < 8; i++) {

        // disk ID if [-4, 1] for -Z side and [1, 4] for +Z side
        int disk = (i < 4) ? i - 4 : i - 3;
        
	//TH2F 2D array
        for (unsigned int r = 1; r < 6; r++) {

               //occupancy plots     
	       m_histoOccAll[i][r-1] = oc.make<TH1F>(std::string("Occupancy for all particles on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for all particles", 48, .5, 48.5); 
     	    
	       // dummy histo for row-column plotting
     	       m_histoRowCol[i][r-1] = oc.make<TH2F>(std::string("TEPX module row-column indices on disk") + disk + " and Ring " + r, "TEPX module row-column indices", 3000, 0, 3000, 3000, 0, 3000);
        }
    }

    // dummy histo for row-column plotting
    //m_histoRowCol = oc.make<TH2F>("TEPX module row-column indices", "TEPX module row-column indices", 3000, 0, 3000, 3000, 0, 3000);
    
    //std::cout << "beginjob finished" << std::endl;
}


// ------------ method called for each event  ------------
void BIBAnalyzerOccupancy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    //std::cout << "BIBAnalyzerOccupancy::analyze" << std::endl;

    //get the digis - COB 26.02.19
    edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
    iEvent.getByToken(m_tokenDigis, tdigis);

    // Get the geometry
    edm::ESHandle<TrackerGeometry> tgeomHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get("idealForDigi", tgeomHandle);

    // Get the topology
    edm::ESHandle<TrackerTopology> tTopoHandle;
    iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

    //get the pointers to geometry, topology and clusters
    tTopo = tTopoHandle.product();
 
    //const TrackerGeometry* tkGeom = &(*geomHandle);
    tkGeom = tgeomHandle.product();
    digis = tdigis.product();
    //std::cout << "number of digis for event " << m_nevents << ": " << digis->size() << std::endl;
    /*
    * loop over pixel digi det sets (=modules)
    * get charge and occupancy
    */

    for (typename edm::DetSetVector<PixelDigi>::const_iterator DSVit = digis->begin(); DSVit != digis->end(); DSVit++) {
        
	//get the detid
        unsigned int rawid(DSVit->detId());
        DetId detId(rawid);

        //figure out the module type using the detID
        TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
        if (mType != TrackerGeometry::ModuleType::Ph2PXF && detId.subdetId() != PixelSubdetector::PixelEndcap)
            continue;

        //find out which layer, side and ring
        unsigned int side = (tTopo->pxfSide(detId));  // values are 1 and 2 for -+Z
        unsigned int disk = (tTopo->pxfDisk(detId)); // values are 1 to 12 for disks TFPX1 to TFPX 8  and TEPX1 to TEPX 4
        unsigned int ring = (tTopo->pxfBlade(detId)); // values are 1 to 5
        unsigned int module = (tTopo->pxfModule(detId)); // values are 1 to 48	
   	//std::cout << "outer loop over detsets, detset(=module) side: " << side << ", disk: " << disk  << ", ring: " << ring << ", module: " << module << std::endl;
	//a TEPX module (9-12 starting from origin on both sides)
        if (disk > 8) {
	
	    //std::cout << "this is a TEPX module" << std::endl;
            //the index in my histogram map
            int disk_id = -1;
            unsigned int ring_id = ring - 1;
            unsigned int module_id = module - 1;
        
	    //this is a TEPX hit on side 1 (-Z)
            if (side == 1) {
                disk_id = 12 - disk; // goes from 0 to 3
            }
            
	    //this is a TEPX hit on side 2 (+Z)
            else if (side == 2) {
                disk_id = 4 + disk - 9; // goes from 4 to 7
            }

            //find the geometry of the module associated to the digi
            const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
            if (!geomDetUnit)
                continue;

            //loop over the digis (the pixels of the module with signal) in TEPX module
            for (edm::DetSet<PixelDigi>::const_iterator digit = DSVit->begin(); digit != DSVit->end(); digit++) {
        
        	// get pixel charge in adc counts
                int digitCharge = digit->adc();
		int digitRow = digit->row();
		int digitCol = digit->column(); 	
		int digitChannel = PixelDigi::pixelToChannel(digitRow, digitCol);

    		//std::cout << "inner loop over digis, digi x: " << digitRow << " digi y: " << digitCol << " digi channel: " << digitChannel << std::endl;
		
		if (digitRow > maxrow){
			maxrow = digitRow;
		}
                if (digitCol > maxcol){
                        maxcol = digitCol;
                }

		// add to occupancy only if channel is occupied (redundant condition as only occupied pixels are recorded in the detset)
		if (digitCharge > 0){
 		
	 	    // get dict key
		    std::ostringstream dict_key_;
		    dict_key_ << "d" << disk_id+1 << "r" << ring_id+1; //<< "m" << module_id+1;
		    std::string dict_key = dict_key_.str(); 
 	            
		    // check if key already exists
		    if (occupancy_dict.find(dict_key) == occupancy_dict.end()) {
               		occupancy_dict[dict_key] = {};
           	    }
		
		    // append channel id only if it is not yet occupied by another particle (redundant too as digis are not tied to particles but to pixels)
            	    if(std::find(occupancy_dict[dict_key].begin(), occupancy_dict[dict_key].end(), digitChannel) == occupancy_dict[dict_key].end()){
                	occupancy_dict[dict_key].push_back(digitChannel);
			m_histoRowCol[disk_id][ring_id]->Fill(digit->row(), digit->column());
			m_histoOccAll[disk_id][ring_id]->AddBinContent(module_id, occupancy);
		    }
		}

            } // end of pixeldigi loop      
    	} // end of single TEPX module if
    } // end of pixel cluster loop
  
    m_nevents++;
}


// ------------ method called once each job just after ending the event loop  ------------
void BIBAnalyzerOccupancy::endJob() {
    std::cout << "max row index: " << maxrow << ", max col index: " << maxcol << std::endl;
    std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
    std::cout << "dict of occupied channels:" << std::endl;
    for (auto const& pair: occupancy_dict) {
        std::cout << "{" << pair.first << ": " << pair.second.size() << "}" << std::endl;
    }
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BIBAnalyzerOccupancy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(BIBAnalyzerOccupancy);
