// -*- C++ -*-
//
// Package:    BRIL_ITsim/BIBAnalyzerRead
// Class:      BIBAnalyzerRead
//
/**\class BIBAnalyzerRead BIBAnalyzerRead.cc BRIL_ITsim/BIBAnalyzer/plugins/BIBAnalyzerRead.cc

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
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
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

//write histo data to hd5 file
#include "H5Cpp.h"
//#include "/afs/cern.ch/work/p/pkicsiny/private/cmssw/CMSSW_11_2_0_pre6/src/BRIL_ITsim/BIBAnalyzer/hd5/toy_hd5_extendible.cpp"
using namespace H5;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class BIBAnalyzerRead : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit BIBAnalyzerRead(const edm::ParameterSet&);
    ~BIBAnalyzerRead();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    uint32_t getModuleID(bool, unsigned int, unsigned int, unsigned int);
    // ----------member data ---------------------------
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> m_tokenClusters;
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis;
    edm::EDGetTokenT<std::vector<PSimHit>> m_tokenSimHits_low;
    edm::EDGetTokenT<std::vector<PSimHit>> m_tokenSimHits_high;
    edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> m_tokenSimLinks;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit>> m_tokenRecHits;
    //TrackerHitAssociator::Config trackerHitAssociatorConfig_;

   
    // the pointers to geometry, topology and clusters
    // these are members so all functions can access them without passing as argument
    const TrackerTopology* tTopo = NULL;
    const TrackerGeometry* tkGeom = NULL;
    const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19
    const edmNew::DetSetVector<SiPixelCluster>* clusters = NULL;
    const std::vector<PSimHit>* simhits_low = NULL;
    const std::vector<PSimHit>* simhits_high = NULL;
    const edmNew::DetSetVector<SiPixelRecHit>* rechits = NULL;

    //max bins of Counting histogram
    
    uint32_t m_nEBins;
    uint32_t m_maxEBin;

    /* histograms total (1) and per disk (8)
    * photon 22
    * electron/positron +/-11
    * muon +/-13
    * neutron/antineutron +/-2112
    * proton/antiproton +/-2212
    * kaon +/-321
    * pion +/-211
    */
   
    //event counter
    uint32_t m_nevents;
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
BIBAnalyzerRead::BIBAnalyzerRead(const edm::ParameterSet& iConfig)
	: m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters")))
        , m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis")))
        , m_tokenSimHits_low(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_low")))
        //, m_tokenSimHits_high(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_high")))
        , m_tokenRecHits(consumes<edmNew::DetSetVector<SiPixelRecHit>>(iConfig.getParameter<edm::InputTag>("rechits")))	
        , m_nEBins(iConfig.getUntrackedParameter<uint32_t>("nEBins")) 
        , m_maxEBin(iConfig.getUntrackedParameter<uint32_t>("maxEBin"))
       {

    std::cout << "BIBAnalyzerRead::BIBAnalyzerRead" << std::endl; 
    //now do what ever initialization is needed
    m_nevents = 0;

    std::cout << "constructor finished" << std::endl;
} // end of constructor


BIBAnalyzerRead::~BIBAnalyzerRead() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//


// ------------ method called once each job just before starting event loop  ------------
void BIBAnalyzerRead::beginJob() {

    std::cout << "BIBAnalyzerRead::beginJob()" << std::endl;
    std::cout << "beginjob finished" << std::endl;
}


// ------------ method called for each event  ------------
void BIBAnalyzerRead::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    std::cout << "BIBAnalyzerRead::analyze" << std::endl;
    //get the clusters
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> tclusters;
    iEvent.getByToken(m_tokenClusters, tclusters);
    
    //get the rechits
    edm::Handle<edmNew::DetSetVector<SiPixelRecHit>> trechits;
    iEvent.getByToken(m_tokenRecHits, trechits);
    
    //get the digis - COB 26.02.19
    edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
    iEvent.getByToken(m_tokenDigis, tdigis);
    
    //get the simhits - COB 30.11.20
    edm::Handle<std::vector<PSimHit>> tsimhits_low;
    iEvent.getByToken(m_tokenSimHits_low, tsimhits_low);
  //  edm::Handle<std::vector<PSimHit>> tsimhits_high;
   // iEvent.getByToken(m_tokenSimHits_high, tsimhits_high);   
    
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
    std::cout << "# digis: " << digis->size() << std::endl;

    simhits_low = tsimhits_low.product();
    std::cout << "# psimhits: " << simhits_low->size() << std::endl;

   // simhits_high = tsimhits_high.product();

    clusters = tclusters.product();
    std::cout << "# clusters: " << clusters->size() << std::endl;

    rechits = trechits.product(); 
    std::cout << "# rechits: " << rechits->size() << std::endl;

    //std::cout << "begin loops" << std::endl;

    /*
    * loop over the modules in the cluster collection (both TEPX and TFPX)
    */
/*
    for (typename edmNew::DetSetVector<SiPixelCluster>::const_iterator DSVit = clusters->begin(); DSVit != clusters->end(); DSVit++)    {
	std::cout << "event " << m_nevents << ", clusters" << std::endl;	
//	break;
    }
*/

    /*
    * loop over the modules in the rechit collection
    */
/*
    for (typename edmNew::DetSetVector<SiPixelRecHit>::const_iterator DSVit = rechits->begin(); DSVit != rechits->end(); DSVit++)    {
	std::cout << "event " << m_nevents << ", rechits" << std::endl;	
//	break;
    }
*/
    /*
    * loop over pixel digis
    * get charge and occupancy
    */
/*
    for (typename edm::DetSetVector<PixelDigi>::const_iterator DSVit = digis->begin(); DSVit != digis->end(); DSVit++) {
	std::cout << "event " << m_nevents << ", digis" << std::endl;
//	break;
    } // end of pixel digi loop
*/
  
    /*
    * loop over PSimHits low
    */
/*
    for (typename std::vector<PSimHit>::const_iterator DSVit = simhits_low->begin(); DSVit != simhits_low->end(); DSVit++) {
	std::cout << "event " << m_nevents << ", simhit low" << std::endl;
	break;
    }
*/
    /*
    * loop over PSimHits high
    */
/*
    for (typename std::vector<PSimHit>::const_iterator DSVit = simhits_high->begin(); DSVit != simhits_high->end(); DSVit++) {
	std::cout << "event " << m_nevents << ", simhit hi" << std::endl;
	break;
    }
*/
    m_nevents++;
}


// ------------ method called once each job just after ending the event loop  ------------
void BIBAnalyzerRead::endJob() {
    std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BIBAnalyzerRead::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

DEFINE_FWK_MODULE(BIBAnalyzerRead);
