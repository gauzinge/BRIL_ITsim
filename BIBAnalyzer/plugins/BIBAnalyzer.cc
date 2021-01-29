// -*- C++ -*-
//
// Package:    BRIL_ITsim/BIBAnalyzer
// Class:      BIBAnalyzer
//
/**\class BIBAnalyzer BIBAnalyzer.cc BRIL_ITsim/BIBAnalyzer/plugins/BIBAnalyzer.cc

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

#include <TH2F.h>
#include <TTree.h>
#include <TStyle.h>

//write histo data to hd5 file
#include <hdf5.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class BIBAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit BIBAnalyzer(const edm::ParameterSet&);
    ~BIBAnalyzer();

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
    edm::DetSetVector<PixelDigiSimLink>::const_iterator findSimLinkDetSet(unsigned int thedetid);
    std::map<std::tuple<float, float>, float>& getSimTrackId(edm::DetSetVector<PixelDigiSimLink>::const_iterator,
				        edm::DetSet<PixelDigi>::const_iterator digit,
					//edmNew::DetSet<SiPixelCluster>::const_iterator, 
					std::map<std::tuple<float, float>, float>& channel_charge_dict,
					unsigned int rawid,
					bool print);

    // the pointers to geometry, topology and clusters
    // these are members so all functions can access them without passing as argument
    const TrackerTopology* tTopo = NULL;
    const TrackerGeometry* tkGeom = NULL;
    const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19
    const edmNew::DetSetVector<SiPixelCluster>* clusters = NULL;
    const std::vector<PSimHit>* simhits_low = NULL;
    const std::vector<PSimHit>* simhits_high = NULL;
    const edm::DetSetVector<PixelDigiSimLink>* simlinks = NULL;

    //max bins of Counting histogram
    
    uint32_t m_nEBins;
    uint32_t m_maxEBin;
 
    //speed of light
    const double c = 1; //3e8;
    const double occupancy = 1.722332451499118e-04;//100/(672*864) [%] 2x2 or 1x4 module size from RD53B manual

    //2d histogram for occupancy
    TH1F* m_histoOccupancy[8][5];

    //array of TH2F for psimhits per disk per ring
    TH2F* m_histosPSimHits_low[8][5]; // 8 disks, 5 rings in each disc
    TH2F* m_histosPSimHits_high[8][5]; // 8 disks, 5 rings in each disc

    //histogram for time of flight vs Z coordinate
    TH2F* m_histoTofZ;
    
    //dict for particle tpye:number of hits
    std::map<int, int> ptypes_dict;

    //dict for simtrackid:charge
    std::map<std::tuple<float, float>, float> channel_charge_dict;

    /* histograms total (1) and per disk (8)
    * photon 22
    * electron/positron +/-11
    * muon +/-13
    * neutron/antineutron +/-2112
    * proton/antiproton +/-2212
    * kaon +/-321
    * pion +/-211
    */
    
    // 1D energy spectra per disc and total
    TH1F* m_spectrumAll[9];
    TH1F* m_spectrumPhot[9];
    TH1F* m_spectrumElPM[9];
    TH1F* m_spectrumMuPM[9];
    TH1F* m_spectrumNeut[9];
    TH1F* m_spectrumProt[9];
    TH1F* m_spectrumKaPM[9];
    TH1F* m_spectrumPiPM[9];
    TH1F* m_spectrumRest[9];

    //2D histos for global X-Y distributions per ptype per disk and total
    TH2F* m_histoXYAll[9];
    TH2F* m_histoXYPhot[9];
    TH2F* m_histoXYElPM[9];
    TH2F* m_histoXYMuPM[9];
    TH2F* m_histoXYNeut[9];
    TH2F* m_histoXYProt[9];
    TH2F* m_histoXYKaPM[9];
    TH2F* m_histoXYPiPM[9];
    TH2F* m_histoXYRest[9];

    //1D histos for Phi distributions per ptype per disk and total
    TH1F* m_histoPhiAll[9];
    TH1F* m_histoPhiPhot[9];
    TH1F* m_histoPhiElPM[9];
    TH1F* m_histoPhiMuPM[9];
    TH1F* m_histoPhiNeut[9];
    TH1F* m_histoPhiProt[9];
    TH1F* m_histoPhiKaPM[9];
    TH1F* m_histoPhiPiPM[9];
    TH1F* m_histoPhiRest[9];

    //2D histos for ToF - E distributions per ptype per disk per ring
    TH2F* m_histoTofAll[8][5];
    TH2F* m_histoTofPhot[8][5];
    TH2F* m_histoTofElPM[8][5];
    TH2F* m_histoTofMuPM[8][5];
    TH2F* m_histoTofNeut[8][5];
    TH2F* m_histoTofProt[8][5];
    TH2F* m_histoTofKaPM[8][5];
    TH2F* m_histoTofPiPM[8][5];
    TH2F* m_histoTofRest[8][5];
    
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
BIBAnalyzer::BIBAnalyzer(const edm::ParameterSet& iConfig)
	: m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters")))
        , m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis")))
        , m_tokenSimHits_low(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_low")))
        , m_tokenSimHits_high(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_high")))
        , m_tokenSimLinks(consumes<edm::DetSetVector<PixelDigiSimLink>>(iConfig.getParameter<edm::InputTag>("simlinks")))
	, m_nEBins(iConfig.getUntrackedParameter<uint32_t>("nEBins")) 
        , m_maxEBin(iConfig.getUntrackedParameter<uint32_t>("maxEBin")) {
    //now do what ever initialization is needed
    m_nevents = 0;
}


BIBAnalyzer::~BIBAnalyzer() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//


// ------------ method called once each job just before starting event loop  ------------
void BIBAnalyzer::beginJob() {

    edm::Service<TFileService> fs;
    
    /*
    * TEPX histograms for hits
    */    
    fs->file().cd("/");
    TFileDirectory td = fs->mkdir("TEPX_tof");
    td = fs->mkdir("TEPX_tof/Hits_tof");

    fs->file().cd("/");
    TFileDirectory pt = fs->mkdir("TEPX_ptype");
    pt = fs->mkdir("TEPX_ptype/Hits_ptype");
 
    fs->file().cd("/");
    TFileDirectory oc = fs->mkdir("TEPX_occupancy");
    oc = fs->mkdir("TEPX_occupancy/Hits_occupancy");

    fs->file().cd("/");
    TFileDirectory ta = fs->mkdir("TEPX_tof_arrays");
    ta = fs->mkdir("TEPX_tof_arrays/Hits_tof_arrays");

    const double phi_min = -3.454;
    const double phi_max = -1*phi_min;
    const unsigned int phi_nbins = 220;
    
    //loop over disks
    for (unsigned int i = 0; i < 8; i++) {

        // disk ID if [-4, 1] for -Z side and [1, 4] for +Z side
        int disk = (i < 4) ? i - 4 : i - 3;
        
        //E spectra
        m_spectrumAll[i] = pt.make<TH1F>(std::string("E spectrum of all particles on disk ") + disk, std::string("E spectrum of all particles on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumPhot[i] = pt.make<TH1F>(std::string("E spectrum of photons on disk ") + disk, std::string("E spectrum of photons on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumElPM[i] = pt.make<TH1F>(std::string("E spectrum of e+/- on disk ") + disk, std::string("E spectrum of e+/- on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumMuPM[i] = pt.make<TH1F>(std::string("E spectrum of muons+/- on disk ") + disk, std::string("E spectrum of muons+/- on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumNeut[i] = pt.make<TH1F>(std::string("E spectrum of neutrons on disk ") + disk, std::string("E spectrum of neutrons on disk ") + disk, m_nEBins, 0, m_maxEBin);
	m_spectrumProt[i] = pt.make<TH1F>(std::string("E spectrum of protons on disk ") + disk, std::string("E spectrum of protons on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumKaPM[i] = pt.make<TH1F>(std::string("E spectrum of kaons+/- on disk ") + disk, std::string("E spectrum of kaons+/- on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumPiPM[i] = pt.make<TH1F>(std::string("E spectrum of pions+/- on disk ") + disk, std::string("E spectrum of pions+/- on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumRest[i] = pt.make<TH1F>(std::string("E spectrum of residuals on disk ") + disk, std::string("E spectrum of residuals on disk ") + disk, m_nEBins, 0, m_maxEBin);

	//X-Y distributions
        m_histoXYAll[i] = pt.make<TH2F>(std::string("X vs Y global coordinates for all particles on disk ") + disk,
				        std::string("X vs Y global coordinates for all particles on disk ") + disk, 1000, -50.0, 50.0, 1000, -50.0, 50.0);
        m_histoXYPhot[i] = pt.make<TH2F>(std::string("X vs Y global coordinates for photons on disk ") + disk,
					 std::string("X vs Y global coordinates for photons on disk ") + disk, 1000, -50.0, 50.0, 1000, -50.0, 50.0);
        m_histoXYElPM[i] = pt.make<TH2F>(std::string("X vs Y global coordinates for e+/- on disk ") + disk,
					 std::string("X vs Y global coordinates for e+/- on disk ") + disk, 1000, -50.0, 50.0, 1000, -50.0, 50.0);
        m_histoXYMuPM[i] = pt.make<TH2F>(std::string("X vs Y global coordinates for muons+/- on disk ") + disk,
					 std::string("X vs Y global coordinates for muons+/- on disk ") + disk, 1000, -50.0, 50.0, 1000, -50.0, 50.0);
        m_histoXYNeut[i] = pt.make<TH2F>(std::string("X vs Y global coordinates for neutrons on disk ") + disk,
					 std::string("X vs Y global coordinates for neutrons on disk ") + disk, 1000, -50.0, 50.0, 1000, -50.0, 50.0);
	m_histoXYProt[i] = pt.make<TH2F>(std::string("X vs Y global coordinates for protons on disk ") + disk,
                                         std::string("X vs Y global coordinates for protons on disk ") + disk, 1000, -50.0, 50.0, 1000, -50.0, 50.0);
        m_histoXYKaPM[i] = pt.make<TH2F>(std::string("X vs Y global coordinates for kaons+/- on disk ") + disk,
                                         std::string("X vs Y global coordinates for kaons+/- on disk ") + disk, 1000, -50.0, 50.0, 1000, -50.0, 50.0);
        m_histoXYPiPM[i] = pt.make<TH2F>(std::string("X vs Y global coordinates for pions+/- on disk ") + disk,
                                         std::string("X vs Y global coordinates for pions+/- on disk ") + disk, 1000, -50.0, 50.0, 1000, -50.0, 50.0);
        m_histoXYRest[i] = pt.make<TH2F>(std::string("X vs Y global coordinates for residuals on disk ") + disk,
                                         std::string("X vs Y global coordinates for residuals on disk ") + disk, 1000, -50.0, 50.0, 1000, -50.0, 50.0);

        //Phi distributions
	m_histoPhiAll[i] = pt.make<TH1F>(std::string("Phi distribution of all particles on disk ") + disk, std::string("Phi distribution of all particles on disk ") + disk, phi_nbins, phi_min, phi_max);
        m_histoPhiPhot[i] = pt.make<TH1F>(std::string("Phi distribution of photons on disk ") + disk, std::string("Phi distribution of photons on disk ") + disk, phi_nbins, phi_min, phi_max);
        m_histoPhiElPM[i] = pt.make<TH1F>(std::string("Phi distribution of e+/- on disk ") + disk, std::string("Phi distribution of e+/- on disk ") + disk, phi_nbins, phi_min, phi_max);
        m_histoPhiMuPM[i] = pt.make<TH1F>(std::string("Phi distribution of muons+/- on disk ") + disk, std::string("Phi distribution of muons+/- on disk ") + disk, phi_nbins, phi_min, phi_max);
        m_histoPhiNeut[i] = pt.make<TH1F>(std::string("Phi distribution of neutrons on disk ") + disk, std::string("Phi distribution of neutrons on disk ") + disk, phi_nbins, phi_min, phi_max);
	m_histoPhiProt[i] = pt.make<TH1F>(std::string("Phi distribution of protons on disk ") + disk, std::string("Phi distribution of protons on disk ") + disk, phi_nbins, phi_min, phi_max);
        m_histoPhiKaPM[i] = pt.make<TH1F>(std::string("Phi distribution of kaons+/- on disk ") + disk, std::string("Phi distribution of kaons+/- on disk ") + disk, phi_nbins, phi_min, phi_max);
        m_histoPhiPiPM[i] = pt.make<TH1F>(std::string("Phi distribution of pions+/- on disk ") + disk, std::string("Phi distribution of pions+/- on disk ") + disk, phi_nbins, phi_min, phi_max);
        m_histoPhiRest[i] = pt.make<TH1F>(std::string("Phi distribution of residuals on disk ") + disk, std::string("Phi distribution of residuals on disk ") + disk, phi_nbins, phi_min, phi_max);


	//TH2F 2D array
        for (unsigned int r = 1; r < 6; r++) {

               //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
               m_histosPSimHits_low[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight (low) distribution for Disc ") + disk + " and Ring " + r + "",
							   std::string("Charge - time of flight (low) distribution for Disc ") + disk + " and Ring " + r,
						 	   550, -50, 500, 100, 0, 3e-3);
	       m_histosPSimHits_high[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight (high) distribution for Disc ") + disk + " and Ring " + r + "", 
                                                           std::string("Charge - time of flight (high) distribution for Disc ") + disk + " and Ring " + r, 
                                                           550, -50, 500, 100, 0, 3e-3);

	      //ToF - E plots
 	      m_histoTofAll[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for all particles on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for all particles on disk ") + disk + " and Ring " + r, 550, -50, 500, 100, 0, 3e-3);
              m_histoTofPhot[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for photons on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for photons on disk ") + disk + " and Ring " + r, 550, -50, 500, 100, 0, 3e-3);
              m_histoTofElPM[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for e+/- on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for e+/- on disk ") + disk + " and Ring " + r, 550, -50, 500, 100, 0, 3e-3);
              m_histoTofMuPM[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for muons+/- on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for muons+/- on disk ") + disk + " and Ring " + r, 550, -50, 500, 100, 0, 3e-3);
              m_histoTofNeut[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for neutrons on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for neutrons on disk ") + disk + " and Ring " + r, 550, -50, 500, 100, 0, 3e-3);
              m_histoTofProt[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for protons on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for protons on disk ") + disk + " and Ring " + r, 550, -50, 500, 100, 0, 3e-3);
              m_histoTofKaPM[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for kaons+/- on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for kaons+/- on disk ") + disk + " and Ring " + r, 550, -50, 500, 100, 0, 3e-3);
              m_histoTofPiPM[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for pions+/- on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for pions+/- on disk ") + disk + " and Ring " + r, 550, -50, 500, 100, 0, 3e-3);
	      m_histoTofRest[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for residuals on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for residuals on disk ") + disk + " and Ring " + r, 550, -50, 500, 100, 0, 3e-3);
             
	      //occupancy 
	      m_histoOccupancy[i][r-1] = oc.make<TH1F>(std::string("Occupancy on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis", 48, .5, 48.5);
        }
    }

    //total E and XY distributions per particle type and for everything
    m_spectrumAll[8] = pt.make<TH1F>("E spectrum of all particles on all disks", "E spectrum of all particles on all disks", m_nEBins, 0, m_maxEBin);
    m_spectrumPhot[8] = pt.make<TH1F>("E spectrum of photons on all disks", "E spectrum of photons on all disks", m_nEBins, 0, m_maxEBin);
    m_spectrumElPM[8] = pt.make<TH1F>("E spectrum of e+/- on all disks", "E spectrum of e+/- on all disks", m_nEBins, 0, m_maxEBin);
    m_spectrumMuPM[8] = pt.make<TH1F>("E spectrum of muons+/- on all disks", "E spectrum of muons+/- on all disks", m_nEBins, 0, m_maxEBin);
    m_spectrumNeut[8] = pt.make<TH1F>("E spectrum of neutrons on all disks", "E spectrum of neutrons on all disks", m_nEBins, 0, m_maxEBin);
    m_spectrumProt[8] = pt.make<TH1F>("E spectrum of protons on all disks", "E spectrum of protons on all disks", m_nEBins, 0, m_maxEBin);
    m_spectrumKaPM[8] = pt.make<TH1F>("E spectrum of kaons+/- on all disks", "E spectrum of kaons+/- on all disks", m_nEBins, 0, m_maxEBin);
    m_spectrumPiPM[8] = pt.make<TH1F>("E spectrum of pions+/- on all disks", "E spectrum of pions+/- on all disks", m_nEBins, 0, m_maxEBin);
    m_spectrumRest[8] = pt.make<TH1F>("E spectrum of residuals on all disks", "E spectrum of residuals on all disks", m_nEBins, 0, m_maxEBin);

    m_histoXYAll[8] = pt.make<TH2F>("X vs Y global coordinates for all particles on all disks",
				    "X vs Y global coordinates for all particles on all disks", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    m_histoXYPhot[8] = pt.make<TH2F>("X vs Y global coordinates for photons on all disks",
				     "X vs Y global coordinates for photons on all disks", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    m_histoXYElPM[8] = pt.make<TH2F>("X vs Y global coordinates for e+/- on all disks",
				     "X vs Y global coordinates for e+/- on all disks", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    m_histoXYMuPM[8] = pt.make<TH2F>("X vs Y global coordinates for muons+/- on all disks",
                                     "X vs Y global coordinates for muons+/- on all disks", 1000, -50.0, 50.0, 1000, -50.0, 50.0);    
    m_histoXYNeut[8] = pt.make<TH2F>("X vs Y global coordinates for neutrons on all disks",
				     "X vs Y global coordinates for neutrons on all disks", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    m_histoXYProt[8] = pt.make<TH2F>("X vs Y global coordinates for protons on all disks",
				     "X vs Y global coordinates for protons on all disks", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    m_histoXYKaPM[8] = pt.make<TH2F>("X vs Y global coordinates for kaons+/- on all disks",
				     "X vs Y global coordinates for kaons+/- on all disks", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    m_histoXYPiPM[8] = pt.make<TH2F>("X vs Y global coordinates for pions+/- on all disks",
				     "X vs Y global coordinates for pions+/- on all disks", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    m_histoXYRest[8] = pt.make<TH2F>("X vs Y global coordinates for residuals on all disks",
				     "X vs Y global coordinates for residuals on all disks", 1000, -50.0, 50.0, 1000, -50.0, 50.0);

    //phi distributions
    m_histoPhiAll[8] = pt.make<TH1F>("Phi distribution of all particles all disks", "Phi distribution of all particles on all disks", phi_nbins, phi_min, phi_max);
    m_histoPhiPhot[8] = pt.make<TH1F>("Phi distribution of photons on all disks", "Phi distribution of photons on all disks", phi_nbins, phi_min, phi_max);
    m_histoPhiElPM[8] = pt.make<TH1F>("Phi distribution of e+/- on all disks", "Phi distribution of e+/- on all disks", phi_nbins, phi_min, phi_max);
    m_histoPhiMuPM[8] = pt.make<TH1F>("Phi distribution of muons+/- on all disks", "Phi distribution of muons+/- on all disks", phi_nbins, phi_min, phi_max);
    m_histoPhiNeut[8] = pt.make<TH1F>("Phi distribution of neutrons on all disks", "Phi distribution of neutrons on all disks", phi_nbins, phi_min, phi_max);
    m_histoPhiProt[8] = pt.make<TH1F>("Phi distribution of protons on all disks", "Phi distribution of protons on all disks", phi_nbins, phi_min, phi_max);
    m_histoPhiKaPM[8] = pt.make<TH1F>("Phi distribution of kaons+/- on all disks", "Phi distribution of kaons+/- on all disks", phi_nbins, phi_min, phi_max);
    m_histoPhiPiPM[8] = pt.make<TH1F>("Phi distribution of pions+/- on all disks", "Phi distribution of pions+/- on all disks", phi_nbins, phi_min, phi_max);
    m_histoPhiRest[8] = pt.make<TH1F>("Phi distribution of residuals on all disks", "Phi distribution of residuals on all disks", phi_nbins, phi_min, phi_max);

    // global coord. vs tof
    m_histoTofZ = td.make<TH2F>("ToF vs Z", "Time of flight vs. global Z coordinate", 6000, -300, 300, 6500, -50, 600);
}


// ------------ method called for each event  ------------
void BIBAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    //get the clusters
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> tclusters;
    iEvent.getByToken(m_tokenClusters, tclusters);

    //get the digis - COB 26.02.19
    edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
    iEvent.getByToken(m_tokenDigis, tdigis);

    //get the simhits - COB 30.11.20
    edm::Handle<std::vector<PSimHit>> tsimhits_low;
    iEvent.getByToken(m_tokenSimHits_low, tsimhits_low);
    edm::Handle<std::vector<PSimHit>> tsimhits_high;
    iEvent.getByToken(m_tokenSimHits_high, tsimhits_high);   

    //get the simlinks
    edm::Handle<edm::DetSetVector<PixelDigiSimLink>> tsimlinks;
    iEvent.getByToken(m_tokenSimLinks, tsimlinks);

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
    simhits_low = tsimhits_low.product();
    simhits_high = tsimhits_high.product();
    simlinks = tsimlinks.product();
    clusters = tclusters.product();

    /*
    * loop over pixel clusters
    * get simtrackid
    * get charge
    * create a dict: {simtrackid: charge}
    *
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
	
	//a TEPX module (9-12 starting from origin on both sides)
        if (disk > 8) {

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

            //loop over the clusters in TEPX module
            for (edm::DetSet<PixelDigi>::const_iterator digit = DSVit->begin(); digit != DSVit->end(); digit++) {
               
        	// get pixel charge in adc counts
                int digitCharge = digit->adc();
		
		// add to occupancy
		if (digitCharge > 0){
			m_histoOccupancy[disk_id][ring_id]->AddBinContent(module_id, occupancy);
		}
            } // end of pixeldigi loop      
    	} // end of single TEPX module if
    } // end of pixel cluster loop

    /*
    * loop over PSimHits low
    */
    for (typename std::vector<PSimHit>::const_iterator DSVit = simhits_low->begin(); DSVit != simhits_low->end(); DSVit++) {

        //get the detid
        unsigned int rawid(DSVit->detUnitId());
        DetId detId(rawid);

        //module type => need phase 2 pixel forward module, in endcap
        TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
        if (mType != TrackerGeometry::ModuleType::Ph2PXF && detId.subdetId() != PixelSubdetector::PixelEndcap) 
            continue;

        //find out which layer, side and ring
        unsigned int side = (tTopo->pxfSide(detId));  // values are 1 and 2 for -+Z
        unsigned int disk = (tTopo->pxfDisk(detId)); // values are 1 to 12 for disks TFPX1 to TFPX 8  and TEPX1 to TEPX 4
        unsigned int ring = (tTopo->pxfBlade(detId)); // values are 1 to 5
        //unsigned int module = (tTopo->pxfModule(detId)); // values are 1 to 48	
	
	//a TEPX module (9-12 starting from origin on both sides)
        if (disk > 8) {

            //the index in my histogram map
            int disk_id = -1;
            unsigned int ring_id = ring - 1;
            //unsigned int module_id = module - 1;
        
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

	    // get local coordinates of PSimHit
	    Point3DBase<float, LocalTag> local_coords = DSVit->localPosition();
	    
	    // PixelGeomDetUnit is derived from GeomDet, transform local to global coords
	    GlobalPoint global_coords = geomDetUnit->toGlobal(local_coords);
      
            //get particle type
            int ptype = DSVit->particleType();
	    if (ptypes_dict.find(ptype) == ptypes_dict.end()) {
 	    	ptypes_dict[ptype] = 1;
	    }
	    else{
		ptypes_dict[ptype] += 1;
	    }

 	    // get energy deposit in electrons and convert to adc charge
	    float signal_in_eloss = DSVit->energyLoss(); // o(1e-3)
	    //float signal_threshold = 1000.0;
            //int signal_in_adc = Phase2TrackerDigitizerAlgorithm::convertSignalToAdc(detId, signal_in_eloss, signal_threshold);
 	    
	    //debug info
            //std::cout << "---- Hit info:" << std::endl;
	    //std::cout << "hit tof: " << DSVit->timeOfFlight() << std::endl;
	    //std::cout << "global position: " << global_coords.z() << std::endl;		
	    //std::cout << "disk id: " <<	disk_id << std::endl;
	    //std::cout << "particle type: " << DSVit->particleType() << std::endl;
	
	    // fill up particle type dependent histos
	    switch (std::abs(ptype)){
		case 22: // photon
                    m_spectrumPhot[8]->Fill(c*DSVit->pabs()); 
	            m_histoXYPhot[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiPhot[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x())); 
        	    m_spectrumPhot[disk_id]->Fill(c*DSVit->pabs());
            	    m_histoXYPhot[disk_id]->Fill(global_coords.x(), global_coords.y());
		    m_histoPhiPhot[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
		    m_histoTofPhot[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break;
                 case 11: // electron/positron   
                    m_spectrumElPM[8]->Fill(c*DSVit->pabs());
                    m_histoXYElPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiElPM[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumElPM[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYElPM[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiElPM[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofElPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break; 
		 case 13: // muon+/-
                    m_spectrumMuPM[8]->Fill(c*DSVit->pabs());
                    m_histoXYMuPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiMuPM[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumMuPM[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYMuPM[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiMuPM[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofMuPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break; 
		case 2112: // neutron/antineutron 
                    m_spectrumNeut[8]->Fill(c*DSVit->pabs());
                    m_histoXYNeut[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiNeut[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumNeut[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYNeut[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiNeut[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofNeut[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break; 
		case 2212: // proton/antiproton
                    m_spectrumProt[8]->Fill(c*DSVit->pabs());
                    m_histoXYProt[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiProt[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
		    m_spectrumProt[disk_id]->Fill(c*DSVit->pabs());
		    m_histoXYProt[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiProt[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
	    	    m_histoTofProt[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break;
                case 321: // kaon+/-
                    m_spectrumKaPM[8]->Fill(c*DSVit->pabs()); 
                    m_histoXYKaPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiKaPM[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumKaPM[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYKaPM[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiKaPM[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
	            m_histoTofKaPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
                    break;
                case 211: // pion+/-
                    m_spectrumPiPM[8]->Fill(c*DSVit->pabs()); 
                    m_histoXYPiPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiPiPM[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumPiPM[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYPiPM[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiPiPM[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofPiPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break; 
	        default: // none of the previous
                    m_spectrumRest[8]->Fill(c*DSVit->pabs());
                    m_histoXYRest[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiRest[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumRest[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYRest[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiRest[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofRest[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
                    break;	
	}
		
	    // fill up all particle histos	
	    m_spectrumAll[8]->Fill(c*DSVit->pabs());	    
            m_histoXYAll[8]->Fill(global_coords.x(), global_coords.y());
            m_histoPhiAll[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
	    m_spectrumAll[disk_id]->Fill(c*DSVit->pabs());
            m_histoXYAll[disk_id]->Fill(global_coords.x(), global_coords.y());
            m_histoPhiAll[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
            m_histoTofAll[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
	    m_histoTofZ->Fill(global_coords.z(), DSVit->timeOfFlight());

            // fill up tof histos	
            m_histosPSimHits_low[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;	

	}
    }

    /*
    * loop over PSimHits high
    */
    for (typename std::vector<PSimHit>::const_iterator DSVit = simhits_high->begin(); DSVit != simhits_high->end(); DSVit++) {

        //get the detid
        unsigned int rawid(DSVit->detUnitId());
        DetId detId(rawid);

        //module type => need phase 2 pixel forward module, in endcap
        TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
        if (mType != TrackerGeometry::ModuleType::Ph2PXF && detId.subdetId() != PixelSubdetector::PixelEndcap) 
            continue;

        //find out which layer, side and ring
        unsigned int side = (tTopo->pxfSide(detId));  // values are 1 and 2 for -+Z
        unsigned int disk = (tTopo->pxfDisk(detId)); // values are 1 to 12 for disks TFPX1 to TFPX 8  and TEPX1 to TEPX 4
        unsigned int ring = (tTopo->pxfBlade(detId)); // values are 1 to 5
        //unsigned int module = (tTopo->pxfModule(detId)); // values are 1 to 48	
	
	//a TEPX module (9-12 starting from origin on both sides)
        if (disk > 8) {

            //the index in my histogram map
            int disk_id = -1;
            unsigned int ring_id = ring - 1;
            //unsigned int module_id = module - 1;
        
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

	    // get local coordinates of PSimHit
	    Point3DBase<float, LocalTag> local_coords = DSVit->localPosition();
	   
	    // PixelGeomDetUnit is derived from GeomDet, transform local to global coords
	    GlobalPoint global_coords = geomDetUnit->toGlobal(local_coords);
      
            //get particle type
            int ptype = DSVit->particleType();
	    if (ptypes_dict.find(ptype) == ptypes_dict.end()) {
 	    	ptypes_dict[ptype] = 1;
	    }
	    else{
		ptypes_dict[ptype] += 1;
	    }

	    // get energy deposit in electrons and convert to adc charge
	    float signal_in_eloss = DSVit->energyLoss(); // o(1e-3)
	    //float signal_threshold = 1000.0;
            //int signal_in_adc = Phase2TrackerDigitizerAlgorithm::convertSignalToAdc(detId, signal_in_eloss, signal_threshold);
 	    
 	    //debug info
            //std::cout << "---- Hit info:" << std::endl;
            //if (std::abs(DSVit->timeOfFlight()-90) < 10 && std::abs(signal_in_eloss - 0.00005) < 0.000025){
		// std::cout << "particle type: " << DSVit->particleType() << " hit tof: " << DSVit->timeOfFlight() << std::endl;
	    //}
	    //std::cout << "hit tof: " << DSVit->timeOfFlight() << std::endl;
	    //std::cout << "global position: " << global_coords.z() << std::endl;		
	    //std::cout << "disk id: " <<	disk_id << std::endl;
	    //std::cout << "particle type: " << DSVit->particleType() << std::endl;
            //std::cout << "momentum: " << DSVit->pabs() << std::endl;
	
	    // fill up particle type dependent histos
	    switch (std::abs(ptype)){
		case 22: // photon
                    m_spectrumPhot[8]->Fill(c*DSVit->pabs()); 
	            m_histoXYPhot[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiPhot[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x())); 
        	    m_spectrumPhot[disk_id]->Fill(c*DSVit->pabs());
            	    m_histoXYPhot[disk_id]->Fill(global_coords.x(), global_coords.y());
		    m_histoPhiPhot[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
		    m_histoTofPhot[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break;
                 case 11: // electron/positron   
                    m_spectrumElPM[8]->Fill(c*DSVit->pabs());
                    m_histoXYElPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiElPM[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumElPM[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYElPM[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiElPM[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofElPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break; 
		 case 13: // muon+/-
                    m_spectrumMuPM[8]->Fill(c*DSVit->pabs());
                    m_histoXYMuPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiMuPM[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumMuPM[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYMuPM[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiMuPM[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofMuPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break; 
		case 2112: // neutron/antineutron 
                    m_spectrumNeut[8]->Fill(c*DSVit->pabs());
                    m_histoXYNeut[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiNeut[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumNeut[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYNeut[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiNeut[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofNeut[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break; 
		case 2212: // proton/antiproton
                    m_spectrumProt[8]->Fill(c*DSVit->pabs());
                    m_histoXYProt[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiProt[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
		    m_spectrumProt[disk_id]->Fill(c*DSVit->pabs());
		    m_histoXYProt[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiProt[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
	    	    m_histoTofProt[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break;
                case 321: // kaon+/-
                    m_spectrumKaPM[8]->Fill(c*DSVit->pabs()); 
                    m_histoXYKaPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiKaPM[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumKaPM[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYKaPM[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiKaPM[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
	            m_histoTofKaPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
                    break;
                case 211: // pion+/-
                    m_spectrumPiPM[8]->Fill(c*DSVit->pabs()); 
                    m_histoXYPiPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiPiPM[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumPiPM[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYPiPM[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiPiPM[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofPiPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
		    break;
                default: // none of the previous
                    m_spectrumRest[8]->Fill(c*DSVit->pabs());
                    m_histoXYRest[8]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiRest[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_spectrumRest[disk_id]->Fill(c*DSVit->pabs());
                    m_histoXYRest[disk_id]->Fill(global_coords.x(), global_coords.y());
                    m_histoPhiRest[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
                    m_histoTofRest[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
                    break;
		}
		
	    // fill up all particle histos	
	    m_spectrumAll[8]->Fill(c*DSVit->pabs());	    
            m_histoXYAll[8]->Fill(global_coords.x(), global_coords.y());
            m_histoPhiAll[8]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
	    m_spectrumAll[disk_id]->Fill(c*DSVit->pabs());
            m_histoXYAll[disk_id]->Fill(global_coords.x(), global_coords.y());
            m_histoPhiAll[disk_id]->Fill(TMath::ATan2(global_coords.y(), global_coords.x()));
            m_histoTofAll[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
	    m_histoTofZ->Fill(global_coords.z(), DSVit->timeOfFlight());

            // fill up tof histos	
            m_histosPSimHits_high[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;	
        }

    }


    m_nevents++;
}


// ------------ method called once each job just after ending the event loop  ------------
void BIBAnalyzer::endJob() {
  /*  
    // close hd5
    delete[] b;
    delete[] buffer;
    H5Sclose(mem_space);
    H5Dclose(dset);
    H5Fclose(file);
*/
    std::cout << "ptypes dict" << std::endl;
    for (auto const& pair: ptypes_dict) {
        std::cout << "{" << pair.first << ": " << pair.second << "}" << std::endl;
    }
/*	
    std::cout << "charge dict" << std::endl;
    int i = 0;
    for (auto const& pair: channel_charge_dict) {
	if (i==100){
	    break;
	}
	i++;
        std::cout << "{X=" << std::get<0>(pair.first) << " Y=" << std::get<1>(pair.first)  << ": " << pair.second << "}" << std::endl;
    }
    std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
*/
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BIBAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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


edm::DetSetVector<PixelDigiSimLink>::const_iterator BIBAnalyzer::findSimLinkDetSet(unsigned int thedetid) {
    ////basic template
    edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDS = simlinks->find(thedetid);
    return simLinkDS;
}

// ----------- method to get a unique list of simtrackids for the pixels in a cluster ------------
/*
std::map<std::tuple<float, float>, float>& BIBAnalyzer::getSimTrackId(edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter,
 						 edm::DetSet<PixelDigi>::const_iterator digit,
						 //edmNew::DetSet<SiPixelCluster>::const_iterator cluster,
						 std::map<std::tuple<float, float>, float>& channel_charge_dict,
						 unsigned int rawid,
						 bool print) {
    //int size = cluster->size();

    // loop over pixels in cluster
    //for (int i = 0; i < size; i++) {

	// get a single pixel and its channel
   // SiPixelCluster::Pixel pix = cluster->pixel(i);
    unsigned int digitChannel = PixelDigi::pixelToChannel(digit->row(), digit->column());
    //std::cout << "pixeldigi x: " << digit->row() << " pixeldigi y: " << digit->column() << std::endl;
	
	// get pixel charge in adc counts
    //int digitCharge = digit->adc();

	// only the simlinkdsv object has simtrackid info
    //if (simLinkDSViter != simlinks->end()) {
		
	// loop through simlink detset vector to find matching simlink matching by CHANNEL (unique)
	for (edm::DetSet<PixelDigiSimLink>::const_iterator it = simLinkDSViter->data.begin(); it != simLinkDSViter->data.end(); it++) {
            
	    //if simlink channel id matches digi channel id
	    if (digitChannel == it->channel()) {
			
		//store pixel charge and channel in dict
	        //channel_charge_dict[digitChannel] = digitCharge;
 //               if (print){
          //          std::cout << "Channel: " << digitChannel << " SimTrack ID: " << it->SimTrackId() << " Geomdetunit: " << rawid << " ADC charge: " << digitCharge << std::endl;
   //             }
	    }
        }
    }
    //}
    //if(simTrackIds.size() != 1){
    //std::cout << "WARNING: have more than 1 simTrackId for this cluster! " << std::endl;
    //return 0;
    //}
    return channel_charge_dict;
}
*/
DEFINE_FWK_MODULE(BIBAnalyzer);
