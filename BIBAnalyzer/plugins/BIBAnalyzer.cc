// -*- C++ -*-
//
// Package:    BRIL_ITsim/BIBAnalyzer
// Class:      BIBAnalyzer
//
/**\class BIBAnalyzer BIBAnalyzer.cc BRIL_ITsim/BIBAnalyzer/plugins/BIBAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
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

#include <TH2F.h>
#include <TTree.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

// a struct to hold the residuals for each matched cluster
struct Residual {

    double dx;
    double dy;
    double dr;

    Residual(double x, double y) : dx(x), dy(y) {
        dr = sqrt(pow(dx, 2) + pow(dy, 2));
    }

    void print() {
        std::cout << "(dx: " << dx << " dy: " << dy << " dr: " << dr << ") ";
    }

};

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
    edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> m_tokenSimLinks;
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis;
    edm::EDGetTokenT<std::vector<PSimHit>> m_tokenSimHits_low;
    edm::EDGetTokenT<std::vector<PSimHit>> m_tokenSimHits_high;

    // the pointers to geometry, topology and clusters
    // these are members so all functions can access them without passing as argument
    const TrackerTopology* tTopo = NULL;
    const TrackerGeometry* tkGeom = NULL;
    const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19
    const std::vector<PSimHit>* simhits_low = NULL;
    const std::vector<PSimHit>* simhits_high = NULL;

    //max bins of Counting histogram
    
    uint32_t m_nEBins;
    uint32_t m_maxEBin;
 
    //speed of light
    const double c = 1; //3e8;

    //array of TH2F for psimhits per disk per ring
    TH2F* m_histosPSimHits_low[8][5]; // 8 discs, 5 rings in each disc
    TH2F* m_histosPSimHits_high[8][5]; // 8 discs, 5 rings in each disc

    //histogram for time of flight vs Z coordinate
    TH2F* m_histoTofZ_low;
    TH2F* m_histoTofZ_high;
    
     //charge of hits per ring per module
    std::map<int, int> ptypes_dict;

    //1D histo for E spectra total and per disc 
    /*
    * photon 22
    * electron/positron +/-11
    * muon +/-13
    * neutron/antineutron +/-2112
    * proton/antiproton +/-2212
    * kaon +/-321
    * pion +/-211
    */
    TH1F* m_spectrumAll[9];
    TH1F* m_spectrumPhot[9];
    TH1F* m_spectrumElPM[9];
    TH1F* m_spectrumMuPM[9];
    TH1F* m_spectrumNeut[9];
    TH1F* m_spectrumProt[9];
    TH1F* m_spectrumKaPM[9];
    TH1F* m_spectrumPiPM[9];

    //2D histos for global X-Y distributions total and per disc
    TH2F* m_histoXYAll[9];
    TH2F* m_histoXYPhot[9];
    TH2F* m_histoXYElPM[9];
    TH2F* m_histoXYMuPM[9];
    TH2F* m_histoXYNeut[9];
    TH2F* m_histoXYProt[9];
    TH2F* m_histoXYKaPM[9];
    TH2F* m_histoXYPiPM[9];

    //cuts for the coincidence
    double m_dx;
    double m_dy;
    double m_dz;

    //event counter
    uint32_t m_nevents;

    // --
    // Variables for cluster parameterization studies
    bool m_storeClusterTree;
    TFile *outFileCluster;
    TTree *outTreeCluster;
    double CluX;
    double CluY;
    double CluZ;
    double CluArea;
    double CluSize;
    double CluTheta;
    double CluPhi;
    double CluCharge;
    unsigned int CluNum;
    unsigned int CluMerge;
    unsigned int mergeClu;

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
        : m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis")))
        , m_tokenSimHits_low(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_low")))
        , m_tokenSimHits_high(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_high")))
	, m_nEBins(iConfig.getUntrackedParameter<uint32_t>("nEBins")) 
        , m_maxEBin(iConfig.getUntrackedParameter<uint32_t>("maxEBin"))
        , m_dx(iConfig.getParameter<double>("dx_cut"))
        , m_dy(iConfig.getParameter<double>("dy_cut"))
        , m_dz(iConfig.getParameter<double>("dz_cut"))
        , m_storeClusterTree(iConfig.getUntrackedParameter<bool>("storeClusterTree")) {
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
    
   //loop over discs
    for (unsigned int i = 0; i < 8; i++) {

        // disc ID if [-4, 1] for -Z side and [1, 4] for +Z side
        int disk = (i < 4) ? i - 4 : i - 3;
        
        //per disc E and XY distributions per particle type and for everything
        m_spectrumAll[i] = pt.make<TH1F>(std::string("E spectrum of all particles on disk ") + disk, std::string("E spectrum of all particles on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumPhot[i] = pt.make<TH1F>(std::string("E spectrum of photons on disk ") + disk, std::string("E spectrum of photons on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumElPM[i] = pt.make<TH1F>(std::string("E spectrum of e+/- on disk ") + disk, std::string("E spectrum of e+/- on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumMuPM[i] = pt.make<TH1F>(std::string("E spectrum of muons+/- on disk ") + disk, std::string("E spectrum of muons+/- on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumNeut[i] = pt.make<TH1F>(std::string("E spectrum of neutrons on disk ") + disk, std::string("E spectrum of neutrons on disk ") + disk, m_nEBins, 0, m_maxEBin);
	m_spectrumProt[i] = pt.make<TH1F>(std::string("E spectrum of protons on disk ") + disk, std::string("E spectrum of protons on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumKaPM[i] = pt.make<TH1F>(std::string("E spectrum of kaons+/- on disk ") + disk, std::string("E spectrum of kaons+/- on disk ") + disk, m_nEBins, 0, m_maxEBin);
        m_spectrumPiPM[i] = pt.make<TH1F>(std::string("E spectrum of pions+/- on disk ") + disk, std::string("E spectrum of pions+/- on disk ") + disk, m_nEBins, 0, m_maxEBin);

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

	//TH2F 2D array
        for (unsigned int r = 0; r < 5; r++) {

               //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
               m_histosPSimHits_low[i][r] = td.make<TH2F>(std::string("Charge - time of flight (low) distribution for Disc ") + disk + " and Ring " + r+1 + ";Time of flight;Total induced charge",
							  std::string("Charge - time of flight (low) distribution for Disc ") + disk + " and Ring " + r+1,
						 	  100, -50, 50, 100, 0, 1e-3);
	       m_histosPSimHits_high[i][r] = td.make<TH2F>(std::string("Charge - time of flight (high) distribution for Disc ") + disk + " and Ring " + r+1 + ";Time of flight;Total induced charge", 
                                                           std::string("Charge - time of flight (high) distribution for Disc ") + disk + " and Ring " + r+1, 
                                                           550, -50, 500, 100, 0, 1e-3);
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

    // global coord. vs tof
    m_histoTofZ_low = td.make<TH2F>("ToF (low) vs Z", "Time of flight vs. global Z coordinate", 6000, -300, 300, 1000, -50, 50);
    m_histoTofZ_high = td.make<TH2F>("ToF (high) vs Z", "Time of flight vs. global Z coordinate", 6000, -300, 300, 5000, 0, 500);
  // ---------------------------------------------------

    if (m_storeClusterTree) {

        outFileCluster = new TFile("Cluster.root","RECREATE");
        outFileCluster->cd();
        outTreeCluster = new TTree("cluster_tree","cluster");

        outTreeCluster->Branch("CluX", &CluX);
        outTreeCluster->Branch("CluY", &CluY);
        outTreeCluster->Branch("CluZ", &CluZ);
        outTreeCluster->Branch("CluTheta", &CluTheta);
        outTreeCluster->Branch("CluPhi", &CluPhi);
        outTreeCluster->Branch("CluCharge", &CluCharge);
        outTreeCluster->Branch("CluArea", &CluArea);
        outTreeCluster->Branch("CluSize", &CluSize);
        outTreeCluster->Branch("CluMerge", &CluMerge);
        outTreeCluster->Branch("CluNum", &CluNum);
    }
}

// ------------ method called for each event  ------------
void BIBAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    //get the digis - COB 26.02.19
    edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
    iEvent.getByToken(m_tokenDigis, tdigis);

    //get the simhits - COB 30.11.20
    edm::Handle<std::vector<PSimHit>> tsimhits_low;
    iEvent.getByToken(m_tokenSimHits_low, tsimhits_low);
    edm::Handle<std::vector<PSimHit>> tsimhits_high;
    iEvent.getByToken(m_tokenSimHits_high, tsimhits_high);   

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
            int hist_id = -1;
            unsigned int ring_id = ring - 1;
            //unsigned int module_id = module - 1;
        
	    //this is a TEPX hit on side 1 (-Z)
            if (side == 1) {
                hist_id = 12 - disk; // goes from 0 to 3
            }
            
	    //this is a TEPX hit on side 2 (+Z)
            else if (side == 2) {
                hist_id = 4 + disk - 9; // goes from 4 to 7
            }

            //find the geometry of the module associated to the digi
            const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
            if (!geomDetUnit)
                continue;

	    // get global coordinates of PSimHit
	    Point3DBase<float, LocalTag> local_coords = DSVit->localPosition();
	    
	    // PixelGeomDetUnit is derived from GeomDet, transform local to global coords
	    GlobalPoint global_coords = geomDetUnit->toGlobal(local_coords);
      
            //get particle type
            int ptype = DSVit->particleType();
	    if (ptypes_dict.find(ptype) == ptypes_dict.end()) {
 	    	ptypes_dict[ptype] = 1;
		std::cout << "new key" << std::endl;
	    }
	    else{
		ptypes_dict[ptype] += 1;
		std::cout << "updating existing key" << std::endl;
	    }

 	    //debug info
            //std::cout << "---- Hit info:" << std::endl;
	    //std::cout << "hit tof: " << DSVit->timeOfFlight() << std::endl;
	    //std::cout << "global position: " << global_coords.z() << std::endl;		
	    //std::cout << "disc id: " <<	hist_id << std::endl;
	    //std::cout << "particle type: " << DSVit->particleType() << std::endl;

	    // fill up particle type dependent histos
	    switch (std::abs(ptype)){
		case 22: // photon
                    m_spectrumPhot[8]->Fill(c*DSVit->pabs()); 
	            m_histoXYPhot[8]->Fill(global_coords.x(), global_coords.y()); 
        	    m_spectrumPhot[hist_id]->Fill(c*DSVit->pabs());
            	    m_histoXYPhot[hist_id]->Fill(global_coords.x(), global_coords.y());
		    break;
                 case 11: // electron/positron   
                    m_spectrumElPM[8]->Fill(c*DSVit->pabs());
                    m_histoXYElPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumElPM[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYElPM[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break; 
		 case 13: // muon+/-
                    m_spectrumMuPM[8]->Fill(c*DSVit->pabs());
                    m_histoXYMuPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumMuPM[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYMuPM[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break; 
		case 2112: // neutron/antineutron 
                    m_spectrumNeut[8]->Fill(c*DSVit->pabs());
                    m_histoXYNeut[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumNeut[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYNeut[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break; 
		case 2212: // proton/antiproton
                    m_spectrumProt[8]->Fill(c*DSVit->pabs());
                    m_histoXYProt[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumProt[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYProt[hist_id]->Fill(global_coords.x(), global_coords.y());
	    	    break;
                case 321: // kaon+/-
                    m_spectrumKaPM[8]->Fill(c*DSVit->pabs()); 
                    m_histoXYKaPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumKaPM[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYKaPM[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break;
                case 211: // pion+/-
                    m_spectrumPiPM[8]->Fill(c*DSVit->pabs()); 
                    m_histoXYPiPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumPiPM[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYPiPM[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break; 
		}
		
	    // fill up all particle histos	
	    m_spectrumAll[8]->Fill(c*DSVit->pabs());	    
            m_histoXYAll[8]->Fill(global_coords.x(), global_coords.y());	
	    m_spectrumAll[hist_id]->Fill(c*DSVit->pabs());
            m_histoXYAll[hist_id]->Fill(global_coords.x(), global_coords.y());

	    // fill up tof histos
            m_histosPSimHits_low[hist_id][ring_id]->Fill(DSVit->timeOfFlight(), DSVit->energyLoss());
	    m_histoTofZ_low->Fill(global_coords.z(), DSVit->timeOfFlight());
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
            int hist_id = -1;
            unsigned int ring_id = ring - 1;
            //unsigned int module_id = module - 1;
        
	    //this is a TEPX hit on side 1 (-Z)
            if (side == 1) {
                hist_id = 12 - disk; // goes from 0 to 3
            }
            
	    //this is a TEPX hit on side 2 (+Z)
            else if (side == 2) {
                hist_id = 4 + disk - 9; // goes from 4 to 7
            }

            //find the geometry of the module associated to the digi
            const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
            if (!geomDetUnit)
                continue;

	    // get global coordinates of PSimHit
	    Point3DBase<float, LocalTag> local_coords = DSVit->localPosition();
	    
	    // PixelGeomDetUnit is derived from GeomDet, transform local to global coords
	    GlobalPoint global_coords = geomDetUnit->toGlobal(local_coords);
      
            //get particle type
            int ptype = DSVit->particleType();
	    if (ptypes_dict.find(ptype) == ptypes_dict.end()) {
 	    	ptypes_dict[ptype] = 1;
		std::cout << "new key" << std::endl;
	    }
	    else{
		ptypes_dict[ptype] += 1;
		std::cout << "updating existing key" << std::endl;
	    }

 	    //debug info
            //std::cout << "---- Hit info:" << std::endl;
	    //std::cout << "hit tof: " << DSVit->timeOfFlight() << std::endl;
	    //std::cout << "global position: " << global_coords.z() << std::endl;		
	    //std::cout << "disc id: " <<	hist_id << std::endl;
	    //std::cout << "particle type: " << DSVit->particleType() << std::endl;
            //std::cout << "momentum: " << DSVit->pabs() << std::endl;

	    // fill up particle type dependent histos
	    switch (std::abs(ptype)){
		case 22: // photon
                    m_spectrumPhot[8]->Fill(c*DSVit->pabs()); 
	            m_histoXYPhot[8]->Fill(global_coords.x(), global_coords.y()); 
        	    m_spectrumPhot[hist_id]->Fill(c*DSVit->pabs());
            	    m_histoXYPhot[hist_id]->Fill(global_coords.x(), global_coords.y());
		    break;
                 case 11: // electron/positron   
                    m_spectrumElPM[8]->Fill(c*DSVit->pabs());
                    m_histoXYElPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumElPM[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYElPM[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break; 
		 case 13: // muon+/-
                    m_spectrumMuPM[8]->Fill(c*DSVit->pabs());
                    m_histoXYMuPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumMuPM[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYMuPM[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break; 
		case 2112: // neutron/antineutron 
                    m_spectrumNeut[8]->Fill(c*DSVit->pabs());
                    m_histoXYNeut[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumNeut[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYNeut[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break; 
		case 2212: // proton/antiproton
                    m_spectrumProt[8]->Fill(c*DSVit->pabs());
                    m_histoXYProt[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumProt[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYProt[hist_id]->Fill(global_coords.x(), global_coords.y());
	    	    break;
                case 321: // kaon+/-
                    m_spectrumKaPM[8]->Fill(c*DSVit->pabs()); 
                    m_histoXYKaPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumKaPM[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYKaPM[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break;
                case 211: // pion+/-
                    m_spectrumPiPM[8]->Fill(c*DSVit->pabs()); 
                    m_histoXYPiPM[8]->Fill(global_coords.x(), global_coords.y());
                    m_spectrumPiPM[hist_id]->Fill(c*DSVit->pabs());
                    m_histoXYPiPM[hist_id]->Fill(global_coords.x(), global_coords.y());
                    break; 
		}
		
	    // fill up all particle histos	
	    m_spectrumAll[8]->Fill(c*DSVit->pabs());	    
            m_histoXYAll[8]->Fill(global_coords.x(), global_coords.y());	
	    m_spectrumAll[hist_id]->Fill(c*DSVit->pabs());
            m_histoXYAll[hist_id]->Fill(global_coords.x(), global_coords.y());

            // fill up tof histos	
            m_histosPSimHits_high[hist_id][ring_id]->Fill(DSVit->timeOfFlight(), DSVit->energyLoss());	
	    m_histoTofZ_high->Fill(global_coords.z(), DSVit->timeOfFlight());
        }

    }


    m_nevents++;
}

// ------------ method called once each job just after ending the event loop  ------------
void BIBAnalyzer::endJob() {

    if (m_storeClusterTree) {
        outFileCluster->Write("",TObject::kOverwrite);
        outFileCluster->Close();
    }
    std::cout << "ptypes dict" << std::endl;
    for (auto const& pair: ptypes_dict) {
        std::cout << "{" << pair.first << ": " << pair.second << "}" << std::endl;
    }
    std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
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


//----------
//Adding function to construct DetId from PXFDetId
//COB - 22.May.2019
uint32_t BIBAnalyzer::getModuleID(bool isTEPX, unsigned int side, unsigned int disk, unsigned int ring) {

    //std::cout << "isTEPX " << isTEPX << " side " << side << " disk " << disk << " ring " << ring << std::endl;

    uint32_t modid = -999;

    if (isTEPX) {
        if (side==1) {
            if (disk==9) {
                if (ring==2) {
                    modid = 346301444;
                } else if (ring==3) {
                    modid = 346305540;
                } else if (ring==4) {
                    modid = 346309636;
                } else if (ring==5) {
                    modid = 346313732;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << " side " << side << "!" << std::endl;
                    return modid;
                }
            } else if (disk==10) {
                if (ring==2) {
                    modid = 346563588;
                } else if (ring==3) {
                    modid = 346567684;
                } else if (ring==4) {
                    modid = 346571780;
                } else if (ring==5) {
                    modid = 346575876;
                } else {
                   std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                   return modid;
                }
            } else if (disk==11) {
                if (ring==2) {
                    modid = 346825732;
                } else if (ring==3) {
                    modid = 346829828;
                } else if (ring==4) {
                    modid = 346833924;
                } else if (ring==5) {
                    modid = 346838020;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==12) {
                if (ring==2) {
                    modid = 347087876;
                } else if (ring==3) {
                    modid = 347091972;
                } else if (ring==4) {
                    modid = 347096068;
                } else if (ring==5) {
                    modid = 347100164;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else {
                std::cout << "Non-existent disk for TEPX!" << std::endl;
                return modid;
            }
        } else if (side==2) {
            if (disk==9) {
                if (ring==2) {
                    modid = 354690052;
                } else if (ring==3) {
                    modid = 354694148;
                } else if (ring==4) {
                    modid = 354698244;
                } else if (ring==5) {
                    modid = 354702340;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==10) {
                if (ring==2) {
                    modid = 354952196;
                } else if (ring==3) {
                    modid = 354956292;
                } else if (ring==4) {
                    modid = 354960388;
                } else if (ring==5) {
                    modid = 354964484;
                } else {
                   std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                   return modid;
                }
            } else if (disk==11) {
                if (ring==2) {
                    modid = 355214340;
                } else if (ring==3) {
                    modid = 355218436;
                } else if (ring==4) {
                    modid = 355222532;
                } else if (ring==5) {
                    modid = 355226628;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==12) {
                if (ring==2) {
                    modid = 355476484;
                } else if (ring==3) {
                    modid = 355480580;
                } else if (ring==4) {
                    modid = 355484676;
                } else if (ring==5) {
                    modid = 355488772;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else {
                std::cout << "Non-existent disk for TEPX!" << std::endl;
                return modid;
            }
        } else {
            std::cout << "Non-existent side!" << std::endl;
            return modid;
        }
    } else {  //TFPX

        if (side==1) {
            if (disk==1) {
                if (ring==2) {
                    modid = 344204292;
                } else if (ring==3) {
                    modid = 344208388;
                } else if (ring==4) {
                    modid = 344212484;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==2) {
                if (ring==2) {
                    modid = 344466436;
                } else if (ring==3) {
                    modid = 344470532;
                } else if (ring==4) {
                    modid = 344474628;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==3) {
                if (ring==2) {
                    modid = 344728580;
                } else if (ring==3) {
                    modid = 344732676;
                } else if (ring==4) {
                    modid = 344736772;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==4) {
                if (ring==2) {
                    modid = 344990724;
                } else if (ring==3) {
                    modid = 344994820;
                } else if (ring==4) {
                    modid = 344998916;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==5) {
                if (ring==2) {
                    modid = 345252868;
                } else if (ring==3) {
                    modid = 345256964;
                } else if (ring==4) {
                    modid = 345261060;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==6) {
                if (ring==2) {
                    modid = 345515012;
                } else if (ring==3) {
                    modid = 345519108;
                } else if (ring==4) {
                    modid = 345523204;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==7) {
                if (ring==2) {
                    modid = 345777156;
                } else if (ring==3) {
                    modid = 345781252;
                } else if (ring==4) {
                    modid = 345785348;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==8) {
                if (ring==2) {
                    modid = 346039300;
                } else if (ring==3) {
                    modid = 346043396;
                } else if (ring==4) {
                    modid = 346047492;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else {
                std::cout << "Non-existent disk for TFPX!" << std::endl;
                return modid;
            }
        } else if (side==2) {
            if (disk==1) {
                if (ring==2) {
                    modid = 352592900;
                } else if (ring==3) { 
                    modid = 352596996;
                } else if (ring==4) {
                    modid = 352601092;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==2) {
                if (ring==2) {
                    modid = 352855044;
                } else if (ring==3) {
                    modid = 352859140;
                } else if (ring==4) {
                    modid = 352863236;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==3) {
                if (ring==2) {
                    modid = 353117188;
                } else if (ring==3) {
                    modid = 353121284;
                } else if (ring==4) {
                    modid = 353125380;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==4) {
                if (ring==2) {
                    modid = 353379332;
                } else if (ring==3) {
                    modid = 353383428;
                } else if (ring==4) {
                    modid = 353387524;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==5) {
                if (ring==2) {
                    modid = 353641476;
                } else if (ring==3) {
                    modid = 353645572;
                } else if (ring==4) {
                    modid = 353649668;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==6) {
                if (ring==2) {
                    modid = 353903620;
                } else if (ring==3) {
                    modid = 353907716;
                } else if (ring==4) {
                    modid = 353911812;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==7) {
                if (ring==2) {
                    modid = 354165764;
                } else if (ring==3) {
                    modid = 354169860;
                } else if (ring==4) {
                    modid = 354173956;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else if (disk==8) {
                if (ring==2) {
                    modid = 354427908;
                } else if (ring==3) {
                    modid = 354432004;
                } else if (ring==4) {
                    modid = 354436100;
                } else {
                    std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
                    return modid;
                }
            } else {
                std::cout << "Non-existent disk for TFPX!" << std::endl;
                return modid;
            }                            
        } else {
            std::cout << "Non-existent side!" << std::endl;
            return modid;
        }

    }

    return modid;

}

DEFINE_FWK_MODULE(BIBAnalyzer);
