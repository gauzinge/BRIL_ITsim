// -*- C++ -*-
//
// Package:    BRIL_ITsim/ITclusterAnalyzer
// Class:      ITclusterAnalyzer
//
/**\class ITclusterAnalyzer ITclusterAnalyzer.cc BRIL_ITsim/ITclusterAnalyzer/plugins/ITclusterAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Georg Auzinger
//         Created:  Wed, 05 Dec 2018 15:10:16 GMT
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
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TTree.h>
#include <TH2F.h>

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

    Residual(double x, double y)
        : dx(x)
        , dy(y)
    {
        dr = sqrt(pow(dx, 2) + pow(dy, 2));
    }
    void print()
    {
        std::cout << "(dx: " << dx << " dy: " << dy << " dr: " << dr << ") ";
    }
};

class ITclusterAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit ITclusterAnalyzer(const edm::ParameterSet&);
    ~ITclusterAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //bool findCoincidence(DetId, Global3DPoint, bool);
    bool findCoincidence2x(DetId, Global3DPoint, unsigned int&, edmNew::DetSet<SiPixelCluster>::const_iterator&);
    bool findCoincidence3x(DetId, Global3DPoint, unsigned int&, edmNew::DetSet<SiPixelCluster>::const_iterator&);
    edm::DetSetVector<PixelDigiSimLink>::const_iterator findSimLinkDetSet(unsigned int thedetid);
    bool isMerged(edm::DetSetVector<PixelDigiSimLink>::const_iterator, edmNew::DetSet<SiPixelCluster>::const_iterator, bool print);
    std::set<unsigned int> getSimTrackId(edm::DetSetVector<PixelDigiSimLink>::const_iterator, edmNew::DetSet<SiPixelCluster>::const_iterator, bool print);
    bool areSameSimTrackId(std::set<unsigned int> first, std::set<unsigned int> second, std::set<unsigned int>&);

    // ----------member data ---------------------------
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> m_tokenClusters;
    edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> m_tokenSimLinks;
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis;

    // the pointers to geometry, topology and clusters
    // these are members so all functions can access them without passing as argument
    const TrackerTopology* tTopo = NULL;
    const TrackerGeometry* tkGeom = NULL;
    const edmNew::DetSetVector<SiPixelCluster>* clusters = NULL;
    const edm::DetSetVector<PixelDigiSimLink>* simlinks = NULL;
    const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19

    //max bins of Counting histogram
    uint32_t m_maxBin;
    //flag for checking coincidences
    bool m_docoincidence;

    //array of TH2F for clusters per disk per ring
    TH2F* m_diskHistosCluster[8];
    //array of TH2F for hits per disk per ring
    TH2F* m_diskHistosHits[8];

    //tracker maps for clusters
    TH2F* m_trackerLayoutClustersZR;
    TH2F* m_trackerLayoutClustersYX;
    //tracker maps for hits
    TH2F* m_trackerLayoutHitsZR;
    TH2F* m_trackerLayoutHitsYX;

    //array of TH2F for 2xcoinc per disk per ring
    //first all coincidences
    TH2F* m_diskHistos2x[8];
    //and the real ones
    TH2F* m_diskHistos2xreal[8];
    //tracker maps for 2xcoinc
    TH2F* m_trackerLayout2xZR;
    TH2F* m_trackerLayout2xYX;

    //array of TH2F for 3xcoinc per disk per ring
    //first all coincidences
    TH2F* m_diskHistos3x[8];
    //and the real ones
    TH2F* m_diskHistos3xreal[8];
    //tracker maps for 3xcoinc
    TH2F* m_trackerLayout3xZR;
    TH2F* m_trackerLayout3xYX;

    //simple residual histograms for the cuts
    TH1F* m_residualX;
    TH1F* m_residualY;
    TH1F* m_residualR;

    //the number of clusters per module
    TH1F* m_nClusters;
    TH1F* m_nClustersV1;
    TH1F* m_nHits;

   
    //Added by Joseph Cordero Mercado
    //---- cluster size
    TH1F* m_HitX;
    TH1F* m_HitY;
    TH1F* m_HitZ;

    TH1F* m_CluX;
    TH1F* m_CluY;
    TH1F* m_CluZ;

    TH1F* m_CluPhi;
    TH1F* m_CluTheta;

    TH1F* m_mergeClus;

    TH2F* m_ClusterXsize[8];
    TH2F* m_ClusterYsize[8];
    TH2F* m_ClusterArea[8];


    TFile *outFileHit;
    TTree *outTreeHit;
    double HitX;
    double HitY;
    double HitZ;
    double HitTheta;
    double HitPhi;
    unsigned int HitNum;
    unsigned int mergeHit;

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
    //cuts for the coincidence
    double m_dx;
    double m_dy;
    double m_dz;

    //event counter
    uint32_t m_nevents;
    //coincidence counter
    uint32_t m_total2xcoincidences;
    uint32_t m_fake2xcoincidences;
    uint32_t m_total3xcoincidences;
    uint32_t m_fake3xcoincidences;
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
ITclusterAnalyzer::ITclusterAnalyzer(const edm::ParameterSet& iConfig)
	: //m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>> ("clusters"))
		m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters")))
		, m_tokenSimLinks(consumes<edm::DetSetVector<PixelDigiSimLink>>(iConfig.getParameter<edm::InputTag>("simlinks")))
		, m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis"))) //adding digis variable - COB 26.02.19
		, m_maxBin(iConfig.getUntrackedParameter<uint32_t>("maxBin"))
		, m_docoincidence(iConfig.getUntrackedParameter<bool>("docoincidence"))
		, m_dx(iConfig.getParameter<double>("dx_cut"))
		, m_dy(iConfig.getParameter<double>("dy_cut"))
		, m_dz(iConfig.getParameter<double>("dz_cut"))
{
	//now do what ever initialization is needed
	m_nevents = 0;
	m_total2xcoincidences = 0;
	m_fake2xcoincidences = 0;
	m_total3xcoincidences = 0;
	m_fake3xcoincidences = 0;
}

ITclusterAnalyzer::~ITclusterAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void ITclusterAnalyzer::beginJob()
{

	outFileHit = new TFile("Hits.root","RECREATE");
	outFileHit->cd();
	outTreeHit = new TTree("hits_tree","hits");

	outTreeHit->Branch("HitX"      , &HitX);
	outTreeHit->Branch("HitY"      , &HitY);
	outTreeHit->Branch("HitZ"      , &HitZ);
	outTreeHit->Branch("HitTheta"  , &HitTheta);
	outTreeHit->Branch("HitPhi"    , &HitPhi);
	outTreeHit->Branch("HitNum"    , &HitNum);

	outFileCluster = new TFile("Cluster.root","RECREATE");
	outFileCluster->cd();
	outTreeCluster = new TTree("cluster_tree","cluster");

	outTreeCluster->Branch("CluX"      , &CluX);
	outTreeCluster->Branch("CluY"      , &CluY);
	outTreeCluster->Branch("CluZ"      , &CluZ);
	outTreeCluster->Branch("CluTheta"  , &CluTheta);
	outTreeCluster->Branch("CluPhi"    , &CluPhi);
	outTreeCluster->Branch("CluCharge" , &CluCharge);
	outTreeCluster->Branch("CluArea"   , &CluArea);
	outTreeCluster->Branch("CluSize"   , &CluSize);
	outTreeCluster->Branch("CluMerge"  , &CluMerge);
	outTreeCluster->Branch("CluNum"    , &CluNum);
	//m_nClustersV1 = td.make<TH1F>("Number of Clusters per module per event", "# of Clusters;# of Clusters; Occurence", 500, 0, 500);
	//////////////////////////////////////////////////////
	edm::Service<TFileService> fs;
	//////////////////////////////////////////////////////
	fs->file().cd("/");
	TFileDirectory td = fs->mkdir("Residuals");

	m_residualX = td.make<TH1F>("ResidualsX", "ResidualsX;deltaX;counts", 1000, 0, 1);
	m_residualY = td.make<TH1F>("ResidualsY", "ResidualsY;deltaY;counts", 1000, 0, 1);
	m_residualR = td.make<TH1F>("ResidualsR", "ResidualsR;deltaR;counts", 1000, 0, 1);

	fs->file().cd("/");
	td = fs->mkdir("perModule");
	m_nClusters = td.make<TH1F>("Number of Clusters per module per event", "# of Clusters;# of Clusters; Occurence", 500, 0, 500);
        m_nHits = td.make<TH1F>("Number of Hits per module per event", "# of Hits; # of Hits; Occurrences", 500, 0, 500);

	//////////////////////////////////////////////////////
	fs->file().cd("/");
	td = fs->mkdir("Clusters");

	m_CluX = td.make<TH1F>("CluX", "CluX;X;counts", 1000, -50, 50);
	m_CluY = td.make<TH1F>("CluY", "CluY;Y;counts", 1000, -50, 50);
	m_CluZ = td.make<TH1F>("CluZ", "CluR;Z;counts", 6000, -300, 300);

	m_CluPhi = td.make<TH1F>("CluPhi", "CluPhi;Phi;counts", 1000, -4, 4);
	m_CluTheta = td.make<TH1F>("CluTheta", "CluTheta;Theta;counts", 1000, -4, 4);

	m_mergeClus = td.make<TH1F>("mergeClu","mergeClu;merge;# of merged events",5,0,5);	

	//now lets create the histograms
	for (unsigned int i = 0; i < 8; i++) {
		int disk = (i < 4) ? i - 4 : i - 3;
		std::stringstream histoname;
		histoname << "Number of clusters for Disk " << disk << ";Ring;# of Clusters per event";
		std::stringstream histotitle;
		histotitle << "Number of clusters for Disk " << disk;
		//name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
		m_diskHistosCluster[i] = td.make<TH2F>(histotitle.str().c_str(), histoname.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);

		
		// Added by Joseph Cordero Mercado
		std::stringstream histonameSizeX;
		histonameSizeX << "SizeX of clusters for Disk " << disk << ";Ring;SizeX of Clusters per event";
		std::stringstream histotitleSizeX;
		histotitleSizeX << "SizeX of clusters for Disk " << disk;
		m_ClusterXsize[i] = td.make<TH2F>(histotitleSizeX.str().c_str(), histonameSizeX.str().c_str(),5,0.5,5.5,15,0,15);

		std::stringstream histonameSizeY;
		histonameSizeY << "SizeY of clusters for Disk " << disk << ";Ring;SizeY of Clusters per event";
		std::stringstream histotitleSizeY;
		histotitleSizeY << "SizeY of clusters for Disk " << disk;
		m_ClusterYsize[i] = td.make<TH2F>(histotitleSizeY.str().c_str(), histonameSizeY.str().c_str(),5,0.5,5.5,15,0,15);

		std::stringstream histonameArea;
		histonameArea << "Area of clusters for Disk " << disk << ";Ring;Area of Clusters per event";
		std::stringstream histotitleArea;
		histotitleArea << "Area of clusters for Disk " << disk;
		m_ClusterArea[i]  = td.make<TH2F>(histotitleArea.str().c_str(), histonameArea.str().c_str(),5,0.5,5.5,100,0,100);
	}
	m_trackerLayoutClustersZR = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 600, 0.0, 30.0);
	m_trackerLayoutClustersYX = td.make<TH2F>("XVsY", "x vs. y position", 1000, -50.0, 50.0, 1000, -50.0, 50.0);

	//////////////////////////////////////////////////////
        fs->file().cd("/");
        td = fs->mkdir("Hits");

	m_HitX = td.make<TH1F>("HitX", "HitX;X;counts", 1000, -50, 50);
	m_HitY = td.make<TH1F>("HitY", "HitY;Y;counts", 1000, -50, 50);
	m_HitZ = td.make<TH1F>("HitZ", "HitR;Z;counts", 6000, -300, 300);
        //histograms
        for (unsigned int i = 0; i < 8; i++) {
            int disk = (i < 4) ? i - 4 : i - 3;
            std::stringstream histoname;
            histoname << "Number of hits for Disk " << disk << ";Ring;# of Hits per event";
            std::stringstream histotitle;
            histotitle << "Number of hits for Disk " << disk;
            m_diskHistosHits[i] = td.make<TH2F>(histotitle.str().c_str(), histoname.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);
        }

        m_trackerLayoutHitsZR = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 600, 0.0, 30.0);
        m_trackerLayoutHitsYX = td.make<TH2F>("XVsY", "x vs. y position", 1000, -50.0, 50.0, 1000, -50.0, 50.0);

	//////////////////////////////////////////////////////
	if (m_docoincidence) {
		fs->file().cd("/");
		td = fs->mkdir("2xCoincidences");
		//now lets create the histograms
		for (unsigned int i = 0; i < 8; i++) {
			int disk = (i < 4) ? i - 4 : i - 3;
			std::stringstream histoname;
			histoname << "Number of 2x Coincidences for Disk " << disk << ";Ring;# of coincidences per event";
			std::stringstream histotitle;
			histotitle << "Number of 2x Coincidences for Disk " << disk;
			//name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
			m_diskHistos2x[i] = td.make<TH2F>(histotitle.str().c_str(), histoname.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);

			std::stringstream histonamereal;
			histonamereal << "Number of real 2x Coincidences for Disk " << disk << ";Ring;# of real coincidences per event";
			std::stringstream histotitlereal;
			histotitlereal << "Number of real 2x Coincidences for Disk " << disk;
			//name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
			m_diskHistos2xreal[i] = td.make<TH2F>(histotitlereal.str().c_str(), histonamereal.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);
		}
		m_trackerLayout2xZR = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 600, 0.0, 30.0);
		m_trackerLayout2xYX = td.make<TH2F>("XVsY", "x vs. y position", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
	}

	//////////////////////////////////////////////////////
	if (m_docoincidence) {
		fs->file().cd("/");
		td = fs->mkdir("3xCoincidences");
		//now lets create the histograms
		for (unsigned int i = 0; i < 8; i++) {
			int disk = (i < 4) ? i - 4 : i - 3;
			std::stringstream histoname;
			histoname << "Number of 3x Coincidences for Disk " << disk << ";Ring;# of coincidences per event";
			std::stringstream histotitle;
			histotitle << "Number of 3x Coincidences for Disk " << disk;
			//name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
			m_diskHistos3x[i] = td.make<TH2F>(histotitle.str().c_str(), histoname.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);

			std::stringstream histonamereal;
			histonamereal << "Number of real 3x Coincidences for Disk " << disk << ";Ring;# of real coincidences per event";
			std::stringstream histotitlereal;
			histotitlereal << "Number of real 3x Coincidences for Disk " << disk;
			//name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
			m_diskHistos3xreal[i] = td.make<TH2F>(histotitlereal.str().c_str(), histonamereal.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);
		}
		m_trackerLayout3xZR = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 600, 0.0, 30.0);
		m_trackerLayout3xYX = td.make<TH2F>("XVsY", "x vs. y position", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
	}
}

// ------------ method called for each event  ------------
void ITclusterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //get the digis - COB 26.02.19
    edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
    iEvent.getByToken(m_tokenDigis, tdigis);

    //get the clusters
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> tclusters;
    iEvent.getByToken(m_tokenClusters, tclusters);

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
    clusters = tclusters.product();
    simlinks = tsimlinks.product();
    digis = tdigis.product();  //pointer to digis - COB 26.02.19

    //a 2D counter array to count the number of clusters per disk and per ring
    unsigned int cluCounter[8][5];
    memset(cluCounter, 0, sizeof(cluCounter));
    //counter array for hits per disk and per ring
    unsigned int hitCounter[8][5];
    memset(hitCounter, 0, sizeof(hitCounter));
    //counter for 2x coincidences
    unsigned int x2Counter[8][5];
    memset(x2Counter, 0, sizeof(x2Counter));
    unsigned int x2Counterreal[8][5];
    memset(x2Counterreal, 0, sizeof(x2Counterreal));

    unsigned int x3Counter[8][5];
    memset(x3Counter, 0, sizeof(x3Counter));
    unsigned int x3Counterreal[8][5];
    memset(x3Counterreal, 0, sizeof(x3Counterreal));

    // Added by Joseph Cordero Mercado 
    unsigned int hitX[8][5];
    memset(hitX, 0, sizeof(hitX));
    unsigned int hitY[8][5];
    memset(hitY, 0, sizeof(hitY));
    unsigned int hitZ[8][5];
    memset(hitZ, 0, sizeof(hitZ));

    unsigned int cluX[8][5];
    memset(cluX, 0, sizeof(cluX));
    unsigned int cluY[8][5];
    memset(cluY, 0, sizeof(cluY));
    unsigned int cluZ[8][5];
    memset(cluZ, 0, sizeof(cluZ));

    unsigned int cluSizeXCounter[8][5];
    memset(cluSizeXCounter, 0, sizeof(cluSizeXCounter));
    unsigned int cluSizeYCounter[8][5];
    memset(cluSizeYCounter, 0, sizeof(cluSizeYCounter));
    unsigned int cluAreaCounter[8][5];
    memset(cluAreaCounter, 0, sizeof(cluAreaCounter));

    //-------------------------------------------------------------
    //loop over digis - COB 26.02.19
    //debugging...
    //std::cout << "--- New event ---" << std::endl;
    for (typename edm::DetSetVector<PixelDigi>::const_iterator DSVit = digis->begin(); DSVit != digis->end(); DSVit++) {

        //get the detid
        unsigned int rawid(DSVit->detId());
        DetId detId(rawid);
     
        //debugging...
        //std::cout << "DetId " << std::hex << "0x" << rawid << std::dec << " " << detId.det() << " " << detId.subdetId() << " " 
        //          << ((rawid >> 23) & 0x3) << " " << ((rawid >> 18) & 0xF) << " " << ((rawid >> 12) & 0x3F) << " " << ((rawid >> 2) & 0xFF) << std::endl;
        //std::cout << "side " << tTopo->pxfSide(detId) << " layer " << tTopo->pxfDisk(detId) << " ring " << tTopo->pxfBlade(detId) 
        //          << " module " << tTopo->pxfModule(detId) << std::endl;       

        //module type => need phase 2 pixel forward module, in endcap
        TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
        if (mType != TrackerGeometry::ModuleType::Ph2PXF && detId.subdetId() != PixelSubdetector::PixelEndcap) 
            continue;

        //obtaining side, disk, ring of module
        unsigned int side = (tTopo->pxfSide(detId));   
        unsigned int disk = (tTopo->pxfDisk(detId));
        unsigned int ring = (tTopo->pxfBlade(detId));

        //look only at TEPX modules
        if (disk < 9)
            continue; 

        //index in histogram map
        int hist_id = -1;
        unsigned int ring_id = ring - 1;
        if (side == 1) {
            //this is a TEPX- hit on side1
            hist_id = disk - 9;
        } else if (side == 2) {
            //this is a TEPX+ hit on side 2
            hist_id = 4 + disk - 9;
        }

        //find the geometry of the module associated to the digi
        const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
        if (!geomDetUnit)
            continue;

        //store the total number of digis in module
        m_nHits->Fill(DSVit->size());
	HitNum = DSVit->size();
 
        //loop over the digis in each module
        for (edm::DetSet<PixelDigi>::const_iterator digit = DSVit->begin(); digit != DSVit->end(); digit++) {

            hitCounter[hist_id][ring_id]++;

            //finding the position of the hit
            MeasurementPoint mpHit(digit->row(), digit->column());
            Local3DPoint localPosHit = geomDetUnit->topology().localPosition(mpHit);
            Global3DPoint globalPosHit = geomDetUnit->surface().toGlobal(localPosHit);

            //debugging...
            //std::cout << "---- Hit info:" << std::endl;
            //std::cout << "measurement point of hit => row " << digit->row() << " column " << digit->column() << std::endl;
            //std::cout << "local position of hit => x " << localPosHit.x() << " y " << localPosHit.y() << " z " << localPosHit.z() << std::endl;
            //std::cout << "global position of hit => x " << globalPosHit.x() << " y " << globalPosHit.y() << " z " << globalPosHit.z() << std::endl;

            //fill layout histograms
            m_trackerLayoutHitsZR->Fill(globalPosHit.z(), globalPosHit.perp());
            m_trackerLayoutHitsYX->Fill(globalPosHit.x(), globalPosHit.y());


	    HitX      = globalPosHit.x();
	    HitY      = globalPosHit.y();
	    HitZ      = globalPosHit.z();
	    HitPhi    = globalPosHit.phi();
	    HitTheta  = globalPosHit.theta();

	    //Added by Joseph Cordero 
	    m_HitX->Fill(globalPosHit.x());
	    m_HitY->Fill(globalPosHit.y());
	    m_HitZ->Fill(globalPosHit.z());

        }

    }
    //-------------------------------------------------------------

    //loop the modules in the cluster collection
    for (typename edmNew::DetSetVector<SiPixelCluster>::const_iterator DSVit = clusters->begin(); DSVit != clusters->end(); DSVit++) {
        //get the detid
        unsigned int rawid(DSVit->detId());
        DetId detId(rawid);

        //figure out the module type using the detID
        TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
        if (mType != TrackerGeometry::ModuleType::Ph2PXF && detId.subdetId() != PixelSubdetector::PixelEndcap)
            continue;
        //std::cout << "DetID " << std::hex << "0x" << rawid << std::dec << " " << detId.det() << " " << detId.subdetId() << " " << ((rawid >> 23) & 0x3) << " " << ((rawid >> 18) & 0xF) << " " << ((rawid >> 12) & 0x3F) << " " << ((rawid >> 2) & 0xFF) << std::endl;

        //find out which layer, side and ring
        unsigned int side = (tTopo->pxfSide(detId));  // values are 1 and 2 for -+Z
        unsigned int layer = (tTopo->pxfDisk(detId)); //values are 1 to 12 for disks TFPX1 to TFPX 8  and TEPX1 to TEPX 4
        unsigned int ring = (tTopo->pxfBlade(detId));

        //now make sure we are only looking at TEPX
        if (layer < 9)
            continue;

        //the index in my histogram map
        int hist_id = -1;
        unsigned int ring_id = ring - 1;

        if (side == 1) {
            //this is a TEPX- hit on side1
            hist_id = layer - 9;
        } else if (side == 2) {
            //this is a TEPX+ hit on side 2
            hist_id = 4 + layer - 9;
        }

        // Get the geomdet
        const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
        if (!geomDetUnit)
            continue;

        unsigned int nClu = 0;

        //fill the number of clusters for this module
        m_nClusters->Fill(DSVit->size());
	CluNum = DSVit->size();
        //m_nClustersV1->Fill(DSVit->size());

        //now loop the clusters for each detector
        for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit = DSVit->begin(); cluit != DSVit->end(); cluit++) {
            //increment the counters
            nClu++;
            cluCounter[hist_id][ring_id]++;
	
            cluSizeXCounter[hist_id][ring_id] += cluit->sizeX();
            cluSizeYCounter[hist_id][ring_id] += cluit->sizeY();
            cluAreaCounter [hist_id][ring_id] += (cluit->sizeY())*(cluit->sizeX());
	

           // std::cout << "------------------------------------" << std::endl;
           // std::cout << "------------------------------------" << std::endl;
	   // std::cout << "------------CLUIT X:: " << cluit->sizeX() << std::endl;
	   // std::cout << "------------CLUIT Y:: " << cluit->sizeY() << std::endl;
	   // std::cout << "------------CLUIT SIZE:: " << cluit->size() << std::endl;
           // std::cout << "------------------------------------" << std::endl;
           // std::cout << "------------------------------------" << std::endl;
	    // determine the position
            MeasurementPoint mpClu(cluit->x(), cluit->y());
            Local3DPoint localPosClu = geomDetUnit->topology().localPosition(mpClu);
            Global3DPoint globalPosClu = geomDetUnit->surface().toGlobal(localPosClu);

            //fill TkLayout histos
            m_trackerLayoutClustersZR->Fill(globalPosClu.z(), globalPosClu.perp());
            m_trackerLayoutClustersYX->Fill(globalPosClu.x(), globalPosClu.y());


	    //std::cout << "----------------------" << std::endl;
	    //std::cout << "----------------------" << std::endl;
	    //std::cout << "--------ColSpan:: " << cluit->colSpan() << " -- sizeY:: " << cluit->sizeY() << std::endl;
	    //std::cout << "--------RowSpan:: " << cluit->rowSpan() << " -- sizeX:: " << cluit->sizeX() << std::endl;
	    //std::cout << "----------------------" << std::endl;
	    //std::cout << "----------------------" << std::endl;

	    //     MERGE SECTION
	    //
	    
	    //std::cout << "XX---  " << cluit << "  --XX" << std::endl;
            edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter = findSimLinkDetSet(rawid);
	    //std::cout << "Side: " << side << " Layer: " << layer << " Ring: " << ring << " rawID: " << rawid << std::endl;
            //mergeClu = this->isMerged(simLinkDSViter, cluit, false);
            std::set<unsigned int> mergeClu = this->getSimTrackId(simLinkDSViter, cluit, false);
	    m_mergeClus->Fill(mergeClu.size());
	    //if(simTrackId.size() > 1 ) 
	    //	std::cout << "-- Merged Cluster\n";
	    //std::cout<< std::endl;
	    
	    //
	    //     MERGE SECTION

	    CluX      = globalPosClu.x();
	    CluY      = globalPosClu.y();
	    CluZ      = globalPosClu.z();
	    CluPhi    = globalPosClu.phi();
	    CluTheta  = globalPosClu.theta();
	    CluCharge = cluit->charge();
	    CluArea   = (cluit->sizeY())*(cluit->sizeX());
	    CluSize   = cluit->size();
	    CluMerge  = mergeClu.size();

	    m_CluX->Fill(CluX);
	    m_CluY->Fill(CluY);
	    m_CluZ->Fill(CluZ);
	    m_CluTheta->Fill(CluTheta);
	    m_CluPhi->Fill(CluPhi);

	    outTreeCluster->Fill();


	    //m_CluX->Fill(globalPosClu.x());
	    //m_CluY->Fill(globalPosClu.y());
	    //m_CluZ->Fill(globalPosClu.z());
	    //m_CluPhi->Fill(globalPosClu.phi());
	    //m_CluTheta->Fill(globalPosClu.theta());

            //std::cout << globalPosClu.x() << " " << globalPosClu.y() << std::endl;
            if (m_docoincidence) {
                unsigned int coincidenceId;
                edmNew::DetSet<SiPixelCluster>::const_iterator coincidenceCluster;
                bool found = this->findCoincidence2x(detId, globalPosClu, coincidenceId, coincidenceCluster);
                if (found) {
                    m_total2xcoincidences++;
                    x2Counter[hist_id][ring_id]++;
                    //I have found a coincidence hit, so now I need to check the sim Track ID of each of the two clusters
                    //so first, let's get the simTrackID of the original cluster-for this I need to loop over all pixels
                    //and check if multiple SimTracks have caused it
                    //if that is not the case, I use the SimLink DSViter and the DSiter (cluit) to get a simTrackID
                    //I need to use the DetId of the coincidence cluster to get another SimLinkDSViter and the coincidence cluster DSiter for the second simTrackID

                    //now get the simlink detset
                    edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter = findSimLinkDetSet(rawid);
                    std::set<unsigned int> simTrackId = this->getSimTrackId(simLinkDSViter, cluit, false);
                    //now get the simlink detset based on the coincidence hit detid
                    simLinkDSViter = findSimLinkDetSet(coincidenceId);
                    std::set<unsigned int> coincidencesimTrackId = this->getSimTrackId(simLinkDSViter, coincidenceCluster, false);
                    std::set<unsigned int> intersection;
                    bool areSame = areSameSimTrackId(simTrackId, coincidencesimTrackId, intersection);

                    if (areSame) {
                        x2Counterreal[hist_id][ring_id]++;
                    } else
                        m_fake2xcoincidences++;

                    m_trackerLayout2xZR->Fill(globalPosClu.z(), globalPosClu.perp());
                    m_trackerLayout2xYX->Fill(globalPosClu.x(), globalPosClu.y());

                    //done with 2 fold coincidences, now 3 fold
                    //only if we have a 2 fold coincidence we can search for a third one in another ring
                    found = false;
                    unsigned int coincidenceId3x;
                    edmNew::DetSet<SiPixelCluster>::const_iterator coincidenceCluster3x;
                    found = this->findCoincidence3x(detId, globalPosClu, coincidenceId3x, coincidenceCluster3x);
                    if (found) {
                        m_total3xcoincidences++;
                        x3Counter[hist_id][ring_id]++;

                        //now get the simlink detset based on the coincidence hit detid
                        edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter3x = findSimLinkDetSet(coincidenceId3x);
                        std::set<unsigned int> coincidencesimTrackId3x = this->getSimTrackId(simLinkDSViter3x, coincidenceCluster3x, false);
                        std::set<unsigned int> intersection3x;
                        //std::cout << "third track ID: ";
                        //for (auto it : coincidencesimTrackId3x)
                        //std::cout << it;
                        //std::cout << std::endl;
                        areSame = areSameSimTrackId(intersection, coincidencesimTrackId3x, intersection3x);
                        //std::cout << "common track ID: ";
                        //for (auto it : intersection3x)
                        //std::cout << it;
                        //std::cout << std::endl;

                        if (areSame)
                            x3Counterreal[hist_id][ring_id]++;
                        else {
                            m_fake3xcoincidences++;
                        }

                        m_trackerLayout3xZR->Fill(globalPosClu.z(), globalPosClu.perp());
                        m_trackerLayout3xYX->Fill(globalPosClu.x(), globalPosClu.y());
                    }
                }
            }
        } //end of cluster loop
    }     //end of module loop
    
    //std::cout << "------WASSA ------------------\n";
    //ok, now I know the number of clusters/hits per ring per disk and should fill the histogram once for this event
    for (unsigned int i = 0; i < 8; i++) {
        //loop the disks
        for (unsigned int j = 0; j < 5; j++) {
            //and the rings
            m_diskHistosCluster[i]->Fill(j + 1, cluCounter[i][j]);
            m_diskHistosHits[i]   ->Fill(j + 1, hitCounter[i][j]);
            if (m_docoincidence) {
                m_diskHistos2x[i]    ->Fill(j + 1, x2Counter[i][j]);
                m_diskHistos2xreal[i]->Fill(j + 1, x2Counterreal[i][j]);
                m_diskHistos3x[i]    ->Fill(j + 1, x3Counter[i][j]);
                m_diskHistos3xreal[i]->Fill(j + 1, x3Counterreal[i][j]);

            }

	    /*
	    std::cout << "---- Area::" << cluAreaCounter[i][j] << "-- Average:: " << cluAreaCounter [i][j]/cluCounter[i][j]  << std::endl;
	    std::cout << "---- SizeX::" <<  cluSizeXCounter[i][j] << "-- Average:: " << cluSizeXCounter [i][j]/cluCounter[i][j]  << std::endl;
	    std::cout << "---- SizeY::" <<  cluSizeYCounter[i][j] << "-- Average:: " << cluSizeXCounter [i][j]/cluCounter[i][j]  << std::endl;
	    */

	    if(cluCounter[i][j] != 0){	
		    /*
		    std::cout << "---- Area::" << cluAreaCounter[i][j] << std::endl;
		    std::cout << "---- SizeX::" <<  cluSizeXCounter[i][j] << std::endl;
		    std::cout << "---- SizeY::" <<  cluSizeYCounter[i][j] << std::endl;
		    */
		    m_ClusterXsize[i]->Fill(j + 1, cluSizeXCounter[i][j]/cluCounter[i][j]);
		    m_ClusterYsize[i]->Fill(j + 1, cluSizeYCounter[i][j]/cluCounter[i][j]);
		    m_ClusterArea[i]->Fill(j + 1, cluAreaCounter [i][j]/cluCounter[i][j]);
	    }
	    //
	    //m_ClusterXsize[i]->Fill(j + 1, cluSizeXCounter[i][j]/cluCounter[i][j] );
	    //m_ClusterYsize[i]->Fill(j + 1, cluSizeYCounter[i][j]/cluCounter[i][j]  );
	    //m_ClusterYsize[i]->Fill(j + 1, cluAreaCounter [i][j]/cluCounter[i][j]  );
	
        }
    }
    m_nevents++;
}

// ------------ method called once each job just after ending the event loop  ------------
void ITclusterAnalyzer::endJob()
{
    outFileCluster->Write();
    outFileCluster->Close();

    std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
    if (m_docoincidence) {
        std::cout << "IT cluster Analyzer found " << m_fake2xcoincidences / (double)m_total2xcoincidences * 100 << "\% fake double coincidences." << std::endl;
        std::cout << "IT cluster Analyzer found " << m_fake3xcoincidences / (double)m_total3xcoincidences * 100 << "\% fake triple coincidences." << std::endl;
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ITclusterAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
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

bool ITclusterAnalyzer::findCoincidence2x(DetId thedetid, Global3DPoint theglobalPosClu, unsigned int& foundDetId, edmNew::DetSet<SiPixelCluster>::const_iterator& foundCluster)
{
    bool found = false;
    uint32_t rawid = thedetid.rawId();
    uint32_t newid = rawid;
    //now I have the raw ID and can mess with the bits
    //the side, layer and ring are the same and I just have to increment or decrement the module number
    unsigned int themodule = (tTopo->pxfModule(thedetid));
    unsigned int thering = (tTopo->pxfBlade(thedetid));

    //in order to avoid duplicates, only look in the next module clockwise
    //depending on the ring, if I am already on the module with the highest module id in the ring, I need to go to the first one
    uint32_t newmodule = themodule + 1;
    if (thering == 1 && themodule == 20)
        newmodule = 1;
    else if (thering == 2 && themodule == 28)
        newmodule = 1;
    else if (thering == 3 && themodule == 36)
        newmodule = 1;
    else if (thering == 4 && themodule == 44)
        newmodule = 1;
    else if (thering == 5 && themodule == 48)
        newmodule = 1;

    //now encode
    newid = (newid & 0xFFFFFC03) | ((newmodule & 0xFF) << 2);

    //now I have a raw id of the module I want to use
    DetId id(newid);
    //unsigned int ring = (tTopo->pxfBlade(id));
    //unsigned int module = (tTopo->pxfModule(id));

    edmNew::DetSetVector<SiPixelCluster>::const_iterator theit = clusters->find(id);
    if (theit == clusters->end()) {
        return false;
    }

    // Get the geomdet
    const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(id));

    unsigned int nClu = 0;
    //at the end of the day, need to find the closest coincidence hit, so store the minimum 2D distance in a temporary variable and a vector for all values
    double r_min = 1000.;
    std::vector<Residual> r_vec;

    //make the return value end();
    foundCluster = theit->end();

    for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit = theit->begin(); cluit != theit->end(); cluit++) {

        // determine the position
        MeasurementPoint mpClu(cluit->x(), cluit->y());
        Local3DPoint localPosClu = geomDetUnit->topology().localPosition(mpClu);
        Global3DPoint globalPosClu = geomDetUnit->surface().toGlobal(localPosClu);

        //now check that the global position is within the cuts
        if (fabs(globalPosClu.x() - theglobalPosClu.x()) < m_dx
            && fabs(globalPosClu.y() - theglobalPosClu.y()) < m_dy
            && fabs(globalPosClu.z() - theglobalPosClu.z()) < m_dz) {
            nClu++;

            double delta_x = fabs(globalPosClu.x() - theglobalPosClu.x());
            double delta_y = fabs(globalPosClu.y() - theglobalPosClu.y());
            Residual r(delta_x, delta_y);
            r_vec.push_back(r);

            if (r.dr < r_min) {
                r_min = r.dr;
                found = true;
                // I assign this here to always have the closest cluster
                foundCluster = cluit;
                foundDetId = newid;
            }

            //std::cout << "Found matching cluster # " << nClu << std::endl;

            //std::cout << "Original x: " << theglobalPosClu.x() << " y: " << theglobalPosClu.y() << " z: " << theglobalPosClu.z() << " ring: " << thering << " module: " << themodule << std::endl;
            //std::cout << "New      x: " << globalPosClu.x() << " y: " << globalPosClu.y() << " z: " << globalPosClu.z() << " ring: " << ring << " module: " << module << std::endl;
        }
    }
    if (nClu > 1) {
        //std::cout << "Warning, found " << nClu << "Clusters within the cuts - the minimum distance is " << r_min << "!" << std::endl;
        //std::cout << "All distances: ";
        for (auto r : r_vec) {
            //r.print();
            if (r.dr == r_min) {
                m_residualX->Fill(r.dx);
                m_residualY->Fill(r.dy);
                m_residualR->Fill(r.dr);
            }
        }
        //std::cout << std::endl;
    }
    return found;
}

bool ITclusterAnalyzer::findCoincidence3x(DetId thedetid, Global3DPoint theglobalPosClu, unsigned int& foundDetId, edmNew::DetSet<SiPixelCluster>::const_iterator& foundCluster)
{
    bool found = false;
    uint32_t rawid = thedetid.rawId();
    uint32_t newid = rawid;
    //now I have the raw ID and can mess with the bits
    //the side and layer are the same and I just have to look in a lower ring
    //unsigned int themodule = (tTopo->pxfModule(thedetid));
    unsigned int thering = (tTopo->pxfBlade(thedetid));

    unsigned int maxmodule = 0;
    unsigned int newring = 0;

    if (thering > 1 && thering <= 5) {
        newring = thering - 1; //look in a lower ring
        //now get the max range of modules for the respective ring
        if (newring == 1)
            maxmodule = 20;
        else if (newring == 2)
            maxmodule = 28;
        else if (newring == 3)
            maxmodule = 36;
        else if (newring == 4)
            maxmodule = 44;
    } else
        return false;

    unsigned int nClu = 0;
    //to make sure we only use the closest hit
    double r_min = 1000.;
    std::vector<Residual> r_vec;

    //make the return value end();
    //foundCluster = theit->end();

    for (uint32_t newmodule = 1; newmodule <= maxmodule; newmodule++) {
        newid = (newid & 0xFFFC0C03) | ((newring & 0x3F) << 12) | ((newmodule & 0xFF) << 2);

        DetId id(newid);
        //unsigned int ring = (tTopo->pxfBlade(id));
        //unsigned int module = (tTopo->pxfModule(id));

        edmNew::DetSetVector<SiPixelCluster>::const_iterator theit = clusters->find(id);
        if (theit == clusters->end()) {
            return false;
        }
        // Get the geomdet
        const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(id));
        if (!geomDetUnit)
            continue;

        for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit = theit->begin(); cluit != theit->end(); cluit++) {

            // determine the position
            MeasurementPoint mpClu(cluit->x(), cluit->y());
            Local3DPoint localPosClu = geomDetUnit->topology().localPosition(mpClu);
            Global3DPoint globalPosClu = geomDetUnit->surface().toGlobal(localPosClu);

            //now check that the global position is within the cuts
            if (fabs(globalPosClu.x() - theglobalPosClu.x()) < m_dx
                && fabs(globalPosClu.y() - theglobalPosClu.y()) < m_dy
                && fabs(globalPosClu.z() - theglobalPosClu.z()) < m_dz) {
                nClu++;
                //found = true;
                //foundCluster = cluit;
                //foundDetId = newid;

                double delta_x = fabs(globalPosClu.x() - theglobalPosClu.x());
                double delta_y = fabs(globalPosClu.y() - theglobalPosClu.y());
                Residual r(delta_x, delta_y);
                r_vec.push_back(r);

                if (r.dr < r_min) {
                    r_min = r.dr;
                    found = true;
                    // i assign this here to be sure to always have the closest cluster
                    foundCluster = cluit;
                    foundDetId = newid;
                    //std::cout << "New det ID: " << newid << std::endl;
                }
                //std::cout << "Found matching cluster # " << nClu << " which is a 3x coincidence" << std::endl;
                //std::cout << "Original x: " << theglobalPosClu.x() << " y: " << theglobalPosClu.y() << " z: " << theglobalPosClu.z() << " ring: " << thering << " module: " << themodule << std::endl;
                //std::cout << "New      x: " << globalPosClu.x() << " y: " << globalPosClu.y() << " z: " << globalPosClu.z() << " ring: " << ring << " module: " << module << std::endl;
            }
        }
    }
    //if (found && nClu > 1) {
    //std::cout << "3X COINCIDENCE - Warning, found " << nClu << "Clusters within the cuts - the minimum distance is " << r_min << "!" << std::endl;
    //std::cout << "All distances: ";
    //for (auto r : r_vec) {
    //r.print();
    //}
    //std::cout << std::endl;
    //}
    return found;
}

edm::DetSetVector<PixelDigiSimLink>::const_iterator ITclusterAnalyzer::findSimLinkDetSet(unsigned int thedetid)
{
    ////basic template
    edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDS = simlinks->find(thedetid);
    return simLinkDS;
}

bool ITclusterAnalyzer::isMerged(edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter, edmNew::DetSet<SiPixelCluster>::const_iterator cluster, bool print)
{
    int size = cluster->size();
    //std::vector<unsigned int> simTrackIds;
    std::set<unsigned int> simTrackIds;

    for (int i = 0; i < size; i++) {

        SiPixelCluster::Pixel pix = cluster->pixel(i);
        unsigned int clusterChannel = PixelDigi::pixelToChannel(pix.x, pix.y);

        if (simLinkDSViter != simlinks->end()) {
            for (edm::DetSet<PixelDigiSimLink>::const_iterator it = simLinkDSViter->data.begin(); it != simLinkDSViter->data.end(); it++) {
                if (clusterChannel == it->channel()) {
                    //simTrackIds.push_back(it->SimTrackId());
                    simTrackIds.insert(it->SimTrackId());
                    if (print)
                        std::cout << "---- "
				  << " Pos (x,y): "        << pix.x << " , " << pix.y 
				  << " Cluster size: "     << size 
				  << " Cluster Channel: "  << clusterChannel 
				  << " SimTrack Channel: " << it->channel()
				  << " SimTrack ID: "      << it->SimTrackId() 
				  << std::endl;
                }
            }
        }
    }

    //std::cout << "---------------------\n";
    //for(std::set<unsigned int>::iterator i = simTrackIds.begin() ;i != simTrackIds.end() ; i++)
    //	std::cout << " Tracks " << *i << std::endl;
    if(simTrackIds.size() != 1){
      //std::cout << "WARNING: have more than 1 simTrackId for this cluster! " << std::endl;
      return true;
    }
    //return simTrackIds;
    return false;
}

std::set<unsigned int> ITclusterAnalyzer::getSimTrackId(edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter, edmNew::DetSet<SiPixelCluster>::const_iterator cluster, bool print)
{
    int size = cluster->size();
    std::set<unsigned int> simTrackIds;

    for (int i = 0; i < size; i++) {

        SiPixelCluster::Pixel pix = cluster->pixel(i);
        unsigned int clusterChannel = PixelDigi::pixelToChannel(pix.x, pix.y);

        if (simLinkDSViter != simlinks->end()) {
            for (edm::DetSet<PixelDigiSimLink>::const_iterator it = simLinkDSViter->data.begin(); it != simLinkDSViter->data.end(); it++) {
                if (clusterChannel == it->channel()) {
                    simTrackIds.insert(it->SimTrackId());
                    if (print)
                        std::cout << "Channel: " << clusterChannel << " SimTrack ID: " << it->SimTrackId() << std::endl;
                }
            }
        }
    }
    //if(simTrackIds.size() != 1){
    //std::cout << "WARNING: have more than 1 simTrackId for this cluster! " << std::endl;
    //return 0;
    //}
    return simTrackIds;
}

bool ITclusterAnalyzer::areSameSimTrackId(std::set<unsigned int> first, std::set<unsigned int> second, std::set<unsigned int>& intersection)
{
    //method to check if the sim Track id is present in both sets
    //std::set<unsigned int> intersection;
    std::set_intersection(first.begin(), first.end(), second.begin(), second.end(), std::inserter(intersection, intersection.begin()));
    if (!intersection.size()) {
        //std::cout << "WARNING, these clusters have not been caused by the same SimTrackID" << std::endl;
        return false;
    } else if (intersection.size() == 1) {
        return true;
    } else {
        //std::cout << "WARNING: both clusters caused by multiple tracks!" << std::endl;
        return true;
    }
}

DEFINE_FWK_MODULE(ITclusterAnalyzer);
