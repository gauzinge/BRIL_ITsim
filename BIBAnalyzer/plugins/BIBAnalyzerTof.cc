// -*- C++ -*-
//
// Package:    BRIL_ITsim/BIBAnalyzerTof
// Class:      BIBAnalyzerTof
//
/**\class BIBAnalyzerTof BIBAnalyzerTof.cc BRIL_ITsim/BIBAnalyzerTof/plugins/BIBAnalyzerTof.cc

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
#include "FWCore/Framework/interface/ConsumesCollector.h"

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

class BIBAnalyzerTof : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit BIBAnalyzerTof(const edm::ParameterSet&);
    ~BIBAnalyzerTof();

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

    // recohits and pshimits matching
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit>> m_tokenRechHits;
    //TrackerHitAssociator::Config trackerHitAssociatorConfig_;

    const double occupancy = 1.722332451499118e-04;//100/(672*864) [%] 2x2 or 1x4 module size from RD53B manual
 
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
 
    // array for particle type names
    const char* ptypes[8] = {"photons", "e+-", "muons+-", "neutrons", "protons", "kaons+-", "pions+-", "residuals"}; 
   
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

    //1D histos for occupancy per ptype per disk per ring
    TH1F* m_histoOccAll[8][5];
    TH1F* m_histoOccPhot[8][5];
    TH1F* m_histoOccElPM[8][5];
    TH1F* m_histoOccMuPM[8][5];
    TH1F* m_histoOccNeut[8][5];
    TH1F* m_histoOccProt[8][5];
    TH1F* m_histoOccKaPM[8][5];
    TH1F* m_histoOccPiPM[8][5];
    TH1F* m_histoOccRest[8][5];

    //event counter
    uint32_t m_nevents;
 
    //hdf5 setup 
    // group in file to group datasets
    Group* group;    
    const char* group_name = "/Tof_q/";
 
    //hdf5 dataset name containing tof and q per particle type
    const H5std_string DATASET_NAME[8];
    
    //initialize memory data space as 1x2 array and unlimited max. rows for each particle type
    const int RANK = 2;
    const hsize_t dims[2] = {1, 2};
    const hsize_t maxdims[2] = {H5S_UNLIMITED, 2};
    
    // memory space for hd5 per ptype per disc per ring
    DataSpace* mspacePhot[8][5];
    DataSpace* mspaceElPM[8][5];
    DataSpace* mspaceMuPM[8][5];
    DataSpace* mspaceNeut[8][5];
    DataSpace* mspaceProt[8][5];
    DataSpace* mspaceKaPM[8][5];
    DataSpace* mspacePiPM[8][5];
    DataSpace* mspaceRest[8][5];

    //enable chunking to have infinite dimension
    DSetCreatPropList* props;
   
    //DSetCreatPropList* props = new DSetCreatPropList;
    const hsize_t chunk_dims[2] = {1, 2};
    
    //Set default values and datatype for the dataset elements
    const float fill_val = 0.0;
 
    //Create hd5 file. If file exists its contents will be overwritten.
    H5File* h5file;
   
    //Create a new dataset in the memory dataspace using properties, for all particle types per disc per ring
    DataSet* datasetPhot[8][5];
    DataSet* datasetElPM[8][5]; 
    DataSet* datasetMuPM[8][5]; 
    DataSet* datasetNeut[8][5]; 
    DataSet* datasetProt[8][5]; 
    DataSet* datasetKaPM[8][5]; 
    DataSet* datasetPiPM[8][5]; 
    DataSet* datasetRest[8][5]; 
   
    //array extension init (always append 1 row), tracked separately for each particle type's dataspace, per disc per ring
    hsize_t sizePhot[8][5][2];
    hsize_t sizeElPM[8][5][2];
    hsize_t sizeMuPM[8][5][2];
    hsize_t sizeNeut[8][5][2];
    hsize_t sizeProt[8][5][2];
    hsize_t sizeKaPM[8][5][2];
    hsize_t sizePiPM[8][5][2];
    hsize_t sizeRest[8][5][2];
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
BIBAnalyzerTof::BIBAnalyzerTof(const edm::ParameterSet& iConfig)
	: m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters")))
        , m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis")))
        , m_tokenSimHits_low(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_low")))
        , m_tokenSimHits_high(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_high")))
        , m_tokenSimLinks(consumes<edm::DetSetVector<PixelDigiSimLink>>(iConfig.getParameter<edm::InputTag>("simlinks")))
	, m_nEBins(iConfig.getUntrackedParameter<uint32_t>("nEBins")) 
        , m_maxEBin(iConfig.getUntrackedParameter<uint32_t>("maxEBin"))
        , m_tokenRechHits(consumes<edmNew::DetSetVector<SiPixelRecHit>>(edm::InputTag("siPixelRecHits")))
//	, trackerHitAssociatorConfig_(iConfig, consumesCollector())
       {
 
    //now do what ever initialization is needed
    m_nevents = 0;

    // init hdf5 file
    h5file = new H5File("/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/bib_simulations/hdf5/tof_q.h5", H5F_ACC_TRUNC);

    // create group hierarchy: one group for everything
    group = new Group( h5file->createGroup(group_name)); 
   
    // hdf5 file properties
    props = new DSetCreatPropList;
    
    //enable chunking to have infinite dimension
    props->setChunk(RANK, chunk_dims);
 
    //init default values and datatype for the dataset elements
    props->setFillValue(PredType::NATIVE_FLOAT, &fill_val);
    
    // declares in heap (8=number of TEPX discs); reserve memory space and create datasets therein
    // loop over discs
    for(int i=0; i<8; i++){

        // disk ID if [-4, 1] for -Z side and [1, 4] for +Z side
        int disk = (i < 4) ? i - 4 : i - 3;
    
        // loop over rings
        for (unsigned int r = 0; r < 5; r++) {

        mspacePhot[i][r] = new DataSpace(RANK, dims, maxdims);
        mspaceElPM[i][r] = new DataSpace(RANK, dims, maxdims);
        mspaceMuPM[i][r] = new DataSpace(RANK, dims, maxdims);
        mspaceNeut[i][r] = new DataSpace(RANK, dims, maxdims);
        mspaceProt[i][r] = new DataSpace(RANK, dims, maxdims);
        mspaceKaPM[i][r] = new DataSpace(RANK, dims, maxdims);
        mspacePiPM[i][r] = new DataSpace(RANK, dims, maxdims);
        mspaceRest[i][r] = new DataSpace(RANK, dims, maxdims);

 	// loop over particle types and cdreate dataset for the current disc and ring
	// specify dataset name for each particle type
        char disk_buffer[100];
	sprintf(disk_buffer, "%d", disk);
	char ring_buffer[100];
	sprintf(ring_buffer, "%d", r+1);

	char dset_name_phot[100];
	strcpy(dset_name_phot, group_name);
        strcat(dset_name_phot, ptypes[0]);
        strcat(dset_name_phot, "d");
        strcat(dset_name_phot, disk_buffer);
        strcat(dset_name_phot, "r");
        strcat(dset_name_phot, ring_buffer);
	
        char dset_name_elpm[100];
        strcpy(dset_name_elpm, group_name);
        strcat(dset_name_elpm, ptypes[1]);
        strcat(dset_name_elpm, "d");
        strcat(dset_name_elpm, disk_buffer);
        strcat(dset_name_elpm, "r");
        strcat(dset_name_elpm, ring_buffer);
	
	char dset_name_mupm[100];
	strcpy(dset_name_mupm, group_name);
        strcat(dset_name_mupm, ptypes[2]);
        strcat(dset_name_mupm, "d");
        strcat(dset_name_mupm, disk_buffer);
        strcat(dset_name_mupm, "r");
        strcat(dset_name_mupm, ring_buffer);

	char dset_name_neut[100];
	strcpy(dset_name_neut, group_name);
        strcat(dset_name_neut, ptypes[3]);
        strcat(dset_name_neut, "d");
        strcat(dset_name_neut, disk_buffer);
        strcat(dset_name_neut, "r");
        strcat(dset_name_neut, ring_buffer);
	
	char dset_name_prot[100];
	strcpy(dset_name_prot, group_name);
        strcat(dset_name_prot, ptypes[4]);
        strcat(dset_name_prot, "d");
        strcat(dset_name_prot, disk_buffer);
        strcat(dset_name_prot, "r");
        strcat(dset_name_prot, ring_buffer);
	
	char dset_name_kapm[100];
	strcpy(dset_name_kapm, group_name);
        strcat(dset_name_kapm, ptypes[5]);
        strcat(dset_name_kapm, "d");
        strcat(dset_name_kapm, disk_buffer);
        strcat(dset_name_kapm, "r");
        strcat(dset_name_kapm, ring_buffer);
	
	char dset_name_pipm[100];
	strcpy(dset_name_pipm, group_name);
        strcat(dset_name_pipm, ptypes[6]);
        strcat(dset_name_pipm, "d");
        strcat(dset_name_pipm, disk_buffer);
        strcat(dset_name_pipm, "r");
        strcat(dset_name_pipm, ring_buffer);
	
	char dset_name_rest[100];
	strcpy(dset_name_rest, group_name);
        strcat(dset_name_rest, ptypes[7]);
        strcat(dset_name_rest, "d");
        strcat(dset_name_rest, disk_buffer);
        strcat(dset_name_rest, "r");
        strcat(dset_name_rest, ring_buffer);
	
 	// create dataset for particel type data
	datasetPhot[i][r] = new DataSet(h5file->createDataSet(dset_name_phot, PredType::NATIVE_FLOAT, *mspacePhot[i][r], *props));
	datasetElPM[i][r] = new DataSet(h5file->createDataSet(dset_name_elpm, PredType::NATIVE_FLOAT, *mspaceElPM[i][r], *props));
	datasetMuPM[i][r] = new DataSet(h5file->createDataSet(dset_name_mupm, PredType::NATIVE_FLOAT, *mspaceMuPM[i][r], *props));
	datasetNeut[i][r] = new DataSet(h5file->createDataSet(dset_name_neut, PredType::NATIVE_FLOAT, *mspaceNeut[i][r], *props));
	datasetProt[i][r] = new DataSet(h5file->createDataSet(dset_name_prot, PredType::NATIVE_FLOAT, *mspaceProt[i][r], *props));
	datasetKaPM[i][r] = new DataSet(h5file->createDataSet(dset_name_kapm, PredType::NATIVE_FLOAT, *mspaceKaPM[i][r], *props));
	datasetPiPM[i][r] = new DataSet(h5file->createDataSet(dset_name_pipm, PredType::NATIVE_FLOAT, *mspacePiPM[i][r], *props));
	datasetRest[i][r] = new DataSet(h5file->createDataSet(dset_name_rest, PredType::NATIVE_FLOAT, *mspaceRest[i][r], *props));

        //array extension init (always append 1 row)
        sizePhot[i][r][0] = 1;
        sizePhot[i][r][1] = 2;
        sizeElPM[i][r][0] = 1;
        sizeElPM[i][r][1] = 2;
        sizeMuPM[i][r][0] = 1;
        sizeMuPM[i][r][1] = 2;
        sizeNeut[i][r][0] = 1;
        sizeNeut[i][r][1] = 2;
        sizeProt[i][r][0] = 1;
        sizeProt[i][r][1] = 2;
        sizeKaPM[i][r][0] = 1;
        sizeKaPM[i][r][1] = 2;
        sizePiPM[i][r][0] = 1;
        sizePiPM[i][r][1] = 2;
        sizeRest[i][r][0] = 1;
        sizeRest[i][r][1] = 2;
        } // end of ring loop
    } // end of disc loop 
} // end of constructor


BIBAnalyzerTof::~BIBAnalyzerTof() {
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
    // close created datasets and data spaces
    for (int i=0; i<8; i++){
	for (int r=0; r<5; r++){
        delete datasetPhot[i][r];
  	delete mspacePhot[i][r];
        delete datasetElPM[i][r];
  	delete mspaceElPM[i][r];
        delete datasetMuPM[i][r];
  	delete mspaceMuPM[i][r];
        delete datasetNeut[i][r];
  	delete mspaceNeut[i][r];
        delete datasetProt[i][r];
  	delete mspaceProt[i][r];
        delete datasetKaPM[i][r];
  	delete mspaceKaPM[i][r];
        delete datasetPiPM[i][r];
  	delete mspacePiPM[i][r];
        delete datasetRest[i][r];
  	delete mspaceRest[i][r];
        }
    }
    props->close();
    h5file->close();
}

//
// member functions
//


// ------------ method called once each job just before starting event loop  ------------
void BIBAnalyzerTof::beginJob() {

    edm::Service<TFileService> fs;
    
    /*
    * TEPX histograms for hits
    */    
    fs->file().cd("/");
    TFileDirectory td = fs->mkdir("TEPX_tof");
    td = fs->mkdir("TEPX_tof/Hits_tof");

    fs->file().cd("/");
    TFileDirectory oc = fs->mkdir("TEPX_occupancy");
    oc = fs->mkdir("TEPX_occupancy/Hits_occupancy");

    //loop over disks
    for (unsigned int i = 0; i < 8; i++) {

        // disk ID if [-4, 1] for -Z side and [1, 4] for +Z side
        int disk = (i < 4) ? i - 4 : i - 3;
        
	//TH2F 2D array
        for (unsigned int r = 1; r < 6; r++) {

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
      
               //occupancy plots     
	       m_histoOccAll[i][r-1] = oc.make<TH1F>(std::string("Occupancy for all particles on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for all particles", 48, .5, 48.5); 
               m_histoOccPhot[i][r-1] = oc.make<TH1F>(std::string("Occupancy for photons on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for photons", 48, .5, 48.5); 
               m_histoOccElPM[i][r-1] = oc.make<TH1F>(std::string("Occupancy for e+/- on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for e+/-", 48, .5, 48.5); 
               m_histoOccMuPM[i][r-1] = oc.make<TH1F>(std::string("Occupancy for muons+/- on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for muons+/-", 48, .5, 48.5); 
               m_histoOccNeut[i][r-1] = oc.make<TH1F>(std::string("Occupancy for neutrons on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for neutrons", 48, .5, 48.5); 
               m_histoOccProt[i][r-1] = oc.make<TH1F>(std::string("Occupancy for protons on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for protons", 48, .5, 48.5); 
               m_histoOccKaPM[i][r-1] = oc.make<TH1F>(std::string("Occupancy for kaons+/- on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for kaons+/-", 48, .5, 48.5); 
               m_histoOccPiPM[i][r-1] = oc.make<TH1F>(std::string("Occupancy for pions+/- on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for pions+/-", 48, .5, 48.5); 
               m_histoOccRest[i][r-1] = oc.make<TH1F>(std::string("Occupancy for residuals on disk ") + disk + " and Ring " + r , "Occupancy of pixel digis for residuals", 48, .5, 48.5); 
        }
    }

    // global coord. vs tof
    m_histoTofZ = td.make<TH2F>("ToF vs Z", "Time of flight vs. global Z coordinate", 6000, -300, 300, 6500, -50, 600);
}


// ------------ method called for each event  ------------
void BIBAnalyzerTof::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    //get the clusters
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> tclusters;
    iEvent.getByToken(m_tokenClusters, tclusters);

    //get rechitcollection
    edm::Handle<SiPixelRecHitCollection> recHitColl;
    iEvent.getByToken(m_tokenRechHits, recHitColl);
    //TrackerHitAssociator associate(iEvent, trackerHitAssociatorConfig_);

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
        
		//get particle type
//                int ptype = DSVit->particleType();
       
        	// get pixel charge in adc counts
                int digitCharge = digit->adc();
		//std::cout << digitCharge << std::endl;
		
		// add to occupancy
		if (digitCharge > 0){
	   	
/*			// fill up particle type dependent histos
	    		switch (std::abs(ptype)){
				case 22: // photon
					m_histoOccPhot[disk_id][ring_id]->AddBinContent(module_id, occupancy);	    			
		    			break;
                 		case 11: // electron/positron   
                    			m_histoOccElPM[disk_id][ring_id]->AddBinContent(module_id, occupancy);
		    			break; 
		 		case 13: // muon+/-
                    			m_histoOccMuPM[disk_id][ring_id]->AddBinContent(module_id, occupancy);
		    			break; 
				case 2112: // neutron/antineutron 
                    			m_histoOccNeut[disk_id][ring_id]->AddBinContent(module_id, occupancy);
		    			break; 
				case 2212: // proton/antiproton
	    	    			m_histoOccProt[disk_id][ring_id]->AddBinContent(module_id, occupancy);
		    			break;
                		case 321: // kaon+/-
	            			m_histoOccKaPM[disk_id][ring_id]->AddBinContent(module_id, occupancy);
                    			break;
                		case 211: // pion+/-
                    			m_histoOccPiPM[disk_id][ring_id]->AddBinContent(module_id, occupancy);
		    			break; 
	        		default: // none of the previous
                    			m_histoOccRest[disk_id][ring_id]->AddBinContent(module_id, occupancy);
         		                break;	
			}
*/			m_histoOccAll[disk_id][ring_id]->AddBinContent(module_id, occupancy);	
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

 	    // get energy deposit in electrons and convert to adc charge
	    float signal_in_eloss = DSVit->energyLoss(); // o(1e-3)
	    //float signal_threshold = 1000.0;
	   
            //parameters for hd5 data writing
            DataSpace fspace;
            hsize_t new_data_offset[2];
            hsize_t new_data_dims[2] = {1, 2};
	    float new_data[1][2] = { {DSVit->timeOfFlight(), signal_in_eloss} };
	
	    // fill up particle type dependent histos
	    switch (std::abs(ptype)){
		case 22: // photon
		    m_histoTofPhot[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetPhot[disk_id][ring_id]->extend(sizePhot[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetPhot[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizePhot[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetPhot[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspacePhot[disk_id][ring_id], fspace);
        	    sizePhot[disk_id][ring_id][0]++;
	
		    break;
                 case 11: // electron/positron   
              	    m_histoTofElPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetElPM[disk_id][ring_id]->extend(sizeElPM[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetElPM[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeElPM[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetElPM[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceElPM[disk_id][ring_id], fspace);
        	    sizeElPM[disk_id][ring_id][0]++;
	
		    break; 
		 case 13: // muon+/-
                    m_histoTofMuPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetMuPM[disk_id][ring_id]->extend(sizeMuPM[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetMuPM[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeMuPM[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetMuPM[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceMuPM[disk_id][ring_id], fspace);
        	    sizeMuPM[disk_id][ring_id][0]++;
	
		    break; 
		case 2112: // neutron/antineutron 
                    m_histoTofNeut[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetNeut[disk_id][ring_id]->extend(sizeNeut[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetNeut[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeNeut[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetNeut[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceNeut[disk_id][ring_id], fspace);
        	    sizeNeut[disk_id][ring_id][0]++;
	
		    break; 
		case 2212: // proton/antiproton
	    	    m_histoTofProt[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetProt[disk_id][ring_id]->extend(sizeProt[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetProt[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeProt[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetProt[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceProt[disk_id][ring_id], fspace);
        	    sizeProt[disk_id][ring_id][0]++;
	
		    break;
                case 321: // kaon+/-
	            m_histoTofKaPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetKaPM[disk_id][ring_id]->extend(sizeKaPM[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetKaPM[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeKaPM[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetKaPM[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceKaPM[disk_id][ring_id], fspace);
        	    sizeKaPM[disk_id][ring_id][0]++;
	
                    break;
                case 211: // pion+/-
                    m_histoTofPiPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetPiPM[disk_id][ring_id]->extend(sizePiPM[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetPiPM[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizePiPM[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetPiPM[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspacePiPM[disk_id][ring_id], fspace);
        	    sizePiPM[disk_id][ring_id][0]++;
	
		    break; 
	        default: // none of the previous
                    m_histoTofRest[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetRest[disk_id][ring_id]->extend(sizeRest[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetRest[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeRest[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetRest[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceRest[disk_id][ring_id], fspace);
        	    sizeRest[disk_id][ring_id][0]++;
	
                    break;	
	    }
		
	    // fill up all particle histos	
            m_histoTofAll[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
	    m_histoTofZ->Fill(global_coords.z(), DSVit->timeOfFlight());
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

	    // get energy deposit in electrons and convert to adc charge
	    float signal_in_eloss = DSVit->energyLoss(); // o(1e-3)
	    //float signal_threshold = 1000.0;

            //parameters for hd5 data writing
            DataSpace fspace;
            hsize_t new_data_offset[2];
            hsize_t new_data_dims[2] = {1, 2};
	    float new_data[1][2] = { {DSVit->timeOfFlight(), signal_in_eloss} };
	
	    // fill up particle type dependent histos
	    switch (std::abs(ptype)){
		case 22: // photon
		    m_histoTofPhot[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetPhot[disk_id][ring_id]->extend(sizePhot[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetPhot[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizePhot[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetPhot[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspacePhot[disk_id][ring_id], fspace);
        	    sizePhot[disk_id][ring_id][0]++;
	
		    break;
                 case 11: // electron/positron   
              	    m_histoTofElPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetElPM[disk_id][ring_id]->extend(sizeElPM[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetElPM[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeElPM[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetElPM[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceElPM[disk_id][ring_id], fspace);
        	    sizeElPM[disk_id][ring_id][0]++;
	
		    break; 
		 case 13: // muon+/-
                    m_histoTofMuPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetMuPM[disk_id][ring_id]->extend(sizeMuPM[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetMuPM[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeMuPM[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetMuPM[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceMuPM[disk_id][ring_id], fspace);
        	    sizeMuPM[disk_id][ring_id][0]++;
	
		    break; 
		case 2112: // neutron/antineutron 
                    m_histoTofNeut[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetNeut[disk_id][ring_id]->extend(sizeNeut[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetNeut[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeNeut[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetNeut[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceNeut[disk_id][ring_id], fspace);
        	    sizeNeut[disk_id][ring_id][0]++;
	
		    break; 
		case 2212: // proton/antiproton
	    	    m_histoTofProt[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetProt[disk_id][ring_id]->extend(sizeProt[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetProt[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeProt[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetProt[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceProt[disk_id][ring_id], fspace);
        	    sizeProt[disk_id][ring_id][0]++;
	
		    break;
                case 321: // kaon+/-
	            m_histoTofKaPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetKaPM[disk_id][ring_id]->extend(sizeKaPM[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetKaPM[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeKaPM[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetKaPM[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceKaPM[disk_id][ring_id], fspace);
        	    sizeKaPM[disk_id][ring_id][0]++;
	
                    break;
                case 211: // pion+/-
                    m_histoTofPiPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetPiPM[disk_id][ring_id]->extend(sizePiPM[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetPiPM[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizePiPM[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetPiPM[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspacePiPM[disk_id][ring_id], fspace);
        	    sizePiPM[disk_id][ring_id][0]++;
	
		    break; 
	        default: // none of the previous
                    m_histoTofRest[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;

        	    //hdf5
		    // make a new row in dataset for new data
		    datasetRest[disk_id][ring_id]->extend(sizeRest[disk_id][ring_id]);

        	    //Select a hyperslab (where to write and what size of data to write)
		    fspace = datasetRest[disk_id][ring_id]->getSpace();
	            new_data_offset[0] = sizeRest[disk_id][ring_id][0] - 1; // write data to last row
	            new_data_offset[1] = 0;
	            fspace.selectHyperslab( H5S_SELECT_SET, new_data_dims, new_data_offset);
	            
	            //Write new_data from memory dataspace to file dataspace.
	  	    datasetRest[disk_id][ring_id]->write(new_data, PredType::NATIVE_FLOAT, *mspaceRest[disk_id][ring_id], fspace);
        	    sizeRest[disk_id][ring_id][0]++;
	
                    break;	
	    }
	
	    // fill up all particle histos	
            m_histoTofAll[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), signal_in_eloss) ;
	    m_histoTofZ->Fill(global_coords.z(), DSVit->timeOfFlight());
        }
    }
    m_nevents++;
}


// ------------ method called once each job just after ending the event loop  ------------
void BIBAnalyzerTof::endJob() {
    std::cout << "ptypes dict" << std::endl;
    for (auto const& pair: ptypes_dict) {
        std::cout << "{" << pair.first << ": " << pair.second << "}" << std::endl;
    }
    std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BIBAnalyzerTof::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

DEFINE_FWK_MODULE(BIBAnalyzerTof);
