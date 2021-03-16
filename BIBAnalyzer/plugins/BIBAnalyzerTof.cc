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

    float getCharge(Point3DBase<float, LocalTag> local_coords_simhit, DetId simhit_detid, unsigned int simtrackid_simhit, const edm::DetSetVector<PixelDigi>* digis, const edmNew::DetSetVector<SiPixelCluster>* clusters, const edm::DetSetVector<PixelDigiSimLink>* simlinks,  const GeomDetUnit* geomDetUnit);

    float cluster_charge;

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
    const edmNew::DetSetVector<SiPixelCluster>* clusters = NULL;
    const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19
    const std::vector<PSimHit>* simhits_low = NULL;
    const std::vector<PSimHit>* simhits_high = NULL;
    const edm::DetSetVector<PixelDigiSimLink>* simlinks = NULL;

    //max bins of Counting histogram    
    uint32_t m_nQBins;
    uint32_t m_maxQBin;

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
        , m_nQBins(iConfig.getUntrackedParameter<uint32_t>("nQBins")) 
        , m_maxQBin(iConfig.getUntrackedParameter<uint32_t>("maxQBin"))
       {

    std::cout << "BIBAnalyzerTof::BIBAnalyzerTof" << std::endl; 
    //now do what ever initialization is needed
    m_nevents = 0;

    // init hdf5 file
    h5file = new H5File("/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/bib_simulations_tkonly/hdf5/tof_q_test.h5", H5F_ACC_TRUNC);

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
	
 	// create dataset for particle type data
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

    std::cout << "constructor finished" << std::endl;
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
    
    std::cout << "BIBAnalyzerTof::beginJob()" << std::endl;
    /*
    * TEPX histograms for hits
    */    
    fs->file().cd("/");
    TFileDirectory td = fs->mkdir("TEPX_tof");
    td = fs->mkdir("TEPX_tof/Hits_tof");

    //loop over disks
    for (unsigned int i = 0; i < 8; i++) {

        // disk ID if [-4, 1] for -Z side and [1, 4] for +Z side
        int disk = (i < 4) ? i - 4 : i - 3;
        
	//TH2F 2D array
        for (unsigned int r = 1; r < 6; r++) {

	      //ToF - E plots
 	      m_histoTofAll[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for all particles on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for all particles on disk ") + disk + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
              m_histoTofPhot[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for photons on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for photons on disk ") + disk + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
              m_histoTofElPM[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for e+/- on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for e+/- on disk ") + disk + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
              m_histoTofMuPM[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for muons+/- on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for muons+/- on disk ") + disk + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
              m_histoTofNeut[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for neutrons on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for neutrons on disk ") + disk + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
              m_histoTofProt[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for protons on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for protons on disk ") + disk + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
              m_histoTofKaPM[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for kaons+/- on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for kaons+/- on disk ") + disk + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
              m_histoTofPiPM[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for pions+/- on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for pions+/- on disk ") + disk + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
	      m_histoTofRest[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for residuals on disk ") + disk + " and Ring " + r + "", std::string("Charge - time of flight distribution for residuals on disk ") + disk + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
      
        }
    }

    // global coord. vs tof
    m_histoTofZ = td.make<TH2F>("ToF vs Z", "Time of flight vs. global Z coordinate", 6000, -300, 300, 6500, -50, 600);

    std::cout << "beginjob finished" << std::endl;
}


// ------------ method called for each event  ------------
void BIBAnalyzerTof::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    std::cout << "BIBAnalyzerTof::analyze" << std::endl;

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
   
    clusters = tclusters.product(); 
    digis = tdigis.product();
    simhits_low = tsimhits_low.product();
    simhits_high = tsimhits_high.product();
    simlinks = tsimlinks.product();   
 
    std::cout << "number of modules with clusters for event " << m_nevents << ": " << clusters->size() << std::endl;
    std::cout << "number of modules with digis for event " << m_nevents << ": " << digis->size() << std::endl;
    std::cout << "number of simhits low for event " << m_nevents << ": " << simhits_low->size() << std::endl;
    std::cout << "number of simhits high for event " << m_nevents << ": " << simhits_high->size() << std::endl;

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
	std::cout << "Processing digis in detset(=module) side: " << side << ", disk: " << disk  << ", ring: " << ring << ", module: " << module << std::endl;
	//a TEPX module (9-12 starting from origin on both sides)
        if (disk > 8) {
	
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

                //finding the position of the digi
                MeasurementPoint mpHit(digitRow, digitCol);
                Local3DPoint local_coords_digi = geomDetUnit->topology().localPosition(mpHit);
                //Global3DPoint globalPosHit = geomDetUnit->surface().toGlobal(local_coords_digi);
                std::cout << "	Unique local coord of digi: x: " << local_coords_digi.x() << ", y:" << local_coords_digi.y() << ", z: " << local_coords_digi.z() << ", channel:" << digitChannel << ", charge: " << digitCharge << std::endl;

            } // end of pixeldigi loop      
    	} // end of single TEPX module if
    } // end of module loop


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
        unsigned int module = (tTopo->pxfModule(detId)); // values are 1 to 48	
        std::cout << "Simhit found at detset(=module) side: " << side << ", disk: " << disk  << ", ring: " << ring << ", module: " << module << std::endl;
	
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

	    // get local coordinates and simtrackid of PSimHit
	    Point3DBase<float, LocalTag> local_coords_simhit = DSVit->localPosition();
	    unsigned int simtrackid_simhit = DSVit->trackId();
	    std::cout << "with unique local coord: x: " << local_coords_simhit.x() << ", y:" << local_coords_simhit.y() << ", z: " << local_coords_simhit.z() << ", eloss: " << DSVit->energyLoss() << std::endl;
            // match digi or cluster charge to simhit
	    cluster_charge = getCharge(local_coords_simhit, detId, simtrackid_simhit, digis, clusters, simlinks, geomDetUnit);

	    // if there is no activated pixel, continue with next simhit
	    if (cluster_charge==0)
	       continue;	 

            // PixelGeomDetUnit is derived from GeomDet, transform local to global coords
	    GlobalPoint global_coords = geomDetUnit->toGlobal(local_coords_simhit);
      
            //get particle type
            int ptype = DSVit->particleType();

            //parameters for hd5 data writing
            DataSpace fspace;
            hsize_t new_data_offset[2];
            hsize_t new_data_dims[2] = {1, 2};
	    float new_data[1][2] = { {DSVit->timeOfFlight(), cluster_charge} };
	
	    // fill up particle type dependent histos/
	    switch (std::abs(ptype)){
		case 22: // photon
		    m_histoTofPhot[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
              	    m_histoTofElPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
                    m_histoTofMuPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
                    m_histoTofNeut[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
	    	    m_histoTofProt[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
	            m_histoTofKaPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
                    m_histoTofPiPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
                    m_histoTofRest[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
            m_histoTofAll[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;
	    m_histoTofZ->Fill(global_coords.z(), DSVit->timeOfFlight());

	} // end TEPX module if
    } // end psimhit loop

    /*
    * loop over PSimHits high
    */

    for (typename std::vector<PSimHit>::const_iterator DSVit = simhits_high->begin(); DSVit != simhits_high->end(); DSVit++) {

	std::cout << "event " << m_nevents << ", simhit hi" << std::endl;
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

	    // get local coordinates and simtrackid of PSimHit
	    Point3DBase<float, LocalTag> local_coords_simhit = DSVit->localPosition();
	    unsigned int simtrackid_simhit = DSVit->trackId();
	    std::cout << "with unique local coord: x: " << local_coords_simhit.x() << ", y:" << local_coords_simhit.y() << ", z: " << local_coords_simhit.z() << ", eloss: " << DSVit->energyLoss() << std::endl;
            // match digi or cluster charge to simhit
	    cluster_charge = getCharge(local_coords_simhit, detId, simtrackid_simhit, digis, clusters, simlinks, geomDetUnit);

	    // if there is no activated pixel, continue with next simhit
	    if (cluster_charge==0)
	       continue;	 

            // PixelGeomDetUnit is derived from GeomDet, transform local to global coords
	    GlobalPoint global_coords = geomDetUnit->toGlobal(local_coords_simhit);
      
            //get particle type
            int ptype = DSVit->particleType();

            //parameters for hd5 data writing
            DataSpace fspace;
            hsize_t new_data_offset[2];
            hsize_t new_data_dims[2] = {1, 2};
	    float new_data[1][2] = { {DSVit->timeOfFlight(), cluster_charge} };
	
	    // fill up particle type dependent histos
	    switch (std::abs(ptype)){
		case 22: // photon
		    m_histoTofPhot[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
              	    m_histoTofElPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
                    m_histoTofMuPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
                    m_histoTofNeut[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
	    	    m_histoTofProt[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
	            m_histoTofKaPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
                    m_histoTofPiPM[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
                    m_histoTofRest[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;

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
            m_histoTofAll[disk_id][ring_id]->Fill(DSVit->timeOfFlight(), cluster_charge) ;
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

// inputs: digi local coord, digi detid, simhits data, topology
// returns charge of simhit taken as the charge of the closest digi or cluster
float BIBAnalyzerTof::getCharge(Point3DBase<float, LocalTag> local_coords_simhit, DetId simhit_detid, unsigned int simtrackid_simhit, const edm::DetSetVector<PixelDigi>* digis, const edmNew::DetSetVector<SiPixelCluster>* clusters, const edm::DetSetVector<PixelDigiSimLink>* simlinks, const GeomDetUnit* geomDetUnit) {

    // init
    float closest_digi_dist = 999.9;
    int closest_digi_channel = -1; 

    // condition for exiting cluster loop    
    bool charge_found = false;

    // charge of psimhit as a sum of cluster pixel charges with the same simtrackid (matched through simlinks)
    float psimhit_charge = 0;  

    // find digis on the same module as the psimhit
    edm::DetSetVector<PixelDigi>::const_iterator digiDSVit = digis->find(simhit_detid);

    // find clusters on the same module as the psimhit
    edmNew::DetSetVector<SiPixelCluster>::const_iterator clusterDSVit = clusters->find(simhit_detid);
    
    // find simlinks on the same module as the psimhit
    edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSVit = simlinks->find(simhit_detid);	


    // loop over the clusters that are on the module where the simhit is
    if (digiDSVit != digis->end()) {
            //loop over the digis (the pixels of the module with signal) in the matching module
            for (edm::DetSet<PixelDigi>::const_iterator digit = digiDSVit->data.begin(); digit != digiDSVit->end(); digit++) {
        
		int digitRow = digit->row();
		int digitCol = digit->column(); 	
		int digitChannel = PixelDigi::pixelToChannel(digitRow, digitCol);

                //finding the position of the digi
                MeasurementPoint mpHit(digitRow, digitCol);
                Local3DPoint local_coords_digi = geomDetUnit->topology().localPosition(mpHit);
 
		// calculate distance between simhit and digi using basic coordinate geo
		float x_dist = local_coords_simhit.x() - local_coords_digi.x();
		float y_dist = local_coords_simhit.y() - local_coords_digi.y();
		float dist = sqrt(x_dist * x_dist + y_dist * y_dist);
		if (dist < closest_digi_dist) {
		      closest_digi_dist = dist;

		      // get digit channel for cluster matching
		      int digitRow = digit->row();
	              int digitCol = digit->column(); 	
		      closest_digi_channel = PixelDigi::pixelToChannel(digitRow, digitCol);
    	    	} // end of dustance update if
            } // end of digi loop in module
    } // end of digi iterator check

    std::cout << "Closest digi channel: " << closest_digi_channel << std::endl;
    
    // loop over the clusters that are on the module where the simhit is
    if (clusterDSVit != clusters->end()) {
        for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit = clusterDSVit->begin(); cluit != clusterDSVit->end(); cluit++) {
        
	    // determine the barycenter in local reference frame
    	    MeasurementPoint mpClu(cluit->x(), cluit->y());

            // loop over the pixels in the cluster and match pixel channels with digi channel			
            int cluster_size = cluit->size();
            std::cout << "    Cluster size in pixels: " << cluster_size << std::endl;
            for (int i = 0; i < cluster_size; i++) {

	        // one single cluster pixel (digi) and its channel
	        SiPixelCluster::Pixel pix_ = cluit->pixel(i);
                int cluster_pixel_channel_ = PixelDigi::pixelToChannel(pix_.x, pix_.y);

		// if the digi of the psimhit is in this cluster add up charges
		if (cluster_pixel_channel_ == closest_digi_channel) {
		    std::cout << "        This cluster matched with psimhit channel!" << std::endl; 
                    
		    // loop over cluster pixels again from the beginning and add up charges if simtrackid matches
                    for (int j = 0; j < cluster_size; j++) {
			
			SiPixelCluster::Pixel pix = cluit->pixel(j);
                	unsigned int cluster_pixel_channel = PixelDigi::pixelToChannel(pix.x, pix.y);

		        int cluster_pixel_adc = pix.adc; // charge in electrons not adc!
	 	        std::cout << "        Cluster pixel: x: " << pix.x << ", y:" << pix.y << ", channel: " << cluster_pixel_channel << ", charge: " << cluster_pixel_adc << std::endl;

		        //loop over simlinks and match channel
	                if (simLinkDSVit != simlinks->end()) {
		            for (edm::DetSet<PixelDigiSimLink>::const_iterator simlinkit = simLinkDSVit->data.begin(); simlinkit != simLinkDSVit->data.end(); simlinkit++) {
                				
			        // compare cluster pixel channel with simlink channel (now both unique since on same module and correct cluster)
			        if (cluster_pixel_channel == simlinkit->channel()) {
							
				    // get simtrackid of psimhit that activated this pixel (every pixel has a simhit that activated it)
                    		    unsigned int simtrackid_pixel = simlinkit->SimTrackId();
				    std::cout << "            Cluster pixel activated by track: " << simtrackid_pixel << " (psimhit track: " << simtrackid_simhit << ")" << std::endl;	
				
				    // if simtrackid is the same as that of the original psimhit, add charge of pixel to total charge
				    // (this is necessary as clusters can be composed of multiple psimhit activations and we dont want to add charge of other psimhits)
				    if (simtrackid_pixel == simtrackid_simhit){
					psimhit_charge += cluster_pixel_adc;
					std::cout << "                Added charge of " << cluster_pixel_adc << ". Total psimhit charge so far: " << psimhit_charge << std::endl;
				    } // end of simtrackid comparison 
           	 	        } // end of pixel channel comparison
       		             } // end of simlink loop
                        } // end of simlink iterator check
		    } // end of pixel loop for adding charge
		
		    // dont continue with other clusters    
		    charge_found = true;
		    break;
                } // end of matched channel if
	    } // end of pixels in a single cluster loop
        
	if (charge_found) {break;}
        } // end of clusters iterator check
    } // end of cluster on module loop
  
    // case when there is no cluster for psimhit thus no charge has to be handled outside this method 
    return psimhit_charge;
}


DEFINE_FWK_MODULE(BIBAnalyzerTof);
