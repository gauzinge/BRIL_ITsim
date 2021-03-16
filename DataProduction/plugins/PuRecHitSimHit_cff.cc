// -*- C++ -*-
//
// Package:    BRIL_ITsim/PuRecHitSimHit_cff
// Class:      PuRecHitSimHit_cff
//
/**\class PuRecHitSimHit_cff PuRecHitSimHit_cff.cc
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

// crossing frame to keep simhits
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

#include <TH2F.h>
#include <TTree.h>
#include <TStyle.h>
#include <TNtuple.h>

//write histo data to hd5 file
#include "H5Cpp.h"
using namespace H5;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class PuRecHitSimHit_cff : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit PuRecHitSimHit_cff(const edm::ParameterSet&);
    ~PuRecHitSimHit_cff();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    uint32_t getModuleID(bool, unsigned int, unsigned int, unsigned int);
    // ----------member data ---------------------------
    edm::EDGetTokenT<std::vector<PSimHit>> m_tokenSimHits_low;
    edm::EDGetTokenT<std::vector<PSimHit>> m_tokenSimHits_high;
    edm::EDGetTokenT<CrossingFrame<PSimHit> > m_tokenSimHitsCrossingFrame_low;
    edm::EDGetTokenT<CrossingFrame<PSimHit>> m_tokenSimHitsCrossingFrame_high; 

    edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelRecHit>> m_tokenRecHits;
    TrackerHitAssociator::Config trackerHitAssociatorConfig_;  
 
 
    // the pointers to geometry, topology and clusters
    // these are members so all functions can access them without passing as argument
    const TrackerTopology* tTopo = NULL;
    const TrackerGeometry* tkGeom = NULL;
    const std::vector<PSimHit>* simhits_low = NULL;
    const std::vector<PSimHit>* simhits_high = NULL;
    const CrossingFrame<PSimHit>* simhitscrossingframe_low = NULL;
    const CrossingFrame<PSimHit>* simhitscrossingframe_high = NULL;
    const std::vector<PSimHit>* simhitscrossingfame_low_signal = NULL;

    const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19
    const edmNew::DetSetVector<SiPixelRecHit>* rechits = NULL;
 
    //max bins of Counting histogram    
    uint32_t m_nQBins;
    uint32_t m_maxQBin;

    //jobid for hdf5 file ordering
    uint32_t jobid;

    //histogram for time of flight vs Z coordinate
    TH2F* m_histoTofZ;
 
    //ntuple for storing psimhits
    TNtuple* ntuple_simhits;
    TNtuple* ntuple;

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

    //event counter
    int m_nevents;

    //total rechit counter
    int m_nrechits; 
    int m_nrechits_d4r1; 

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
PuRecHitSimHit_cff::PuRecHitSimHit_cff(const edm::ParameterSet& iConfig)
          : m_tokenSimHits_low(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_low")))
          , m_tokenSimHits_high(consumes<std::vector<PSimHit>>(iConfig.getParameter<edm::InputTag>("simhits_high")))
          , m_tokenSimHitsCrossingFrame_low(consumes<CrossingFrame<PSimHit>>(iConfig.getParameter<edm::InputTag>("crossingframe_simhits_low")))
          , m_tokenSimHitsCrossingFrame_high(consumes<CrossingFrame<PSimHit>>(iConfig.getParameter<edm::InputTag>("crossingframe_simhits_high")))
          , m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis")))
          , m_tokenRecHits(consumes<edmNew::DetSetVector<SiPixelRecHit>>(iConfig.getParameter<edm::InputTag>("rechits")))
	  , trackerHitAssociatorConfig_(iConfig, consumesCollector())
          , m_nQBins(iConfig.getUntrackedParameter<uint32_t>("nQBins")) 
          , m_maxQBin(iConfig.getUntrackedParameter<uint32_t>("maxQBin"))
          , jobid(iConfig.getUntrackedParameter<uint32_t>("jobid"))
       {
   
    std::cout << "PuRecHitSimHit_cff::PuRecHitSimHit_cff" << std::endl; 
    //now do what ever initialization is needed
    m_nevents = 0;
    m_nrechits = 0;
    m_nrechits_d4r1 = 0;

    // init hdf5 file
    std::ostringstream h5fname;
    h5fname << "/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/pu_simulations_fullgeo/hdf5/tof_q_PU_" << jobid << ".h5";
    std::string h5fname_str =  h5fname.str();
    const char* h5fname_chr = h5fname_str.c_str();
    h5file = new H5File(h5fname_chr, H5F_ACC_TRUNC);

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
        int det_disk = (i < 4) ? i - 4 : i - 3;
    
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
	sprintf(disk_buffer, "%d", det_disk);
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


} // end of constructor


PuRecHitSimHit_cff::~PuRecHitSimHit_cff() {
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
void PuRecHitSimHit_cff::beginJob() {

    edm::Service<TFileService> fs;
    
//    std::cout << "PuRecHitSimHit_cff::beginJob()" << std::endl;

    fs->file().cd("/");
    TFileDirectory nt = fs->mkdir("PSimHit");
    
    ntuple_simhits = nt.make<TNtuple>("ntuple","PU simhit data","detid:simtrackid:tof");

    /*
    * TEPX histograms for hits
    */    
    fs->file().cd("/");
    TFileDirectory td = fs->mkdir("TEPX_tof");
    td = fs->mkdir("TEPX_tof/Hits_tof");

    ntuple = td.make<TNtuple>("No. of events processed","No. of events processed","m_nevents");

    //loop over disks
    for (unsigned int i = 0; i < 8; i++) {

        // disk ID if [-4, 1] for -Z side and [1, 4] for +Z side
        int det_disk_hist = (i < 4) ? i - 4 : i - 3;
        
	//TH2F 2D array
        for (unsigned int r = 1; r < 6; r++) {

	      //ToF - E plots
 	      m_histoTofAll[i][r-1] = td.make<TH2F>(std::string("Charge - time of flight distribution for all particles on disk ") + det_disk_hist + " and Ring " + r + "", std::string("Charge - time of flight distribution for all particles on disk ") + det_disk_hist + " and Ring " + r, 550, -50, 500, m_nQBins, 0, m_maxQBin);
        }
    }

    // global coord. vs tof
    m_histoTofZ = td.make<TH2F>("ToF vs Z", "Time of flight vs. global Z coordinate", 6000, -300, 300, 6500, -50, 600);



//    std::cout << "beginjob finished" << std::endl;
}


// ------------ method called for each event  ------------
void PuRecHitSimHit_cff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    //std::cout << "PuRecHitSimHit_cff::analyze" << std::endl;

    //get the simhits - COB 30.11.20
    edm::Handle<std::vector<PSimHit>> tsimhits_low;
    iEvent.getByToken(m_tokenSimHits_low, tsimhits_low);
    edm::Handle<std::vector<PSimHit>> tsimhits_high;
    iEvent.getByToken(m_tokenSimHits_high, tsimhits_high);   

    edm::Handle<CrossingFrame<PSimHit>> tsimhitscrossingframe_low;
    iEvent.getByToken(m_tokenSimHitsCrossingFrame_low, tsimhitscrossingframe_low);
    edm::Handle<CrossingFrame<PSimHit>> tsimhitscrossingframe_high;
    iEvent.getByToken(m_tokenSimHitsCrossingFrame_high, tsimhitscrossingframe_high);   


    //get the digis - COB 26.02.19
    edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
    iEvent.getByToken(m_tokenDigis, tdigis);

    // get rechits
    edm::Handle<SiPixelRecHitCollection> recHitColl;
    iEvent.getByToken(m_tokenRecHits, recHitColl);

    // for finding matched simhit
    TrackerHitAssociator associate(iEvent, trackerHitAssociatorConfig_);

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
   
    simhits_low = tsimhits_low.product();
    simhits_high = tsimhits_high.product();
    simhitscrossingframe_low = tsimhitscrossingframe_low.product();
    simhitscrossingframe_high = tsimhitscrossingframe_high.product();
   
    std::unique_ptr<MixCollection<PSimHit>> simhits(new MixCollection<PSimHit>(tsimhitscrossingframe_low.product()));

    digis = tdigis.product();
    rechits = recHitColl.product();  
 
    // this points the iterators to the beginning and end of the signal vector in the crossing frame 
    const std::vector<const PSimHit*>& simhitscrossingfame_signal = simhitscrossingframe_low->getPileups();

    //std::cout << "number of simhits low for event " << m_nevents+1 << ": " << simhits_low->size() << std::endl;
    //std::cout << "number of simhits high for event " << m_nevents+1 << ": " << simhits_high->size() << std::endl;
    //std::cout << "number of simhits low accessed through crossing frame for event " << m_nevents+1 << ": " << simhitscrossingfame_signal.size() << std::endl;
 
    //std::cout << "number of modules with digis for event " << m_nevents << ": " << digis->size() << std::endl;
    //std::cout << "number of rechit modules for event " << m_nevents << ": " << rechits->size() << std::endl;
    //std::cout << "number of rechits for event " << m_nevents << ": " << rechits->dataSize() << std::endl;


    //for (typename std::vector<PSimHit*>::const_iterator DSVit = simhitscrossingfame_signal.begin(); DSVit != simhitscrossingfame_signal.end(); DSVit++) {
    for (MixCollection<PSimHit>::MixItr ihit = simhits->begin(); ihit != simhits->end(); ihit++) {
 	const PSimHit* simhit = &(*ihit);

	//get the detid
	unsigned int rawid(simhit->detUnitId());
	DetId detId(rawid);

        //find out which layer, side and ring
        unsigned int side = (tTopo->pxfSide(detId));  // values are 1 and 2 for -+Z
        unsigned int disk = (tTopo->pxfDisk(detId)); // values are 1 to 12 for disks TFPX1 to TFPX 8  and TEPX1 to TEPX 4
        unsigned int ring = (tTopo->pxfBlade(detId)); // values are 1 to 5
        unsigned int module = (tTopo->pxfModule(detId)); // values are 1 to 48  

	//std::cout << "side: " << side << ", disk: " << disk << ", ring: " << ring << ", module: " << module << ", signal tof: " << simhit->timeOfFlight() << std::endl;
    
	ntuple_simhits->Fill(detId, simhit->trackId(), simhit->timeOfFlight());
    } 

    std::vector<PSimHit> matched;
    const PSimHit* closest_simhit = NULL;

    // count rechits in current event
    int rechit_counter = 0;

    // Loop over Detector IDs
    for (auto recHitIdIterator : *(rechits)) {
            SiPixelRecHitCollection::DetSet DSVit = recHitIdIterator;
  
            if (DSVit.empty())
                continue;
            DetId detId = DetId(DSVit.detId());  // Get the Detid object

            //find out which layer, side and ring
            unsigned int side = (tTopo->pxfSide(detId));  // values are 1 and 2 for -+Z
            unsigned int disk = (tTopo->pxfDisk(detId)); // values are 1 to 12 for disks TFPX1 to TFPX 8  and TEPX1 to TEPX 4
            unsigned int ring = (tTopo->pxfBlade(detId)); // values are 1 to 5
            unsigned int module = (tTopo->pxfModule(detId)); // values are 1 to 48	
            //std::cout << "Rechits found at detset(=module) side: " << side << ", disk: " << disk  << ", ring: " << ring << ", module: " << module << std::endl;
	
            //a TEPX module (9-12 starting from origin on both sides)
            if (disk > 8 && (side == 1 || side == 2)) {

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
   
                // Loop over rechits for this module
                //for (edmNew::DetSet<SiPixelRecHit>::const_iterator recit = DSVit->begin(); recit != DSVit->end(); recit++){
		for (auto recit : DSVit) { 
                    // clear simhit vector
                    matched.clear();

		    // get vector of simhits corresponding to this rechit
                    matched = associate.associateHit(recit);
		    //std::cout << "size: " << matched.size() << std::endl;
                    if (!matched.empty()) {
                        float closest = 9999.9;
                        LocalPoint local_coords_rechit = recit.localPosition();
                        float rechit_x = local_coords_rechit.x();
                        float rechit_y = local_coords_rechit.y();
   			
			//std::cout << "    Rechit: x: " << rechit_x << ", y: " << rechit_y << std::endl;
 
                        //loop over simhits and find closest
                        for (auto const& m : matched) {
   		 	    Local3DPoint local_coords_simhit = m.localPosition();
		  	    float simhit_x = local_coords_simhit.x();
		   	    float simhit_y = local_coords_simhit.y();
			   // std::cout << "        Simhit: x: " << simhit_x << ", y: " << simhit_y << ", tof: " << m.timeOfFlight() << std::endl;

                            float x_dist = simhit_x - rechit_x;
                            float y_dist = simhit_y - rechit_y;
                            float dist = sqrt(x_dist * x_dist + y_dist * y_dist);
                            if (dist < closest) {
                                closest = dist;
                                closest_simhit = &m;
		   	    }
		        }  // end of simhit loop
	
			// return a cluster pixel and its charge in electrons
			SiPixelRecHit::ClusterRef const& cluster_pixel = recit.cluster();
		        float rechit_charge = cluster_pixel->charge();
			
			// tof in seconds
			float simhit_tof = closest_simhit->timeOfFlight();
		
			//Local3DPoint local_coords_closest_simhit = closest_simhit->localPosition();			

			//std::cout << "    Closest simhit: x: " << local_coords_closest_simhit.x() << ", y: " << local_coords_closest_simhit.y() << ", tof: " << simhit_tof << std::endl;
			//std::cout << "    Rechit charge: " << rechit_charge << std::endl;

            		// PixelGeomDetUnit is derived from GeomDet, transform local to global coords
	    		GlobalPoint global_coords = geomDetUnit->toGlobal(closest_simhit->localPosition());
      
            		//get particle type
            		int ptype = closest_simhit->particleType();
			
                        //parameters for hd5 data writing
            		DataSpace fspace;
            		hsize_t new_data_offset[2];
            		hsize_t new_data_dims[2] = {1, 2};
	    		float new_data[1][2] = { {simhit_tof, rechit_charge} };
	
	    		// fill up all particle histos	
	                m_histoTofAll[disk_id][ring_id]->Fill(simhit_tof, rechit_charge) ;
	    		m_histoTofZ->Fill(global_coords.z(), simhit_tof);

	   	        // fill up particle type dependent histos/
	                switch (std::abs(ptype)){
				case 22: // photon
		    		    //m_histoTofPhot[disk_id][ring_id]->Fill(simhit_tof, rechit_charge) ;

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
		              	    //m_histoTofElPM[disk_id][ring_id]->Fill(simhit_tof, rechit_charge) ;

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
		                    //m_histoTofMuPM[disk_id][ring_id]->Fill(simhit_tof, rechit_charge) ;

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
		                    //m_histoTofNeut[disk_id][ring_id]->Fill(simhit_tof, rechit_charge) ;

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
			    	    //m_histoTofProt[disk_id][ring_id]->Fill(simhit_tof, rechit_charge) ;

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
			            //m_histoTofKaPM[disk_id][ring_id]->Fill(simhit_tof, rechit_charge) ;

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
		                    //m_histoTofPiPM[disk_id][ring_id]->Fill(simhit_tof, rechit_charge) ;

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
		               	    //m_histoTofRest[disk_id][ring_id]->Fill(simhit_tof, rechit_charge) ;

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
			    } // end of switch case
			    m_nrechits++;
			    rechit_counter++;
			    if (side == 2 && disk == 12 && ring == 1){
				m_nrechits_d4r1++;	
			    }
           	   } // end of mathced if 	
               } // end of rechit loop
      }  // end of TEPX module
  }  // end of rechit modules loop
  std::cout << "Number of rechits found in event " << m_nevents+1 << " : " << rechit_counter << std::endl;	
  std::cout << "Number of rechits so far on D4R1: " << m_nrechits_d4r1 << std::endl;
  m_nevents++;
}


// ------------ method called once each job just after ending the event loop  ------------
void PuRecHitSimHit_cff::endJob() {
    std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
    ntuple->Fill(m_nevents);
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PuRecHitSimHit_cff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

DEFINE_FWK_MODULE(PuRecHitSimHit_cff);
