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
    uint32_t m_maxBin;

    //array of TH2F for psimhits per disk per ring
    TH2F* m_histosPSimHits[8][5]; // 8 discs, 5 rings in each disc

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
        , m_maxBin(iConfig.getUntrackedParameter<uint32_t>("maxBin"))
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
    TFileDirectory td = fs->mkdir("TEPX");
    td = fs->mkdir("TEPX/Hits");

    for (unsigned int i = 0; i < 8; i++) {

        // disc ID if [-4, 1] for -Z side and [1, 4] for +Z side
        int disk = (i < 4) ? i - 4 : i - 3;
        
	//TH1F 2D array
        for (unsigned int r = 0; r < 5; r++) {

               std::stringstream histoname;
               histoname << "Charge - time distribution for Disc " << disk << " and Ring " << r+1 << ";Time of flight;Total induced charge";
               std::stringstream histotitle;
               histotitle << "Charge - time distribution for Disc " << disk << " and Ring " << r+1;

               //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
               m_histosPSimHits[i][r] = td.make<TH2F>(histotitle.str().c_str(), histoname.str().c_str(), 50, -25, 25, 100, 0, 1e-3);
        }
    }

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

    //charge of hits per ring per module
    unsigned int hitCounter[8][5][48];
    memset(hitCounter, 0, sizeof(hitCounter));

    //loop over the entries in the hit collection for TEPX
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
	
	//a TEPX module
        if (disk > 8) {

            //the index in my histogram map
            int hist_id = -1;
            unsigned int ring_id = ring - 1;
            unsigned int module_id = module - 1;
        
	    //this is a TEPX hit on side 1 (-Z)
            if (side == 1) {
                hist_id = disk - 9; // goes from 0 to 3
            }
            
	    //this is a TEPX hit on side 2 (+Z)
            else if (side == 2) {
                hist_id = 4 + disk - 9; // goes from 4 to 7
            }

            //find the geometry of the module associated to the digi
            const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
            if (!geomDetUnit)
                continue;

	    // loop over the simhits in TEPX module
           // for (edm::DetSet<PSimHit>::const_iterator simhit_low = DSVit->begin(); simhit_low != DSVit->end(); simhit_low++) {

            std::cout << "---- Hit info:" << std::endl;
	    std::cout << "hit tof: " << DSVit->timeOfFlight() << std::endl;
            std::cout << "hit eloss: " << DSVit->energyLoss() << std::endl;

            m_histosPSimHits[hist_id][ring_id]->Fill(DSVit->timeOfFlight(), DSVit->energyLoss());
        //    }

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
