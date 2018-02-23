#ifndef HNL_H
#define HNL_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "HNLAnalyzer/PatAnalyzer/interface/GenParticleManager.h"
#include "HNLAnalyzer/PatAnalyzer/interface/Statistics.h"
#include "HNLAnalyzer/PatAnalyzer/interface/Tools.h"
#include "HNLAnalyzer/PatAnalyzer/interface/OnTheFlyCorrections.hh"


//#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"


//Root Classes

#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include "TClonesArray.h"

//Standard C++ classes
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <ostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <memory>
#include <iomanip>

using namespace std;

const int nLeptonsMax = 10; //zhud: originally nLeptonsMax=10;
const int nJetsMax = 30;

const char* _triggersTauNames8[6] =
{"HLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1_v",
 "HLT_Ele24_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1_v",
 "HLT_IsoTkMu22_v",
  "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v",
    "HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v",
    "HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"
};

const char* _triggers1lNames8[8] =
//{"HLT_IsoMu20_v","HLT_IsoTkMu20_v","HLT_IsoMu18_v",
// "HLT_Ele23_WPLoose_Gsf_v","HLT_Ele27_WP85_Gsf_v",
// "HLT_IsoMu22_v","HLT_IsoTkMu22_v"
//};
{"HLT_IsoMu20_v","HLT_IsoTkMu20_v",
"HLT_IsoMu22_v","HLT_IsoMu24_v",
"HLT_Ele27_WPLoose_Gsf_v","HLT_Ele27_WPTight_Gsf_v",
"HLT_Ele27_eta2p1_WPLoose_Gsf_v","HLT_Mu17_Photon22_CaloIdL_L1ISO_v"};

const char* _triggers3lNames8[4] =
{"HLT_TripleMu_12_10_5", "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
"HLT_Mu8_DiEle12_CaloIdL_TrackIdL", "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"
};

const char* _triggers2lNames8[2][6] =
{ {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
   "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
   "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v", "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"},

{"HLT_DoubleMu8_Mass8_PFHT300","empty",
 "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300","empty",
 "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300","empty"}
};


const char* _triggers2lNames8Bkp[2][6] =
{ {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"},
    
    {"HLT_DoubleMu8_Mass8_PFHT250","empty",
        "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250","empty",
        "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250","empty"}
};

const char* _triggersCSNames8[5][5] =
{
    {"HLT_Mu8_TrkIsoVVL_v*","HLT_Mu17_TrkIsoVVL_v*","HLT_Mu24_TrkIsoVVL_v*","HLT_Mu34_TrkIsoVVL_v*","empty"},
    {"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*","HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v*","HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*","HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v*"},
    {"HLT_Mu8_v*","HLT_Mu17_v*","HLT_Mu24_v*","HLT_Mu34_v*","empty"},
    {"HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v*","HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v*","HLT_Ele18_CaloIdM_TrackIdM_PFJet30_v*","HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v*","HLT_Ele33_CaloIdM_TrackIdM_PFJet30_v*"},
    {"HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*","HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*","empty","empty","empty"}
};

const char* _triggersCSbNames8[2] = {"HLT_Mu10_CentralPFJet30_BTagCSV0p5PF_v*","HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5PF_v*"};

const char* _postICHEP8[20] =
{
    "HLT_IsoMu24",
    "HLT_IsoTkMu24",
    
    "HLT_Ele27_WPTight_Gsf",
    
    "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
    "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
    
    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
    "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
    
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
    
    "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300",
    "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300",
    "HLT_DoubleMu8_Mass8_PFHT300",
    
    "HLT_TripleMu_12_10_5",
    "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
    "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
    "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"
};



class HNL : public edm::EDAnalyzer {
public:
 
    explicit HNL(const edm::ParameterSet & iConfig);
    ~HNL(){};
    
private:
    
    //virtual void analyze(edm::Event & iEvent, const edm::EventSetup & iSetup);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob(void);
    
    void fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu);
    void fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* el);
    
    void fillMCVars(const GenParticle* mc, const int leptonCounter);
    void fillCloseJetVars(const int leptonCounter);
    void matchCloseJet(const int leptonCounter);
    void matchAnyJet(const int index);
    double matchCloseParticle(const int leptonCounter, const int pdgIDm);

    void fillIsoMCVars(const int leptonCounter);
    
    double lepMVAvalue(const int leptonCounter);
    
    int photonOrigin(const GenParticle* mc);
    
    void bookTree();
    
    std::vector<const pat::Jet* > SelectedJetsAll;
    edm::Handle<GenParticleCollection> TheGenParticles;

    Vertex::Point PVmc;
    
    std::string Sample;
    edm::InputTag IT_muon;
    edm::InputTag IT_displacedStandAloneMuons;
    edm::InputTag IT_electron;
    edm::InputTag IT_tau;
    edm::InputTag IT_tauDiscriminator;
    edm::InputTag IT_htt;
    edm::InputTag IT_jet;
    edm::InputTag IT_pfmet;
    edm::InputTag IT_beamspot;
    edm::InputTag IT_hltresults;
    edm::InputTag IT_METFilters;
    edm::InputTag IT_genParts;
    
    
    edm::Service<TFileService> fs;
    FILE *outfile;
    
    TH1F *Nvtx;
    
    //desired output variables
    TTree* outputTree;
    
    string _corrLevel;
    
    bool firstEvent_;
  
    double _miniIsoCut;
    double _relIsoCut;
    double _absIsoCut;
    bool _chargeConsistency;
    
    double _minPt0;
    double _minPt1;
    double _tightD0Mu;
    double _tightD0E;
    double _looseD0Mu;
    double _looseD0E;
    double _looseSIP;
    
    double _jetPtCut;
    double _jetEtaCut;
    
    double _tauPt;
    double _tauEta;
    
    bool _regression;
    
    // MVA values and categories (optional)
    edm::EDGetTokenT<bool>                                  ifilterbadChCandToken;
    edm::EDGetTokenT<bool>                                  ifilterbadPFMuonToken;
    
    
    edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
    edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;
    edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesHZZMapToken_;
    edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesHZZMapToken_;
    edm::EDGetTokenT<reco::BeamSpot>                  beamspotToken_;
    edm::EDGetTokenT<std::vector<Vertex> >                  vtxToken_           ;
    edm::EDGetTokenT<double>                                  rhoToken_           ;
    edm::EDGetTokenT<double>                                  fixedGridRhoToken_  ;
    edm::EDGetTokenT< std::vector<PileupSummaryInfo> >        puinfoToken_        ;
    edm::EDGetTokenT< GenEventInfoProduct >                   geneventToken_      ;
    edm::EDGetTokenT< LHEEventProduct >                   lheeventToken_      ;
    
    
    edm::EDGetTokenT<reco::GenParticleCollection>             genparticleToken_   ;
    
    edm::EDGetTokenT<pat::PackedCandidateCollection>          pfcToken_           ;
    edm::EDGetTokenT<pat::JetCollection>                      jetToken_           ;
    edm::EDGetTokenT<pat::MuonCollection>     		            muonToken_  	;
    edm::EDGetTokenT<std::vector<reco::Track>>                       displacedStandAloneMuonsToken_;
    edm::EDGetTokenT<pat::TauCollection>     		            tauToken_  	;
    edm::EDGetTokenT<pat::ElectronCollection>		            electronToken1_	;
    edm::EDGetTokenT<edm::View<pat::Electron> >		            electronToken_	;
    edm::EDGetTokenT< std::vector<reco::Conversion> >     		            convToken_  	;
    edm::EDGetTokenT<pat::METCollection> 	    		            metToken_		;
    edm::EDGetTokenT<pat::JetCollection>                      jetForMetCorrToken_ ;
    
    edm::EDGetTokenT<edm::TriggerResults>                     triggerToken_       ;
    edm::EDGetTokenT<edm::TriggerResults>                     filterToken_       ;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>  triggerObjects_     ;
    //edm::EDGetTokenT<pat::PackedTriggerPrescales>             triggerPrescales_   ;
    
    edm::EDGetTokenT<edm::TriggerResults>                     noiseFilterToken_;
    edm::EDGetTokenT<bool>                                    HBHENoiseFilterLooseResultToken_;
    edm::EDGetTokenT<bool>                                    HBHENoiseFilterTightResultToken_;
    

    
    double _miniisocut[2];
    double _ptratiocut[2];
    double _ptrelcut[2];
    double _multiConst[3][3];

    bool _triggers2l[2][6];
    bool _triggers2lbkp[2][6];
    bool _triggersCS[5][5];
    bool _triggersCSb[2];
    bool _triggers1l[8];
    bool _triggersTau[6];
    bool _triggersPost[20];

    Vertex::Point PV;
    
    
    Float_t LepGood_pt, LepGood_eta, LepGood_jetNDauChargedMVASel,
    LepGood_miniRelIsoCharged, LepGood_miniRelIsoNeutral,
    LepGood_jetPtRelv2, LepGood_jetPtRatio,
    LepGood_jetBTagCSV,
    LepGood_sip3d, LepGood_dxy, LepGood_dz,
    LepGood_segmentCompatibility,
    LepGood_mvaIdSpring15;
    
    TMVA::Reader *reader[2];
    
    double _mvaCut[2];
    
    int _genPhotOrigin;
    
    
    bool Flag_HBHENoiseFilter;
    bool Flag_HBHENoiseIsoFilter;
    bool Flag_globalTightHalo2016Filter;
    bool Flag_CSCTightHalo2015Filter;
    bool Flag_EcalDeadCellTriggerPrimitiveFilter;
    bool Flag_goodVertices;
    bool Flag_eeBadScFilter;
    bool filterbadChCandidate;
    bool filterbadPFMuon;
    
    //std::vector<std::string> myManualCatWeigths;
    //vector<string> myManualCatWeigthsTrig;
    //EGammaMvaEleEstimatorCSA14* myMVATrig;
    double looseMVA[3][2]; //{{0.35, 0.20, -0.52}, {0.73, 0.57, 0.05}};//{0.8, 1.479, };
    

    double myRhoJetsNC;
    double myRhoJets;

    //genlevel particles
    GenParticleManager GPM;
    OnTheFlyCorrections* fMetCorrector;
    
    int _n_bJets;
    int _n_Jets;
    int _n_Jets30;
    
    double _l1HTT;
    
    double mChi20;
    double mChi10;
    
    double _jetEta[nJetsMax];
    double _jetPhi[nJetsMax];
    double _jetPt[nJetsMax];
    double _jetM[nJetsMax];
    double _jetE[nJetsMax];
    int _jethFlav[nJetsMax];
    int _jetpFlav[nJetsMax];
    
    bool _bTagged[nJetsMax];
    double _csv[nJetsMax];
    double _jetDeltaR[nJetsMax][nLeptonsMax];
    double _jetDeltaRloose[nJetsMax];
    
    int _leptonIndex;
    int _closeIndex[nLeptonsMax];
    
    TH1D* flavComp;
    
    
    TClonesArray* _leptonP4;
    TClonesArray* _jetP4;
    TClonesArray* _jetAllP4;
    
    
    //zhud: Displaced Stand Alone Muons
    int _nDSAMu;



    int _nLeptons;
    int _nLooseLeptons;
    int _nEle;
    int _nMu;
    int _nTau;
    double _weight;
    double _genHT;
    
    int _eventType; //ee,mm,em
    bool _sb;
    bool _doubleF;
    int _index1 = -1;
    int _index2 = -1;

    
    int _indeces[nLeptonsMax];
    int _flavors[nLeptonsMax];
    int _charges[nLeptonsMax];
    int _chargesMC[nLeptonsMax];
    double _isolation[nLeptonsMax][2];
    double _miniisolation[nLeptonsMax][3];
    double _miniisolationCharged[nLeptonsMax][3];
    bool _multiisolation[nLeptonsMax][3];//3 WP
    double _isolationMC[nLeptonsMax][4]; //all daughters; all visible daughters; all daughters in the cone; all visible daughters in the cone
    double _ptrel[nLeptonsMax];
    double _ptrel2[nLeptonsMax];
    double _ptratio[nLeptonsMax];
    
    bool _chargeCons[nLeptonsMax];
    
    double _lepDeltaRloose[nLeptonsMax];
    double _lepDeltaR[nLeptonsMax][nLeptonsMax];
    
    int _trackSelectionMultiplicity[nLeptonsMax];
    
    double _dptoverpt[nLeptonsMax];
    double _mvaValue[nLeptonsMax];
    double _muonSegmentComp[nLeptonsMax];
    double _lepMVA[nLeptonsMax];
    
    double _mt[nLeptonsMax];
    double _mllZ[nLeptonsMax];
    double _mllG[nLeptonsMax];
    double _mllZj[nLeptonsMax];
    double _mllGj[nLeptonsMax];
    double _mll[nLeptonsMax][nLeptonsMax];
    
    int _origin[nLeptonsMax];
    int _originPhot[nLeptonsMax];
    int _originDetailed[nLeptonsMax];
    bool _isPromptFinalState[nLeptonsMax];
    bool _fromHardProcessFinalState[nLeptonsMax];
    double _PVchi2;
    double _PVerr[3];
    double _ipPV[nLeptonsMax];
    double _ipPVerr[nLeptonsMax];
    double _ipPVmc[nLeptonsMax];
    
    double _ipZPV[nLeptonsMax];
    double _ipZPVerr[nLeptonsMax];
    
    double _3dIP[nLeptonsMax];
    double _3dIPerr[nLeptonsMax];
    double _3dIPsig[nLeptonsMax];
    
    int _missingHits[nLeptonsMax];
    
    
    double _closeJetPtAll[nLeptonsMax];
    double _closeJetEtaAll[nLeptonsMax];
    double _closeJetPhiAll[nLeptonsMax];
    double _closeJetMAll[nLeptonsMax];
    double _closeJetEAll[nLeptonsMax];
    double _closeJetCSVAll[nLeptonsMax];
    int _closeJetNconstAll[nLeptonsMax];
    double _closeJetAngAll[nLeptonsMax];
    double _ptRelAll[nLeptonsMax];
    
    double _closeJetPtAllMC[nLeptonsMax];
    double _closeJetPtAllstatus[nLeptonsMax];
    int _partonIdMatched[nLeptonsMax];
    bool _sameParton[nLeptonsMax];
    
    bool _isloose[nLeptonsMax];
    bool _islooseID[nLeptonsMax];
    bool _istight[nLeptonsMax];
    bool _istightIso[nLeptonsMax];
    bool _istightID[nLeptonsMax];
    bool _triggerMatch[nLeptonsMax];
    bool _triggerIsoMatch[nLeptonsMax];
    
    bool _tauIDold[nLeptonsMax][5];
    bool _tauIDnew[nLeptonsMax][5];
    bool _tauID03[nLeptonsMax][4];
    
    bool _tauIDeleVeto[nLeptonsMax][5];
    bool _tauIDmuVeto[nLeptonsMax][2];

    bool _tauDecMode[nLeptonsMax][2];
    
    double _tauDeltaR[nLeptonsMax];

    int _n_PV;
    int _n_Interactions;
    double _n_trueInteractions;
    
    double _ptSystem;
    double _ptTTSystem;
    
    int _n_electrons;
    int _n_muons;
    int _n_taus;
    
    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;

    int _nGenLep;
    int _nGenE;
    int _nGenMu;
    int _nGenHNL;
    int _nGenStatusNot1;
    
    //zhud: Begin efficiency KPIs
    int _nGenHNLMu, _nlHNLMu;
    double _GenHNLMuPt[nLeptonsMax], _GenHNLMuEta[nLeptonsMax], _GenHNLMuPhi[nLeptonsMax], _GenHNLMuE[nLeptonsMax];
    double _lHNLMuPtmc[nLeptonsMax];

    double _lPt[nLeptonsMax], _lEta[nLeptonsMax], _lPhi[nLeptonsMax], _lE[nLeptonsMax];
    double _lPtmc[nLeptonsMax], _lEtamc[nLeptonsMax], _lPhimc[nLeptonsMax], _lEmc[nLeptonsMax];
    double _nuPtmc[nLeptonsMax], _nuEtamc[nLeptonsMax], _nuPhimc[nLeptonsMax], _nuEmc[nLeptonsMax];
    double _lMuPt[nLeptonsMax], _lMuEta[nLeptonsMax], _lMuPhi[nLeptonsMax], _lMuE[nLeptonsMax];
    double _GenMuPt[nLeptonsMax], _GenMuEta[nLeptonsMax], _GenMuPhi[nLeptonsMax], _GenMuE[nLeptonsMax];
    double _GenHNLP, _GenHNLEta, _GenHNLPhi, _GenHNLE, _GenHNLMass, _GenHNLGamma, _GenHNLBeta, _GenHNLBetaGamma;

    double _foundGenMuPt[nLeptonsMax];
    bool _isdetectedMu[nLeptonsMax], _isfromHNLMu[nLeptonsMax], _ispromptMu[nLeptonsMax];

    double _GenMuVx[nLeptonsMax], _GenMuVy[nLeptonsMax], _GenMuVxy[nLeptonsMax], _GenMuVz[nLeptonsMax], _GenMuVxyz[nLeptonsMax];
    double _GenHNLVxProd, _GenHNLVyProd, _GenHNLVxyProd, _GenHNLVzProd, _GenHNLVxyzProd;
    double _GenHNLVxDecay[nLeptonsMax], _GenHNLVyDecay[nLeptonsMax], _GenHNLVxyDecay[nLeptonsMax], _GenHNLVzDecay[nLeptonsMax], _GenHNLVxyzDecay[nLeptonsMax], _GenHNLMu3DDisplacement[nLeptonsMax];
    //zhud: End efficiency KPIs
    
    int _pdgmc[nLeptonsMax];

    double _mtmc[nLeptonsMax];
    
    double _mompt[nLeptonsMax];
    double _momphi[nLeptonsMax];
    double _mometa[nLeptonsMax];
    int _mompdg[nLeptonsMax];
    
    
    double _met;
    double _met_phi;
    double _HT;
    
    double _genmet;
    double _genmet_phi;
    
    long _nEventsTotal;
    long _nEventsTotalCounted;
    long _nEventsFiltered;
    
    TH1D* _hCounter;
    
    double _regVars[15];
    double hJet_ptRaw;
    double hJet_genPt;
    double hJet_pt;
    double hJet_phi;
    double hJet_eta;
    double hJet_e;
    
    double hJet_ptLeadTrack;
    
    double hJet_vtx3dL;
    double hJet_vtx3deL;
    double hJet_vtxMass;
    double hJet_vtxPt;
    
    double hJet_cef;
    double hJet_nconstituents;
    double hJet_JECUnc;
    
    double hJet_SoftLeptptRel;
    double hJet_SoftLeptPt;
    double hJet_SoftLeptdR;
    
    double hJet_SoftLeptIdlooseMu;
    double hJet_SoftLeptId95;
};

#endif
