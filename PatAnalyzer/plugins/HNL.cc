#include "HNL.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/CandAlgos/interface/CandMatcher.h"
#include "PhysicsTools/HepMCCandAlgos/interface/MCTruthPairSelector.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include <DataFormats/L1Trigger/interface/L1EtMissParticle.h>



#include "PhysicsTools/CandUtils/interface/CandMatcherNew.h"



using namespace std;
using namespace edm;
using namespace reco;
using namespace tools;
using namespace math;
using namespace reco::tau;

HNL::HNL(const edm::ParameterSet & iConfig) :
    _miniIsoCut(0.4),//loose: abs iso < 30 || reliso < 1
    _relIsoCut(1),
    _absIsoCut(1),
    _chargeConsistency(false),
    _minPt0(1.), //default = 5.
    _minPt1(10.),
    _tightD0Mu(0.01),
    _tightD0E(0.02),
    _looseD0Mu(9999999.),
    _looseD0E(9999999.),
    _looseSIP(9999999.),
    //_looseD0Mu(0.2),
    //_looseD0E(9999999.),
    //_jetPtCut(40.),
    _jetPtCut(25.),
    _jetEtaCut(2.4),
    _tauPt(20),
    _tauEta(2.3),
    ifilterbadChCandToken(consumes<bool>                                (iConfig.getParameter<edm::InputTag>("BadChCandFilter"))),
    ifilterbadPFMuonToken(consumes<bool>                                (iConfig.getParameter<edm::InputTag>("BadPFMuon"))),
    mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
    mvaCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"))),
    mvaValuesHZZMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMapHZZ"))),
    mvaCategoriesHZZMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMapHZZ"))),

    beamspotToken_                    (consumes< reco::BeamSpot >                (iConfig.getParameter<edm::InputTag>("BeamSpotLabel"))),
    vtxToken_                         (consumes<std::vector<Vertex> >                (iConfig.getParameter<edm::InputTag>("vtxLabel"))),
    rhoToken_                         (consumes<double>                                (iConfig.getParameter<edm::InputTag>("rhoLabel"))),

    fixedGridRhoToken_                (consumes<double>                                (iConfig.getParameter<edm::InputTag>("rhoLabelCN"))),
    puinfoToken_                      (consumes<std::vector<PileupSummaryInfo> >       (iConfig.getParameter<edm::InputTag>("PUInfoLabel"))),
    geneventToken_                    (consumes<GenEventInfoProduct>                   (iConfig.getParameter<edm::InputTag>("generatorLabel"))),
    lheeventToken_                     (consumes<LHEEventProduct>                        (iConfig.getParameter<edm::InputTag>("lheevent"))),
    genparticleToken_                 (consumes<reco::GenParticleCollection>           (iConfig.getParameter<edm::InputTag>("genPartsLabel"))),

    pfcToken_                         (consumes<pat::PackedCandidateCollection>                    (iConfig.getParameter<edm::InputTag>("pfcLabel"))),
    jetToken_                         (consumes<pat::JetCollection>                    (iConfig.getParameter<edm::InputTag>("JetLabel"))),
    muonToken_                        (consumes<pat::MuonCollection>                   (iConfig.getParameter<edm::InputTag>("MuonLabel"))),
    displacedStandAloneMuonsToken_    (consumes<std::vector<reco::Track>>                 (iConfig.getParameter<edm::InputTag>("displacedStandAloneMuonsLabel"))),
    tauToken_                        (consumes<pat::TauCollection>                   (iConfig.getParameter<edm::InputTag>("TauLabel"))),
    electronToken1_                    (consumes<pat::ElectronCollection>             (iConfig.getParameter<edm::InputTag>("ElectronLabel"))),
    electronToken_                    (consumes<edm::View<pat::Electron> >             (iConfig.getParameter<edm::InputTag>("ElectronLabel"))),
    convToken_                        (consumes< std::vector<reco::Conversion> >                   (iConfig.getParameter<edm::InputTag>("convLabel"))),
    metToken_                         (consumes<pat::METCollection>                    (iConfig.getParameter<edm::InputTag>("METLabel"))),
    triggerToken_                     (consumes<edm::TriggerResults>                   (iConfig.getParameter<edm::InputTag>("HLTResultsLabel"))),
    filterToken_                     (consumes<edm::TriggerResults>                   (iConfig.getParameter<edm::InputTag>("filterResultsLabel")))

    {
        Sample                         = iConfig.getUntrackedParameter<std::string>("SampleLabel") ;
        IT_muon                        = iConfig.getParameter<edm::InputTag>("MuonLabel") ;
        IT_displacedStandAloneMuons    = iConfig.getParameter<edm::InputTag>("displacedStandAloneMuonsLabel") ;
        IT_electron                    = iConfig.getParameter<edm::InputTag>("ElectronLabel") ;
        IT_tau                         = iConfig.getParameter<edm::InputTag>("TauLabel") ;
        IT_htt                         = iConfig.getParameter<edm::InputTag>("L1httLabel");
        IT_jet                         = iConfig.getParameter<edm::InputTag>("JetLabel");
        IT_pfmet                       = iConfig.getParameter<edm::InputTag>("METLabel")  ;
        IT_beamspot                    = iConfig.getParameter<edm::InputTag>("BeamSpotLabel");
        IT_hltresults                  = iConfig.getParameter<edm::InputTag>("HLTResultsLabel");
        IT_genParts                    = iConfig.getParameter<edm::InputTag>("genPartsLabel");
    }


void HNL::beginJob()
{

    Nvtx      = fs->make<TH1F>("N_{vtx}"        , "Number of vertices;N_{vtx};events / 1"  ,    40, 0., 40.);
    _hCounter = fs->make<TH1D>("hCounter", "Events counter", 5,0,5);
    
    flavComp = fs->make<TH1D>("flavComp", "Flavor; flavor", 10,0,10);
    
    outputTree = new TTree("fakeTree","fakeTree");
    bookTree();
    
    _leptonP4 = new TClonesArray("TLorentzVector", nLeptonsMax);
    for (int i=0; i!=nLeptonsMax; ++i) {
        new ( (*_leptonP4)[i] ) TLorentzVector();
    }
    
    _jetP4 = new TClonesArray("TLorentzVector", nJetsMax);
    for (int i=0; i!=nJetsMax; ++i) {
        new ( (*_jetP4)[i] ) TLorentzVector();
    }
    
    GPM = GenParticleManager();
    
    
    bool isData = !(Sample=="ElectronsMC");
    if (isData)
        fMetCorrector = new OnTheFlyCorrections("Spring16_25nsV6_DATA", isData); //isData = true
    else
        fMetCorrector = new OnTheFlyCorrections("Spring16_25nsV6_MC", isData); //isData = true
    _corrLevel = "L3Absolute";
    if (isData) _corrLevel = "L2L3Residual";
    
    // pt < 15
    looseMVA[0][0] = 0.77;
    looseMVA[1][0] = 0.55;
    looseMVA[2][0] = 0.48;
    
    // pt > 25
    looseMVA[0][1] = 0.51;
    looseMVA[1][1] = 0.10;
    looseMVA[2][1] = -0.01;
    
    //A for pT < 15 GeV, B for pT > 25 GeV, and linear slope C = (A-B)/10 in between
    //cut = std::min( A, std::max( B, A - C*(pT-15)))
    //cut = std::min( looseMVA[etak][0], std::max( looseMVA[etak][1], looseMVA[etak][0] - (looseMVA[etak][0] - looseMVA[etak][1])/10.*(pt-15));
    
    _nEventsTotal = 0;
    _nEventsFiltered = 0;
    _nEventsTotalCounted = 0;
    
    _miniisocut[0] = 0.16;
    _ptratiocut[0] = 0.76;
    _ptrelcut[0] = 7.2;
    
    _miniisocut[1] = 0.2;
    _ptratiocut[1] = 0.69;
    _ptrelcut[1] = 6.0;
    
    _multiConst[0][0] = 0.2;
    _multiConst[0][1] = 0.69;
    _multiConst[0][2] = 6.0;
    
    _multiConst[1][0] = 0.16;
    _multiConst[1][1] = 0.76;
    _multiConst[1][2] = 7.2;
    
    _multiConst[2][0] = 0.12;
    _multiConst[2][1] = 0.80;
    _multiConst[2][2] = 7.2;
    
    firstEvent_ = true;
    
    reader[0] = new TMVA::Reader( "!Color:!Silent" );
    reader[1] = new TMVA::Reader( "!Color:!Silent" );
    
    for (int i=0; i!=2; ++i) {
        reader[i]->AddVariable( "LepGood_pt", &LepGood_pt );
        reader[i]->AddVariable( "LepGood_eta", &LepGood_eta );
        reader[i]->AddVariable( "LepGood_jetNDauChargedMVASel", &LepGood_jetNDauChargedMVASel );
        reader[i]->AddVariable( "LepGood_miniRelIsoCharged", &LepGood_miniRelIsoCharged );
        reader[i]->AddVariable( "LepGood_miniRelIsoNeutral", &LepGood_miniRelIsoNeutral );
        reader[i]->AddVariable( "LepGood_jetPtRelv2", &LepGood_jetPtRelv2 );
        reader[i]->AddVariable( "min(LepGood_jetPtRatiov2,1.5)", &LepGood_jetPtRatio );
        reader[i]->AddVariable( "max(LepGood_jetBTagCSV,0)", &LepGood_jetBTagCSV );
        reader[i]->AddVariable( "LepGood_sip3d", &LepGood_sip3d );
        reader[i]->AddVariable( "log(abs(LepGood_dxy))", &LepGood_dxy );
        reader[i]->AddVariable( "log(abs(LepGood_dz))", &LepGood_dz );
    }
    reader[0]->AddVariable( "LepGood_mvaIdSpring16GP", &LepGood_mvaIdSpring15 );
    reader[1]->AddVariable( "LepGood_segmentCompatibility", &LepGood_segmentCompatibility );
    
    //reader[0]->BookMVA( "BDTG method", "forMoriond16_el_sigTTZ_bkgTT_BDTG.weights.xml" );
    //reader[1]->BookMVA( "BDTG method", "forMoriond16_mu_sigTTZ_bkgTT_BDTG.weights.xml" );

    reader[0]->BookMVA( "BDTG method", "el_BDTG.weights.xml" );
    reader[1]->BookMVA( "BDTG method", "mu_BDTG.weights.xml" );

    
    _mvaCut[0] = 0.5;
    _mvaCut[1] = -0.2;
    
}

void HNL::endJob() {
    //outputTree -> Write();
    // store nEventsTotal and nEventsFiltered in preferred way
    std::cout<<_nEventsTotal<<std::endl;
    std::cout<<_nEventsFiltered<<std::endl;
    
    delete fMetCorrector;
}

void HNL::analyze(const edm::Event& iEvent, const edm::EventSetup& iEventSetup){
    //if (iEvent.id().event()!=100) return;
    // cout << "--------------------------------------------------------------------------------------------------------------" << endl;
    // cout << "---zhud: doing HNL::analyze() with event: " << iEvent.id().event() << endl;;
    //bool islepton;
    if (Sample=="ElectronsMC") {
        _genHT = 0;


        //******************************************************************************************************************
        // begin Gen level particles                  ****************************************************************************
        //******************************************************************************************************************
        //
        
            //bool exists = false;

            //initialize counters
            _nGenLep = 0;
            _nGenE = 0;
            _nGenMu = 0;
            _nGenStatusNot1 = 0;
            _nGenHNL = 0;
            _nGenHNLMu = 0;
            _nlHNLMu = 0;

            _nMu=0;
            
            iEvent.getByToken(genparticleToken_, TheGenParticles);
            std::vector<const GenParticle*> vGenElectrons, vGenMuons, vGenNPElectrons, vGenNPMuons, vGenW;
            
            if(TheGenParticles.isValid()) {
                GPM.SetCollection(TheGenParticles);
                GPM.Classify();
                vGenMuons = GPM.filterByStatus(GPM.getPromptMuons(),1);
                vGenElectrons = GPM.filterByStatus(GPM.getPromptElectrons(),1);
                vGenNPMuons = GPM.filterByStatus(GPM.getNonPromptMuons(),1);
                vGenNPElectrons = GPM.filterByStatus(GPM.getNonPromptElectrons(),1);
                //std::cout<<"*************"<<std::endl;
                
                TLorentzVector Gen0;
                Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
                
                TLorentzVector Gen0s;
                Gen0s.SetPtEtaPhiE( 0, 0, 0, 0);
                
                TLorentzVector Gen0s2;
                Gen0s2.SetPtEtaPhiE( 0, 0, 0, 0);
                
                int cnt = 0;
                double leadPt = 0;
                

                //============ Pat Muons ============
                edm::Handle< std::vector<pat::Muon> > thePatMuons;
                iEvent.getByToken( muonToken_, thePatMuons );
                if( ! thePatMuons.isValid() )  ERR(IT_muon) ;
                //==================================


                //============ Pat displacedStandAloneMuons ============
                edm::Handle< std::vector<reco::Track> > theDisplacedStandAloneMuons;
                iEvent.getByToken( displacedStandAloneMuonsToken_, theDisplacedStandAloneMuons );
                if( ! theDisplacedStandAloneMuons.isValid() )  ERR(IT_displacedStandAloneMuons) ;
                //======================================================

                int GenParticlecounter = 0;            
                


                for(GenParticleCollection::const_reverse_iterator q = TheGenParticles->rbegin() ; q != TheGenParticles->rend() ; q++ ){
                   if (q->pdgId() == 9900012){
                        // cout << "Found 1 HNL! ";
                        _nGenHNL++;
                        
                        _GenHNLP              = q->p();
                        _GenHNLEta            = q->eta();
                        _GenHNLPhi            = q->phi();
                        _GenHNLE              = q->energy();
                        _GenHNLMass           = q->mass();
                        _GenHNLBetaGamma      = _GenHNLP / _GenHNLMass;
                        _GenHNLGamma          = _GenHNLE  / _GenHNLMass;
                        _GenHNLBeta           = _GenHNLBetaGamma / _GenHNLGamma;

                        // cout << "Mass = "<< _GenHNLMass << "; Gamma = " << _GenHNLGamma << endl; 
                                 
                        _GenHNLVxProd         = q->vx();
                        _GenHNLVyProd         = q->vy();
                        _GenHNLVxyProd        = abs(sqrt(pow(q->vx(),2) + pow(q->vy(),2)));
                        _GenHNLVzProd         = q->vz();
                        _GenHNLVxyzProd       = abs(sqrt(pow(q->vx(),2) + pow(q->vy(),2) + pow(q->vz(),2)));
                    } 
                }
                
                std::vector<const pat::Muon* > sMu = SelectAllPatMuons( *thePatMuons, _minPt0, PV, _looseD0Mu, true);
                std::vector<const reco::Track* > sDSAMu = SelectAllTrackMuons( *theDisplacedStandAloneMuons);

                for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
                    
                    //cout << "------------------------------------------------------------------"<< endl;
                    // cout << "---zhud: starting for loop through all TheGenParticles; Number = " << GenParticlecounter << endl;
    		        int id = TMath::Abs(p->pdgId());
                    // std::cout<<"gen particle pdg ID "<<id<<std::endl;
                    //if (id == 1000023) mChi20 = p->mass(); // SUSY stuff
                    //if (id == 1000022) mChi10 = p->mass(); // SUSY stuff


                    if ( (id == 12 || id == 14 || id == 16 ) && (p->status() == 1) ) { // 12 = nu_e; 14 = nu_mu; 16 = nu_tau
                        TLorentzVector Gen;
                        //cout << "---zhud: id = 12 | 14 | 16; doing Gen.SetPtEtaPhiE( p->pt(), p->eta(), p->phi(), p->energy() ): ";
                        Gen.SetPtEtaPhiE( p->pt(), p->eta(), p->phi(), p->energy() );
                        Gen0 += Gen;
                    }
                    if ( (id == 11 || id == 13  ) && (p->status() == 1) ) { // 11 = e, 13 = mu
                        
                        _nGenLep++;
                        if (p->status()!=1){_nGenStatusNot1++;}
                        if (id == 11 && p->status()==1){_nGenE++;}
                        if (id == 13 && p->status()==1){
                            
                            _isdetectedMu[_nGenMu]=false;
                            _isfromHNLMu[_nGenMu]=false;
                            _ispromptMu[_nGenMu]=false;

                            _isdetectedDSAMu[_nGenMu]=false;
                            // cout << "---zhud: doing Gen Muon number " << GenParticlecounter << ": ";
                            
                            _GenMuPt[_nGenMu] = p->pt();
                            _GenMuEta[_nGenMu] = p->eta();
                            _GenMuPhi[_nGenMu] = p->phi();
                            _GenMuE[_nGenMu] = p->energy();
                            
                            _GenMuVx[_nGenMu]   = p->vx();
                            _GenMuVy[_nGenMu]   = p->vy();
                            _GenMuVxy[_nGenMu]  = abs(sqrt(pow(p->vx(),2) + pow(p->vy(),2)));
                            _GenMuVz[_nGenMu]   = p->vz();
                            _GenMuVxyz[_nGenMu] = abs(sqrt(pow(p->vx(),2) + pow(p->vy(),2) + pow(p->vz(),2)));


                            // cout << "pt = " << p->pt() << " eta = " << p->eta() << endl;
                            
                            // cout << "---zhud: vx = " << p->vx() << " vy = " << p->vy() << " vz = " << p->vz() << endl;


                            
                            TLorentzVector vGen, vReco;
                            vGen.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
                            double deltaR = 9999.;
			                double deltaRDSA = 9999.;

                            for(unsigned int i = 0 ; i < sMu.size() ;i++ ){
                                const pat::Muon *iM = sMu[i];
                                vReco.SetPtEtaPhiE(iM->pt(),iM->eta(),iM->phi(),iM->energy());
                                double deltaRcur = vGen.DeltaR(vReco);
                                if (deltaRcur <deltaR){
                                    deltaR = deltaRcur;
                                }    
                            }

            			    for(unsigned int i = 0; i<sDSAMu.size(); i++){
                			    const reco::Track *iDSAMu = sDSAMu[i];
                				vReco.SetPtEtaPhiE(iDSAMu->pt(),iDSAMu->eta(),iDSAMu->phi(),0); // energy by default set to 0 since energy is not a function of <reco::Track>
                				double deltaRcur = vGen.DeltaR(vReco);
                				if(deltaRcur < deltaRDSA){
                					deltaRDSA = deltaRcur;
                				}
            			    }
                            
                            if (deltaR < 0.2){
                                _foundGenMuPt[_nGenMu] = p->pt();
                                _isdetectedMu[_nGenMu] = true;
                                // cout << "---zhud: Found a corresponding Reco Muon!!!" << endl;      
                            }

                            
                            if (deltaR >= 0.2){
                                _foundGenMuPt[_nGenMu] = -1;
                                _isdetectedMu[_nGenMu] = false;
                                //cout << "---zhud: NO corresponding Reco Muon founded..." << endl;
                            }

                            if (deltaRDSA < 0.2){
                                _isdetectedDSAMu[_nGenMu] = true;
                                // cout << "---zhud: Found a corresponding Reco DSA Muon!!!" << endl; 
                            }

                            if (deltaRDSA >= 0.2){
                                _isdetectedDSAMu[_nGenMu] = false;
                            }

                            
                            //GPM.printInheritance(&*p);

                            if (GPM.fromID(&*p,9900012) && p->pt()>1) {
                                
                                _GenHNLMuPt[_nGenHNLMu] = p->pt();
                                _GenHNLMuEta[_nGenHNLMu] = p->eta();
                                _GenHNLMuPhi[_nGenHNLMu] = p->phi();
                                _GenHNLMuE[_nGenHNLMu] = p->energy();

                                _isfromHNLMu[_nGenMu] = true;
                                
                                // cout<< "---zhud: Found Muon<-HNL"<<endl;
                                
                                
                                _GenHNLVxDecay[_nGenHNLMu]          = _GenMuVx[_nGenMu] - _GenHNLVxProd;
                                _GenHNLVyDecay[_nGenHNLMu]          = _GenMuVy[_nGenMu] - _GenHNLVyProd;
                                _GenHNLVxyDecay[_nGenHNLMu]         = sqrt(pow(_GenMuVx[_nGenMu],2) + pow(_GenMuVy[_nGenMu],2));
                                _GenHNLVzDecay[_nGenHNLMu]          = _GenMuVz[_nGenMu] - _GenHNLVzProd;
                                _GenHNLVxyzDecay[_nGenHNLMu]        = sqrt(pow(_GenMuVx[_nGenMu],2) + pow(_GenMuVy[_nGenMu],2) + pow(_GenMuVz[_nGenMu],2));
                                _GenHNLMu3DDisplacement[_nGenHNLMu] = _GenHNLVxyzDecay[_nGenHNLMu] / _GenHNLBetaGamma;
                                
                                // cout << "---zhud: 3D Displacement: " << _GenHNLVxyzDecay[_nGenHNLMu] << " / " << _GenHNLBetaGamma << " = " << _GenHNLMu3DDisplacement[_nGenHNLMu] << endl;

                                _nGenHNLMu++;    
                            }

                            if (GPM.fromID(&*p,9900012)==false && p->isPromptFinalState() == true){
                                _ispromptMu[_nGenMu] = true;
                                //cout << "---zhud: Found Muon<-PROMPT" <<endl;
                            }
                        
                            _nGenMu++;
                        }

                    }
                    

                    if ( ( id == 15 ) && (p->status() == 2) ) { //15 = tau
                        //cout<< "---zhud: id = 15; ";
    		            // std::cout<<"---pt = "<<p->pt()<<std::endl;
                        //GPM.printInheritance(&*p);
                    }
                    //if ( (id == 6 || id == 24 || id == 23 || id == 25 ) && (p->isLastCopy()) && (p->status() == 62) ) {
                    if ( (id == 23) && (p->isLastCopy()) && (p->status() == 62) ) { // 23 = Z_0
                        //cout << "---zhud: id = 23; "<< endl;
                        TLorentzVector Gen;
                        Gen.SetPtEtaPhiE( p->pt(), p->eta(), p->phi(), p->energy() );
                        Gen0s += Gen;
                        cnt++;
                        //std::cout<<id<<" "<<p->status()<<" "<<p->pt()<<std::endl;
                        if (id == 6) {
                            Gen0s2 += Gen;
                        }
                    }
                    //zhud: Here is the most printed output line in terminal
    		        //std::cout<<"---Most printed version: All: id = "<<id<<"; status = "<<p->status()<<"; last copy = "<<p->isLastCopy()<<"; mass = "<<p->mass()<<std::endl;
                    if (id == 22 && (p->status() == 23 || p->status() == 1)) { // 22 = gamma
                        if (p->pt() > leadPt) {
                            // cout<<"---zhud: id = 22; "<< endl;
                            _genPhotOrigin = photonOrigin(&*p);
                            leadPt = p->pt();
                        }
                    }
                    
                    if ((id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 21 || id == 22 ) && (p->status() == 23)){
                        _genHT += p->pt();
                    }
                    
            
                    GenParticlecounter ++;
                    //cout <<"------------------------------------------------------------------"<<endl;
                }
                

                if (Gen0.E()!=0) {
                    _genmet = Gen0.Pt();
                    _genmet_phi = Gen0.Phi();
                } else {
                    _genmet = 0;
                    _genmet_phi = 0;
                }
                
                _ptSystem = Gen0s.Pt();
                _ptTTSystem = Gen0s2.Pt();
                
                
                //std::cout<<cnt<<" particles: "<<_ptSystem<<" "<<_ptTTSystem<<std::endl;
                //std::cout<<"masses "<<mChi20<<" "<<mChi10<<std::endl;           
            }
        //std::cout<<"photon "<<int(exists)<<std::endl;
        //******************************************************************************************************************
        // end Gen Level Particles *****************************************************************************************
        //******************************************************************************************************************
        
        //**************************************************************************************
        // MC
        //**************************************************************************************
    }
    

    
    
    
    //============ Begin Counter ============    
    //============ Total number of events is the sum of the events ============
    //============ in each of these luminosity blocks ============
        _runNb = iEvent.id().run();
        _eventNb = iEvent.id().event();
        _lumiBlock = iEvent.luminosityBlock();
        _weight = 1;
        
        if (Sample=="ElectronsMC") {    
            edm::Handle<GenEventInfoProduct> pdfvariables;
            iEvent.getByToken(geneventToken_, pdfvariables);
            _weight=pdfvariables->weight();
            
            if (_weight > 0) _weight = 1;
            else _weight = -1;
            
            edm::Handle<std::vector<PileupSummaryInfo> > pileupInfo;
            iEvent.getByToken(puinfoToken_, pileupInfo); 
            //iEvent.getByToken("AddPileupInfo", pileupInfo);
            std::vector<PileupSummaryInfo>::const_iterator PVI;
            
            //std::cout<<"got pu"<<std::endl;
            for (PVI = pileupInfo->begin(); PVI !=pileupInfo->end(); ++PVI) {
                if( PVI->getBunchCrossing() == 0 ) { // in-time PU
                    _n_Interactions     = PVI->getPU_NumInteractions();
                    _n_trueInteractions = PVI->getTrueNumInteractions();
                }
            }
            //std::cout<<"got pu "<<_n_Interactions<<" "<<_n_trueInteractions<<std::endl;        
        }

        _nEventsTotalCounted+=_weight;
        _hCounter->Fill(0.,_weight);
    //============ End Counter ============


    filterbadChCandidate = false;
    filterbadPFMuon = false;
    if (Sample!="ElectronsMC") {
    
        Handle<bool> ifilterbadChCand;
        iEvent.getByToken(ifilterbadChCandToken, ifilterbadChCand);
        filterbadChCandidate = *ifilterbadChCand;
        
        Handle<bool> ifilterbadPFMuon;
        iEvent.getByToken(ifilterbadPFMuonToken, ifilterbadPFMuon);
        filterbadPFMuon = *ifilterbadPFMuon;
    }

    //============ Beamspot ============
        edm::Handle< reco::BeamSpot > theBeamSpot;
        iEvent.getByToken( beamspotToken_, theBeamSpot );
        if( ! theBeamSpot.isValid() ) ERR( IT_beamspot ) ;
        BeamSpot::Point  BS= theBeamSpot->position();
    //==================================
    
    //============ Primary vertices ============
        edm::InputTag IT_goodVtx = edm::InputTag("offlineSlimmedPrimaryVertices");
        edm::Handle<std::vector<Vertex> > theVertices;
        iEvent.getByToken( vtxToken_, theVertices) ;
        if( ! theVertices.isValid() ) ERR(IT_goodVtx ) ;
        int nvertex = theVertices->size();
        
        _n_PV = nvertex;
        
        Nvtx->Fill(TMath::Min(nvertex,39));
        if(! nvertex ){
            cout << "[WARNING]: No candidate primary vertices passed the quality cuts, so skipping event" << endl;
            return ;
        }
        
        PV = theVertices->begin()->position();
        const Vertex* PVtx = &((*theVertices)[0]);
        _PVchi2 = PVtx->chi2();
        _PVerr[0] = PVtx->xError();
        _PVerr[1] = PVtx->yError();
        _PVerr[2] = PVtx->zError();
    //==================================
    
    for (int i=0; i!=nLeptonsMax; ++i) {
        _tauDeltaR[i] = 9999.;
    }
    
    //============= trigger ==============
        edm::Handle<TriggerResults> trigResults;
        iEvent.getByToken(triggerToken_, trigResults);
        
        for (int i=0; i!=20; ++i) {
            _triggersPost[i] = 0;
        }
        
        for (int i=0; i!=5; ++i) {
            for (int j=0; j!=5; ++j) { _triggersCS[j][i] = 0;}
        }
        for (int i=0; i!=8; ++i) {
            _triggers1l[i] = 0;
        }
        for (int i=0; i!=6; ++i) {
            _triggersTau[i] = 0;
            for (int j=0; j!=2; ++j) _triggers2l[j][i] = 0;
            for (int j=0; j!=2; ++j) _triggers2lbkp[j][i] = 0;
        }
        
        for (int i=0; i!=2; ++i) _triggersCSb[i] = 0;
        
        if( trigResults.failedToGet() && firstEvent_) cout << "--- NO TRIGGER RESULTS !! ---" << endl;
        if( !trigResults.failedToGet() ) {
            
            unsigned int n_Triggers = trigResults->size();
            
            const edm::TriggerNames & triggerNames = iEvent.triggerNames(*trigResults);
            
            // if ( firstEvent_ ) {
            //     edm::TriggerNames::Strings allTriggers( triggerNames.triggerNames() );
            //     std::cout << "--- Trigger Menu --- " << std::endl;
            //     for ( unsigned int i_Name = 0; i_Name < n_Triggers; ++i_Name ) {
            //         std::cout << allTriggers.at( i_Name ) << std::endl;
            //     }
            //     std::cout << "-------------------- " << std::endl;
            //     firstEvent_ = false;
            // }
            
            for( unsigned int i_Trig = 0; i_Trig < n_Triggers; ++i_Trig ) {
                
                if (trigResults.product()->accept(i_Trig)) {
                    
                    TString TrigPath = triggerNames.triggerName(i_Trig);
                    //std::cout<<TrigPath<<std::endl;
                    for (int i=0; i!=6; ++i) {
                        if (TrigPath.Contains(_triggersTauNames8[i]))
                            _triggersTau[i] = 1;
                        for (int j=0; j!=2; ++j) {
                            if (TrigPath.Contains(_triggers2lNames8[j][i]))
                                _triggers2l[j][i] = 1;
                        }
                        for (int j=0; j!=2; ++j) {
                            if (TrigPath.Contains(_triggers2lNames8Bkp[j][i]))
                                _triggers2lbkp[j][i] = 1;
                        }
                    }
                    
                    for (int i=0; i!=5; ++i) {
                        for (int j=0; j!=5; ++j) {
                            if (TrigPath.Contains(_triggersCSNames8[j][i]))
                                _triggersCS[j][i] = 1;
                        }
                    }
                    for (int i=0; i!=8; ++i) {
                        if (TrigPath.Contains(_triggers1lNames8[i]))
                            _triggers1l[i] = 1;
                    }
                    for (int i=0; i!=2; ++i) {
                        if (TrigPath.Contains(_triggersCSbNames8[i]))
                            _triggersCSb[i] = 1;
                    }
                    for (int i=0; i!=20; ++i) {
                        if (TrigPath.Contains(_postICHEP8[i]))
                            _triggersPost[i] = 1;
                    }
                }
            }
        }
    //==================================


    //============= filters ==============
        edm::Handle<TriggerResults> filtResults;
        iEvent.getByToken(filterToken_, filtResults);
        
        Flag_eeBadScFilter = 0;
        Flag_eeBadScFilter = 0;
        Flag_HBHENoiseFilter = 0;
        Flag_HBHENoiseIsoFilter = 0;
        Flag_globalTightHalo2016Filter = 0;
        Flag_CSCTightHalo2015Filter = 0;
        Flag_EcalDeadCellTriggerPrimitiveFilter = 0;
        Flag_goodVertices = 0;
        if( filtResults.failedToGet() ) cout << "--- NO FILTER RESULTS !! ---" << endl;
        if( !filtResults.failedToGet() ) {
            unsigned int n_Triggers = filtResults->size();
            const edm::TriggerNames & triggerNames = iEvent.triggerNames(*filtResults);
            for( unsigned int i_Trig = 0; i_Trig < n_Triggers; ++i_Trig ) {
                if (filtResults.product()->accept(i_Trig)) {
                    
                    TString TrigPath = triggerNames.triggerName(i_Trig);
                    //std::cout<<TrigPath<<std::endl;
                    
                    if (TrigPath.Contains("Flag_eeBadScFilter"))
                        Flag_eeBadScFilter = true;
                    if (TrigPath.Contains("Flag_HBHENoiseFilter"))
                        Flag_HBHENoiseFilter = true;
                    if (TrigPath.Contains("Flag_HBHENoiseIsoFilter"))
                        Flag_HBHENoiseIsoFilter = true;
                    if (TrigPath.Contains("Flag_globalTightHalo2016Filter"))
                        Flag_globalTightHalo2016Filter = true;
                    if (TrigPath.Contains("Flag_CSCTightHalo2015Filter"))
                        Flag_CSCTightHalo2015Filter = true;
                    if (TrigPath.Contains("Flag_EcalDeadCellTriggerPrimitiveFilter"))
                        Flag_EcalDeadCellTriggerPrimitiveFilter = true;
                    if (TrigPath.Contains("Flag_goodVertices"))
                        Flag_goodVertices = true;
                }
            }
        }
    //=================================
    
   
    //============ Pat MET ============
        edm::Handle< vector<pat::MET> > ThePFMET;
        iEvent.getByToken(metToken_, ThePFMET);
        if( ! ThePFMET.isValid() ) ERR( IT_pfmet );
        const vector<pat::MET> *pfmetcol = ThePFMET.product();
        const pat::MET *pfmet;
        pfmet = &(pfmetcol->front());
        _met = pfmet->pt();
        _met_phi = pfmet->phi();
        //if (_met < 50) return;
        //std::cout<<"raw met "<<_met<<std::endl;
    //==================================
    
    //============ Pat Jets ============
        edm::Handle< std::vector< pat::Jet> > thePatJets;
        iEvent.getByToken(jetToken_ , thePatJets );
        if( ! thePatJets.isValid() ) ERR(IT_jet);
    //==================================
    
    //============ PF cand ============
        edm::Handle<pat::PackedCandidateCollection> pfcands;
        iEvent.getByToken(pfcToken_, pfcands);
    //==================================
    
    //============ Pat Muons ============
        edm::Handle< std::vector<pat::Muon> > thePatMuons;
        iEvent.getByToken( muonToken_, thePatMuons );
        if( ! thePatMuons.isValid() )  ERR(IT_muon) ;
    //==================================
    
    //============ Pat displacedStandAloneMuons ============
    edm::Handle< std::vector<reco::Track> > theDisplacedStandAloneMuons;
    iEvent.getByToken( displacedStandAloneMuonsToken_, theDisplacedStandAloneMuons );
    if( ! theDisplacedStandAloneMuons.isValid() )  ERR(IT_displacedStandAloneMuons) ;
    //======================================================

    //============ Pat Electrons ============
        edm::Handle< std::vector<pat::Electron> > thePatElectrons;
        iEvent.getByToken( electronToken1_, thePatElectrons );
        if( ! thePatElectrons.isValid() ) ERR( IT_electron );
        //==================================
        edm::Handle< edm::View<pat::Electron> > thePatElectronsView;
        iEvent.getByToken( electronToken_, thePatElectronsView );
        //if( ! thePatElectrons.isValid() ) ERR( IT_electron );
    //==================================
    
    //============ Pat Taus ============
        edm::Handle< std::vector<pat::Tau> > thePatTaus;
        iEvent.getByToken( tauToken_, thePatTaus );
        if( ! thePatTaus.isValid() )  ERR(IT_tau) ;
    //==================================
    
    
    //============ Rho ============
        edm::Handle<double> rhoJets;
        iEvent.getByToken(rhoToken_ , rhoJets);//kt6PFJets//edm::InputTag("fixedGridRhoFastjetAll","")
        myRhoJets = *rhoJets;
    //==================================
    
    //============ RhoNC ============
        edm::Handle<double> rhoJetsNC;
        iEvent.getByToken(fixedGridRhoToken_ , rhoJetsNC);//kt6PFJets //edm::InputTag("fixedGridRhoFastjetCentralNeutral","")
        myRhoJetsNC = *rhoJetsNC;
    //==================================
    
    
    //============ Conversions ============
        edm::Handle< std::vector<reco::Conversion> > theConversions;
        iEvent.getByToken(convToken_, theConversions);
    //==================================
    
    // Get MVA values and categories (optional)
    //============ MVA ============
        edm::Handle<edm::ValueMap<float> > mvaValues;
        edm::Handle<edm::ValueMap<int> > mvaCategories;
        iEvent.getByToken(mvaValuesMapToken_,mvaValues);
        iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);
        
        edm::Handle<edm::ValueMap<float> > mvaValuesHZZ;
        edm::Handle<edm::ValueMap<int> > mvaCategoriesHZZ;
        iEvent.getByToken(mvaValuesHZZMapToken_,mvaValuesHZZ);
        iEvent.getByToken(mvaCategoriesHZZMapToken_,mvaCategoriesHZZ);    
    //==================================
    
    enum decay {//origindetailed
        W_L,  // 0
        W_T_L, // 1
        W_B_L, // 2
        W_B_D_L, //3
        W_B_D_T_L, // 4
        W_B_T_L, // 5
        W_D_L, // 6
        W_D_T_L, //7
        B_L, // 8
        B_D_L, //9
        B_D_T_L, //10
        B_T_L,  // 11
        D_L, //12
        D_T_L, //13
        B_Baryon, // 14
        C_Baryon, //15
        pi_0, //16
        photon_, //17
        F_L, //18
        N_U_L_L // 19
    };
    
    //zhud: original
    //std::vector<const pat::Muon* > sMu = ssbLooseMuonSelector( *thePatMuons, _minPt0, PV, _looseD0Mu, true);
    
    //zhud: all Reco without preselection
    std::vector<const pat::Muon* > sMu = SelectAllPatMuons( *thePatMuons, _minPt0, PV, _looseD0Mu, true);
    //std::vector<const pat::Muon* > sDSAMu = SelectAllPatMuons( *theDisplacedStandAloneMuons, _minPt0, PV, _looseD0Mu, true);
    std::vector<const reco::Track* > sDSAMu = SelectAllTrackMuons(*theDisplacedStandAloneMuons);



    std::vector<const pat::Electron* > sEl = ssbMVAElectronSelector( *thePatElectrons, _minPt0, PV, _looseD0E, _chargeConsistency, theConversions, BS,
                                                                    false, false);
    
    SelectedJetsAll = JetSelectorAll(*thePatJets, 5., 3.0);
    //std::vector<const pat::Jet* > SelectedJets = JetSelector(*thePatJets, _jetPtCut, _jetEtaCut);
    std::vector<const pat::Jet* > SelectedJets = JetSelector(*thePatJets, 15, _jetEtaCut);
    
    //std::cout<<sEl.size() + sMu.size()<<std::endl;
    
    //zhud: for efficiency, comment this line out to also have events with < 3
    //if (sEl.size() + sMu.size() + sTau.size() < 3) return;
    //std::cout<<sMu.size()<<" "<<sEl.size()<<std::endl;
    


    int leptonCounter = 0;
    _nDSAMu = 0;
    
    //***********************************************************
    //***zhud: Analyzing Displaced Stand Alone Muons*************
    //***********************************************************
    // cout << "---zhud: Doing DSA Muons"<<endl; 
    for(unsigned int i=0; i<sDSAMu.size(); i++){
        const reco::Track *iM = sDSAMu[i];
        // cout<<"---zhud: DSA 1: pt = " << iM->pt() <<endl;
        _nDSAMu++;
    }


    //***********************************************************
    //***zhud: Finished Analyzing Displaced Stand Alone Muons****
    //***********************************************************

    

    //**************************************************
    //**********zhud: Analyzing Reco Muons**************
    //**************************************************
    for(unsigned int i = 0 ; i < sMu.size() ;i++ ){
        //for (std::vector<pat::Muon>::const_iterator mu = thePatMuons->begin() ; mu != thePatMuons->end() ; mu++){
        //const pat::Muon *iM = &*mu;
        const pat::Muon *iM = sMu[i];
        
        //if((iM->innerTrack()->ptError())/(iM->innerTrack()->pt()) > 0.2) continue;
        
        
        _leptonIndex = i;
        //i++;
        
        if (leptonCounter == nLeptonsMax) continue;
        
        _islooseID[leptonCounter] = true;

        _dptoverpt[leptonCounter] = (iM->innerTrack()->ptError())/(iM->innerTrack()->pt());
        
        _flavors[leptonCounter] = 1;
        _charges[leptonCounter] = iM->charge();
        _isolation[leptonCounter][0] = pfRelIso(iM,myRhoJetsNC);
        _isolation[leptonCounter][1] = pfRelIso(iM);
        
        _ipPV[leptonCounter] = TMath::Abs(iM->innerTrack()->dxy(PV));
        _ipPVerr[leptonCounter] = iM->innerTrack()->dxyError();
        
        _ipZPV[leptonCounter] = iM->innerTrack()->dz(PV);
        _ipZPVerr[leptonCounter] = iM->innerTrack()->dzError();
        
        _3dIP[leptonCounter]    = iM->dB(pat::Muon::PV3D);
        _3dIPerr[leptonCounter] = iM->edB(pat::Muon::PV3D);
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);
        
        _mvaValue[leptonCounter] = -1;
        _muonSegmentComp[leptonCounter] = iM->segmentCompatibility();
        
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Loose_Muon
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#MediumID2016_to_be_used_with_Run
        bool goodGlb = iM->isGlobalMuon() && iM->globalTrack()->normalizedChi2() < 3
        && iM->combinedQuality().chi2LocalPosition < 12 && iM->combinedQuality().trkKink < 20;
        bool good = iM->innerTrack()->validFraction() >= 0.8 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
    
        
        _istightID[leptonCounter] = good;
        
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iM->pt(), iM->eta(), iM->phi(), iM->energy());
        
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lMuPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lMuEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lMuPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lMuE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();

        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        
        _miniisolation[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, false, myRhoJetsNC);
        _miniisolation[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.3, 10., false, false, myRhoJetsNC);
        _miniisolation[leptonCounter][2] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.4, 10., false, false, myRhoJetsNC);
        
        _miniisolationCharged[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, true, myRhoJetsNC);
        _miniisolationCharged[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.3, 10., false, true, myRhoJetsNC);
        _miniisolationCharged[leptonCounter][2] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.4, 10., false, true, myRhoJetsNC);
        
        
        fillCloseJetVars(leptonCounter);
        _lepMVA[leptonCounter] = lepMVAvalue(leptonCounter);
        
        _isloose[leptonCounter] =
        ( //_istightID[leptonCounter] &&
         //_miniisolation[leptonCounter][0] < _miniIsoCut)
         (_isolation[leptonCounter][0] < _relIsoCut || _isolation[leptonCounter][0]*_lPt[leptonCounter] < _absIsoCut)
        && ( fabs(_ipPV[leptonCounter])<_looseD0Mu )
        && ( fabs(_3dIPsig[leptonCounter]) < _looseSIP)
         );
        if (_isloose[leptonCounter])
            _istight[leptonCounter] = (_miniisolation[leptonCounter][0] < _miniisocut[_flavors[leptonCounter]]
                                       && (_ptratio[leptonCounter] > _ptratiocut[_flavors[leptonCounter]]
                                           || _ptrel[leptonCounter] > _ptrelcut[_flavors[leptonCounter]]))
            && ( fabs(_3dIPsig[leptonCounter]) < _looseSIP);
        else _istight[leptonCounter] = false;
        
        if (!_isloose[leptonCounter]) continue;
        //if (!_istight[leptonCounter]) continue;
        
        
        _triggerMatch[leptonCounter] = true;
        _triggerIsoMatch[leptonCounter] = isoTriggerEmulator(iM);
        
        _istightIso[leptonCounter] = _istight[leptonCounter]
        && (_miniisolation[leptonCounter][0] < _miniisocut[_flavors[leptonCounter]]
            && (_ptratio[leptonCounter] > _ptratiocut[_flavors[leptonCounter]]
                || _ptrel[leptonCounter] > _ptrelcut[_flavors[leptonCounter]]));
        
        
        //cout<<"---zhud: found a reco muon with pt = "<< iM->pt() <<endl;

        if (Sample=="ElectronsMC") {
            //**************************************************************************************
            // begin MC
            //**************************************************************************************
            
            const GenParticle* mc = GPM.matchedMC(iM, 13);
            _originPhot[leptonCounter] = -2;
            if ( mc!=0 ) {
                fillMCVars(mc, leptonCounter);
                _ipPVmc[leptonCounter] = TMath::Abs(iM->innerTrack()->dxy(PVmc));

                //zhud: for the efficiency plots
                if (fabs(mc->pdgId()) == 13){
                    // _lMuPt[leptonCounter] = mc -> pt();
                    if (GPM.fromID(&*mc,9900012)){
                        _lHNLMuPtmc[leptonCounter] = mc->pt();
                        //cout<<"---zhud: found a reco HNL->Mu with McPt = "<< mc->pt()<<endl;
                        _nlHNLMu++;
                    }       
                }
            }           
            else {
                mc = GPM.matchedMC(iM);
                if ( mc!=0 ) {
                    fillMCVars(mc, leptonCounter);
                    if (mc->pdgId() != 22) //not a photon
                        _origin[leptonCounter] = GPM.originReduced(_originDetailed[leptonCounter]);
                    else {
                        _originPhot[leptonCounter] = photonOrigin(mc);
                    }
                } else {
                    //std::cout<<"No match mu"<<std::endl;
                    _originDetailed[leptonCounter] = -1;
                    _origin[leptonCounter] = 4;
                    _mompt[leptonCounter] = 0;
                    _momphi[leptonCounter] = 0;
                    _mometa[leptonCounter] = 0;
                    _mompdg[leptonCounter] = 0;
                }
            }
            
            //if (_origin[leptonCounter]!=0) continue;
            
            fillIsoMCVars(leptonCounter);
            if (_closeIndex[leptonCounter] >=0)
                matchCloseJet(leptonCounter);
            
            
            
            //**************************************************************************************
            // end MC *
            //**************************************************************************************
        }
        
        
        _mllZ[leptonCounter] = 9999.;
        _mllG[leptonCounter] = 9999.;
        _mllZj[leptonCounter] = 9999.;
        _mllGj[leptonCounter] = 9999.;
        _lepDeltaRloose[i] = 9999.;
        _missingHits[leptonCounter] = 0;
        leptonCounter++;
        
    }

    //**************************************************
    //******zhud: Finished Analyzing Reco Muons*********
    //**************************************************
    
    _nMu = leptonCounter;
    //std::cout<<leptonCounter<<" "<<((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt() <<std::endl;
    
    
    //for(unsigned int i = 0 ; i < sEl.size() ;i++ ){

    //**************************************************
    //**********zhud: Analyzing Reco Electrons**********
    //**************************************************

    for (size_t i = 0; i < thePatElectronsView->size(); ++i){
        const auto iE = thePatElectronsView->ptrAt(i);
        if (!ssbMVAElectronSelectorPassed( iE, _minPt0, PV, _looseD0E, _chargeConsistency, theConversions, BS, false, false, false)) continue;
        
        const reco::GsfTrackRef gsfTrack = iE->gsfTrack();
        _missingHits[leptonCounter] = gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

        
        _muonSegmentComp[leptonCounter] = -1;
        if (iE->pt() > 10)
            _mvaValue[leptonCounter] = (*mvaValues)[iE];
        else
            _mvaValue[leptonCounter] = (*mvaValuesHZZ)[iE];
        
        /*bool passed = false;
        _islooseID[leptonCounter] = passed;
        if (TMath::Abs(iE->eta()) < 0.8 ) {
            passed = _mvaValue[leptonCounter] > looseMVA[0][0];
        } else if (TMath::Abs(iE->eta()) < 1.479 ) {
            passed = _mvaValue[leptonCounter] > looseMVA[1][0];
        } else {
            passed = _mvaValue[leptonCounter] > looseMVA[2][0];
        }*/
        _islooseID[leptonCounter] = true;
        
        //if (!passed) continue;
        
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iE->pt(), iE->eta(), iE->phi(), iE->energy());

        double minDeltaR = 9999;
        for (int j=0; j!=_nMu; ++j) {
            if (fabs(_ipPV[leptonCounter]) > 0.05 || fabs(_3dIPsig[leptonCounter]) > 4) continue;
            if (_miniisolation[leptonCounter][0] > 0.4) continue;
            _lepDeltaR[leptonCounter][j] = ((TLorentzVector *)_leptonP4->At(j))->DeltaR( *((TLorentzVector *)_leptonP4->At(leptonCounter)) );
            if (_lepDeltaR[leptonCounter][j] < minDeltaR)
                minDeltaR = _lepDeltaR[leptonCounter][j];
        }
        if (minDeltaR < 0.05) continue;
        
        _chargeCons[leptonCounter] = iE->isGsfCtfScPixChargeConsistent();
        
        _leptonIndex = i;
        if (leptonCounter == nLeptonsMax) continue;
        _flavors[leptonCounter] = 0;
        _charges[leptonCounter] = iE->charge();
        //_isolation = pfRelIso(iE, myRho);
        _isolation[leptonCounter][0] = pfRelIso(&*iE, myRhoJetsNC);
        _ipPV[leptonCounter] = iE->gsfTrack()->dxy(PV);
        _ipZPV[leptonCounter] = iE->gsfTrack()->dz(PV);
        
        _3dIP[leptonCounter]    = iE->dB(pat::Electron::PV3D);
        _3dIPerr[leptonCounter] = iE->edB(pat::Electron::PV3D);
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);
        
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        
        bool passed = false;
        int etak=0;
        //cut = std::min( looseMVA[etak][0], std::max( looseMVA[etak][1], looseMVA[etak][0] - (looseMVA[etak][0] - looseMVA[etak][1])/10.*(pt-15));

        if (TMath::Abs(iE->eta()) < 0.8 ) {
            etak = 0;
        } else if (TMath::Abs(iE->eta()) < 1.479 ) {
            etak = 1;
        } else {
            etak = 2;
        }
        passed = _mvaValue[leptonCounter] > std::min( looseMVA[etak][0], std::max( looseMVA[etak][1], looseMVA[etak][0] - (looseMVA[etak][0] - looseMVA[etak][1])/10.*(_lPt[leptonCounter]-15)));
        
        _istightID[leptonCounter] = passed;
        
        
        _miniisolation[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.2, 10., false, false, myRhoJetsNC);
        _miniisolation[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.3, 10., false, false, myRhoJetsNC);
        _miniisolation[leptonCounter][2] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.4, 10., false, false, myRhoJetsNC);
        
        _miniisolationCharged[leptonCounter][0] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.2, 10., false, true, myRhoJetsNC);
        _miniisolationCharged[leptonCounter][1] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.3, 10., false, true, myRhoJetsNC);
        _miniisolationCharged[leptonCounter][2] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.4, 10., false, true, myRhoJetsNC);
        
        fillCloseJetVars(leptonCounter);
        _lepMVA[leptonCounter] = lepMVAvalue(leptonCounter);
        
        
        _triggerMatch[leptonCounter] = triggerEmulator(&*iE);
        _triggerIsoMatch[leptonCounter] = isoTriggerEmulator(&*iE);
        
        _isloose[leptonCounter] =
        //triggerEmulator(&*iE) &&
        //&& (_miniisolation[leptonCounter][0] < _miniIsoCut)
        (_isolation[leptonCounter][0] < _relIsoCut || _isolation[leptonCounter][0]*_lPt[leptonCounter] < _absIsoCut)
        && (fabs(_3dIPsig[leptonCounter]) < _looseSIP)
        && fabs(_ipPV[leptonCounter]) < _looseD0E ;
        //&& _triggerMatch[leptonCounter];
        if (_isloose[leptonCounter])
            _istight[leptonCounter] = _istightID[leptonCounter] &&
            (_miniisolation[leptonCounter][0] < _miniisocut[_flavors[leptonCounter]]
             && (_ptratio[leptonCounter] > _ptratiocut[_flavors[leptonCounter]]
                 || _ptrel[leptonCounter] > _ptrelcut[_flavors[leptonCounter]]))
            && (fabs(_3dIPsig[leptonCounter]) < _looseSIP)
            ;
        else _istight[leptonCounter] = false;
        
        if (!_isloose[leptonCounter]) continue;
        //if (!_istightID[leptonCounter]) continue;
        //if (!_istight[leptonCounter]) continue;
        
        _istightIso[leptonCounter] = _istight[leptonCounter] && _istightID[leptonCounter]
        && (_miniisolation[leptonCounter][0] < _miniisocut[_flavors[leptonCounter]]
            && (_ptratio[leptonCounter] > _ptratiocut[_flavors[leptonCounter]]
                || _ptrel[leptonCounter] > _ptrelcut[_flavors[leptonCounter]]));
        
        //if (!_triggerMatch[leptonCounter]) continue;
        //if (!_istightID[leptonCounter]) continue;
        //if (!_istightIso[leptonCounter]) continue;
        
        
        if (Sample=="ElectronsMC") {
            //**************************************************************************************
            // MC
            //**************************************************************************************
            //std::cout<<_eventNb<<std::endl;
            const GenParticle* mc = GPM.matchedMC(&*iE, 11);
            _originPhot[leptonCounter] = -2;
            if ( mc!=0 ) {
                fillMCVars(mc, leptonCounter);
                //Vertex::Point PVmc = mcMom->vertex();
                _ipPVmc[leptonCounter] = TMath::Abs(iE->gsfTrack()->dxy(PVmc));
                
                //if (_origin[leptonCounter]==0 && (_charges[leptonCounter]!=_chargesMC[leptonCounter])) {
                /*std::cout<<"charge "<<_charges[leptonCounter]<<"; origin "<<_origin[leptonCounter]<<"; originReduced "<<
                 _originDetailed[leptonCounter]<<"; "<<
                 mc->isPromptFinalState()<<" "<<mc->fromHardProcessFinalState()
                 <<std::endl;
                 std::cout<<_lPt[leptonCounter]<<" "<<mc->pt()<<std::endl;
                 GPM.printInheritance(mc);*/
                 //}
            }
            else {
                mc = GPM.matchedMC(&*iE);
                if ( mc!=0 ) {
                    fillMCVars(mc, leptonCounter);
                    if (mc->pdgId() != 22) //not a photon
                        _origin[leptonCounter] = GPM.originReduced(_originDetailed[leptonCounter]);
                    else {
                        _originPhot[leptonCounter] = photonOrigin(mc);
                    }
                    /*std::cout<<"did not find normal: photon orig "<<_originPhot[leptonCounter]<<std::endl;
                     std::cout<<"charge "<<_charges[leptonCounter]<<"; origin "<<_origin[leptonCounter]<<"; originReduced "<<
                     _originDetailed[leptonCounter]<<"; "<<
                     mc->isPromptFinalState()<<" "<<mc->fromHardProcessFinalState()
                     <<std::endl;
                    std::cout<<"photon match "<<_originPhot[leptonCounter]<<std::endl;
                     std::cout<<_lPt[leptonCounter]<<" "<<mc->pt()<<std::endl;
                     GPM.printInheritance(mc);*/
                    /*
                     did not find normal:
                     charge -1; origin 3; originReduced 17; 1 1
                     30.4368 31.791
                     gamma (1)     <--  gamma (23)     <--   MANY 21(0), -4(0), )
                     */
                    
                } else {
                    //std::cout<<"No match e"<<std::endl;
                    _originDetailed[leptonCounter] = -1;
                    _origin[leptonCounter] = 4;
                    _mompt[leptonCounter] = 0;
                    _momphi[leptonCounter] = 0;
                    _mometa[leptonCounter] = 0;
                    _mompdg[leptonCounter] = 0;
                }

            }
            
            //if (_origin[leptonCounter]!=0) continue;
            
            fillIsoMCVars(leptonCounter);
            if (_closeIndex[leptonCounter] >=0)
                matchCloseJet(leptonCounter);
            
            //**************************************************************************************
            // MC *
            //**************************************************************************************
        }
        
        
        
        _mllZ[leptonCounter] = 9999.;
        _mllG[leptonCounter] = 9999.;
        _mllZj[leptonCounter] = 9999.;
        _mllGj[leptonCounter] = 9999.;
        _lepDeltaRloose[i] = 9999.;

        leptonCounter++;     
    }
    //**************************************************
    //*******zhud: Finished Analyzing Reco Electrons****
    //**************************************************

    _nEle = leptonCounter - _nMu;
    
    
    
    _nTau = leptonCounter - _nMu - _nEle;
    
    _nLeptons = leptonCounter;
    //if (leptonCounter < 3) return;     //zhud: for efficiency calculations, comment this line out to also have events with < 3
    
    _n_Jets = 0;
    _n_bJets = 0;
    _n_Jets30 = 0;
    _HT = 0;
    TLorentzVector jt;
    
    for(unsigned int i = 0 ; i < SelectedJets.size() ;i++ ){
        if (fabs(SelectedJets[i]->eta()) > 2.4) continue;
        
        double uncPt = (SelectedJets[i]->correctedP4("Uncorrected")).Pt();
        double uncEta = (SelectedJets[i]->correctedP4("Uncorrected")).Eta();
        //double uncPhi = (SelectedJets[i]->correctedP4("Uncorrected")).Phi();
        
        double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJets[i]->jetArea(),_corrLevel);
        
        //std::cout<<SelectedJets[i]->pt()<<" "<<uncPt*corr<<" "<<uncPt<<" "<<uncEta<<std::endl;
        //std::cout<<myRhoJets<<" "<<fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJets[i]->jetArea(),"L1FastJet")<<" "<<
        //fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJets[i]->jetArea(),"L2L3Residual")<<" "<<
        //fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJets[i]->jetArea(),"L3Absolute")<<std::endl;
        
        //if (SelectedJets[i]->pt() < _jetPtCut) continue;
        if (uncPt*corr < _jetPtCut) continue;
        
        _jetEta[_n_Jets] = SelectedJets[i]->eta();
        _jetPhi[_n_Jets] = SelectedJets[i]->phi();
        //_jetPt[_n_Jets] = SelectedJets[i]->pt();
        _jetPt[_n_Jets] = uncPt*corr;
        _jetE[_n_Jets] = (SelectedJets[i]->correctedP4("Uncorrected")).E()*corr;
        //_jetM[_n_Jets] = SelectedJets[i]->mass();
        
        
        jt.SetPtEtaPhiE(_jetPt[_n_Jets],_jetEta[_n_Jets],_jetPhi[_n_Jets],_jetE[_n_Jets]);
        _jetM[_n_Jets] = jt.M();
        
        _jetDeltaRloose[_n_Jets] = 99999;
        for (int j=0; j!=_nLeptons; ++j) {
            _jetDeltaR[_n_Jets][j] = ((TLorentzVector *)_leptonP4->At(j))->DeltaR( jt );
            if (_lPt[j] > 10 && ((j < _nMu+_nEle && _islooseID[j]) || (_isloose[j] && _flavors[j] == 2)) && _jetDeltaR[_n_Jets][j] < _jetDeltaRloose[_n_Jets]) {
                _jetDeltaRloose[_n_Jets] = _jetDeltaR[_n_Jets][j];
            }
        }
        
        _csv[_n_Jets] = SelectedJets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        
        if(_csv[_n_Jets] > 0.800) {
            _bTagged[_n_Jets] = true;
            _n_bJets++;
        } else _bTagged[_n_Jets] = false;
        
        if (_jetPt[_n_Jets] > _jetPtCut) {
            _HT+= _jetPt[_n_Jets];
            _n_Jets30++;
        }
        
        if (Sample=="ElectronsMC") {
            //matchAnyJet(_n_Jets);
            _jetpFlav[_n_Jets] = SelectedJets[i]->partonFlavour();
            _jethFlav[_n_Jets] = SelectedJets[i]->hadronFlavour();
        }
        
        _n_Jets++;
    }
    
    for (int i1=0; i1!=_nLeptons; ++i1) {
        _lepDeltaRloose[i1] = 9999;
        for (int j1=0; j1!=_nLeptons; ++j1) {
            _lepDeltaR[i1][j1] = ((TLorentzVector *)_leptonP4->At(j1))->DeltaR( *((TLorentzVector *)_leptonP4->At(i1)) );
            if (_flavors[i1]!=_flavors[j1] && j1 < _nMu+_nEle && _lepDeltaR[i1][j1] < _lepDeltaRloose[i1]) {
                _lepDeltaRloose[i1] = _lepDeltaR[i1][j1];
            }
            if ((i1!=j1) && (_flavors[i1]==_flavors[j1]) && (_charges[i1]!=_charges[j1])) {
                double mllLoc = Mll_calc( *((TLorentzVector*)_leptonP4->At(i1)), *((TLorentzVector*)_leptonP4->At(j1)));
                if (mllLoc < _mllG[i1])
                    _mllG[i1] = mllLoc;
            }
        }
    }
    
    
    for (int i=0; i!=_nLeptons; ++i) for (int j=0; j!=_nLeptons; ++j) {
        _mll[i][j] =  Mll_calc( *((TLorentzVector*)_leptonP4->At(i)), *((TLorentzVector*)_leptonP4->At(j)));
    }

    outputTree->Fill();
    
    
}

void HNL::fillMCVars(const GenParticle* mc, const int leptonCounter) {
    
    _lPtmc[leptonCounter] = mc->pt();
    _lEmc[leptonCounter] = mc->energy();
    _lPhimc[leptonCounter] = mc->phi();
    _lEtamc[leptonCounter] = mc->eta();
    _pdgmc[leptonCounter] = mc->pdgId();

    
    _chargesMC[leptonCounter] = mc->charge();
    
    _isPromptFinalState[leptonCounter] = mc->isPromptFinalState();
    _fromHardProcessFinalState[leptonCounter] = mc->fromHardProcessFinalState();
    
    _originDetailed[leptonCounter] = GPM.origin(mc);
    _origin[leptonCounter] = GPM.originReduced(_originDetailed[leptonCounter]);
    
    if (_isPromptFinalState[leptonCounter] || _fromHardProcessFinalState[leptonCounter]) {
        _origin[leptonCounter] = 0;
    }
    
    const GenParticle* mcMom = GPM.getMotherParton(mc);
    if (mcMom!=NULL) {
        
        //std::cout<<"Mother: "<<std::endl;
        //std::cout<<mcMom->pdgId()<<" "<<mcMom->pt()<<std::endl;
        //std::cout<<mc->pt()<<std::endl;
        
        _mompt[leptonCounter] = mcMom->pt();
        _momphi[leptonCounter] = mcMom->phi();
        _mometa[leptonCounter] = mcMom->eta();
        _mompdg[leptonCounter] = mcMom->pdgId();
        
        PVmc = mcMom->vertex();
        
        //std::cout<<"d0: "<<_ipPVmc<<" "<<_ipPV<<std::endl;
        
        TLorentzVector Gen0;
        Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
        
        
        _isolationMC[leptonCounter][0] = 0;
        _isolationMC[leptonCounter][1] = 0;
        std::vector<GenParticleRef> dauts;
        for (unsigned int dau=0; dau!=mcMom->numberOfDaughters(); ++dau) {
            GenParticleRef daut = mcMom->daughterRef(dau);
            dauts.push_back(daut);
        }
        unsigned int counterD = 0;
        //if (_isolation == 0)
        //    std::cout<<"==== "<<chargedHadronIso*iM->pt()<<" "<<neutralHadronIso*iM->pt()<<" "<<photonIso*iM->pt()<<" "<<beta*iM->pt()<<std::endl;
        while (counterD < dauts.size()) {
            if (dauts.at(counterD)->status() == 1) {
                _isolationMC[leptonCounter][0]+=dauts.at(counterD)->pt();
                if ((fabs(dauts.at(counterD)->pdgId())!= 12) && (fabs(dauts.at(counterD)->pdgId())!= 14) && (fabs(dauts.at(counterD)->pdgId())!= 16)) {
                    _isolationMC[leptonCounter][1]+=dauts.at(counterD)->pt();
                    /*if (_isolation == 0) {
                     TLorentzVector dauV;
                     dauV.SetPtEtaPhiE(dauts.at(counterD)->pt(), dauts.at(counterD)->eta(), dauts.at(counterD)->phi(), dauts.at(counterD)->energy());
                     double deltaR1 = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR(dauV);
                     std::cout<<counterD<<" "<<dauts.at(counterD)->pdgId()<<" "<<dauts.at(counterD)->pt()<<
                     " "<<dauts.at(counterD)->eta()<<" "<<dauts.at(counterD)->phi()<<" "<<dauts.at(counterD)->energy()<<" in "<<deltaR1<<std::endl;
                     }*/
                } else {
                    TLorentzVector Gen;
                    Gen.SetPtEtaPhiE( dauts.at(counterD)->pt(), dauts.at(counterD)->eta(), dauts.at(counterD)->phi(), dauts.at(counterD)->energy() );
                    Gen0 += Gen;
                }
            } else {
                for (unsigned int dau=0; dau!=dauts.at(counterD)->numberOfDaughters(); ++dau) {
                    GenParticleRef daut = dauts.at(counterD)->daughterRef(dau);
                    dauts.push_back(daut);
                }
            }
            counterD++;
        }
        
        //std::cout<<"PTs: "<<mc->pt()<<" "<<iM->pt()<<std::endl;
        _isolationMC[leptonCounter][0]/=mc->pt();
        _isolationMC[leptonCounter][0]-=1;
        _isolationMC[leptonCounter][1]/=mc->pt();
        _isolationMC[leptonCounter][1]-=1;
        
        _nuPtmc[leptonCounter] = Gen0.Pt();
        _nuPhimc[leptonCounter] = Gen0.Phi();
        _nuEtamc[leptonCounter] = Gen0.Eta();
        _nuEmc[leptonCounter] = Gen0.E();
        
        TLorentzVector lmc;
        lmc.SetPtEtaPhiE(_lPtmc[leptonCounter],_lEtamc[leptonCounter],_lPhimc[leptonCounter],_lEmc[leptonCounter]);
        _mtmc[leptonCounter] = MT_calc(lmc, _nuPtmc[leptonCounter], _nuPhimc[leptonCounter]);
    } else {
        _mompt[leptonCounter] = 0;
        _momphi[leptonCounter] = 0;
        _mometa[leptonCounter] = 0;
        _mompdg[leptonCounter] = 0;
    }
    
}

void HNL::fillCloseJetVars(const int leptonCounter) {
    _closeJetPtAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
    _closeJetAngAll[leptonCounter] = 10000;
    _ptRelAll[leptonCounter] = 0;
    _ptrel[leptonCounter] = 0;
    _ptrel2[leptonCounter] = 0;
    _ptratio[leptonCounter] = 1.;
    _trackSelectionMultiplicity[leptonCounter] = 0.;
    
    _closeIndex[leptonCounter] = 0;
    TLorentzVector pJet;
    for(unsigned int k = 0 ; k < SelectedJetsAll.size() ;k++ ){
        double uncorrPt = (SelectedJetsAll[k]->correctedP4("Uncorrected")).Pt();
        double corr = fMetCorrector->getJetCorrectionRawPt( uncorrPt, (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta(), myRhoJets, SelectedJetsAll[k]->jetArea(),"L1FastJet");
        double corr2 = fMetCorrector->getJetCorrectionRawPt( uncorrPt, (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta(), myRhoJets, SelectedJetsAll[k]->jetArea(),_corrLevel);
        pJet.SetPtEtaPhiE( corr*uncorrPt, SelectedJetsAll[k]->eta(), SelectedJetsAll[k]->phi(), corr*(SelectedJetsAll[k]->correctedP4("Uncorrected")).E());
        pJet-=*((TLorentzVector *)_leptonP4->At(leptonCounter));
        pJet*=corr2/corr;
        pJet+=*((TLorentzVector *)_leptonP4->At(leptonCounter));
        
        double ang = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pJet );
        
        //std::cout<<uncorrPt<<" "<<corr<<" "<<corr2<<" "<<corr2/corr<<std::endl;
        //std::cout<<pJet.Pt()<<std::endl;
        //std::cout<<pJet.Pt()<<std::endl;
        //std::cout<<pJet.Pt()<<std::endl;
        //std::cout<<pJet.Pt()<<std::endl;
        //std::cout<<ang<<std::endl;
        //std::cout<<std::endl;
        
        if (ang < _closeJetAngAll[leptonCounter]) {
            _closeJetAngAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pJet );
            _closeJetPtAll[leptonCounter] = pJet.Pt();
            _closeJetEtaAll[leptonCounter] = pJet.Eta();
            _closeJetPhiAll[leptonCounter] = pJet.Phi();
            _closeJetEAll[leptonCounter] = pJet.E();
            _closeJetMAll[leptonCounter] = pJet.M();
            _closeJetNconstAll[leptonCounter] = SelectedJetsAll[k]->numberOfDaughters();
            _closeJetCSVAll[leptonCounter] = SelectedJetsAll[k]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
            //std::cout<<_closeJetPtAll[leptonCounter]<<" "<<_closeJetEtaAll[leptonCounter]<<" "<<_closeJetPhiAll[leptonCounter]<<std::endl;
            
            //_closeJetPtAll[leptonCounter] = corr2/corr*(uncorrPt*corr - ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()) + ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
            //_closeJetEAll[leptonCounter] = corr2/corr*((SelectedJetsAll[k]->correctedP4("Uncorrected")).E()*corr - ((TLorentzVector *)_leptonP4->At(leptonCounter))->E()) + ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
            
            //std::cout<<uncorrPt<<" "<<corr<<" "<<corr2<<std::endl;
            //std::cout<<(uncorrPt*corr - ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt())*corr2 + ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()<<std::endl;
            //std::cout<<_closeJetPtAll[leptonCounter]<<std::endl;
            
            
            _ptRelAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect());
            _ptrel[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect() - ((TLorentzVector *)_leptonP4->At(leptonCounter))->Vect());
            //std::cout<<_ptrel[leptonCounter]<<std::endl;
            pJet-=*((TLorentzVector *)_leptonP4->At(leptonCounter));
            _ptrel2[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect());
            
            //_ptratio[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()/_closeJetPtAll[leptonCounter];
            _ptratio[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()/(_closeJetPtAll[leptonCounter]);
            //std::cout<<_ptrel[leptonCounter]<<std::endl;
            _closeIndex[leptonCounter] = k;
        }
    }
    
    if (_closeJetAngAll[leptonCounter] > 0.4) {
        _closeJetPtAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _closeJetPhiAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _closeJetEtaAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _closeJetNconstAll[leptonCounter] = 1;
        _closeJetCSVAll[leptonCounter] = 0;
        _closeJetAngAll[leptonCounter] = 0;
        _ptRelAll[leptonCounter] = 0;
        _ptrel[leptonCounter] = 0;
        _ptratio[leptonCounter] = 1.;
        _closeIndex[leptonCounter] = -1;
    } else {
        
        int k = _closeIndex[leptonCounter];
        
        int trackSelectionMult = 0;
        for(int unsigned it=0; it!=SelectedJetsAll[k]->numberOfDaughters(); ++it) {
            const pat::PackedCandidate * icand = dynamic_cast<const pat::PackedCandidate *> (SelectedJetsAll[k]->daughter(it));
            const reco::Track& trk = icand->pseudoTrack();
            Double_t deta = trk.eta() - SelectedJetsAll[k]->eta();
            Double_t dphi = TVector2::Phi_mpi_pi(trk.phi() - SelectedJetsAll[k]->phi());
            bool trackSelection =  trk.charge() != 0 && (TMath::Sqrt( deta*deta+dphi*dphi ) < 0.4) && (icand->fromPV() > 1) && (trk.hitPattern().numberOfValidHits()>=8) && (trk.hitPattern().numberOfValidPixelHits()>=2) && (trk.normalizedChi2()<5) && std::fabs(trk.dxy(PV))<0.2 && std::fabs(trk.dz(PV))<17 && trk.pt() > 1;
            if(trackSelection) trackSelectionMult++;
        }
        _trackSelectionMultiplicity[leptonCounter] = trackSelectionMult;
        
    }
    
    //std::cout<<_closeIndex[leptonCounter]<<" "<<_closeJetPtAll[leptonCounter]<<" "<<_closeJetAngAll[leptonCounter]<<" "<<
    //_ptrel[leptonCounter]<<" "<<_ptratio[leptonCounter]<<std::endl;
    for (int mi=0; mi!=3; ++mi) {
        if (_miniisolation[leptonCounter][0] < _multiConst[mi][0] && (_ptratio[leptonCounter] > _multiConst[mi][1] || _ptrel[leptonCounter] > _multiConst[mi][2]) )
            _multiisolation[leptonCounter][mi] = true;
        else
            _multiisolation[leptonCounter][mi] = false;
        
    }
    
}
void HNL::matchCloseJet(const int leptonCounter) {
    double minDeltaR3 = 9999;
    double minDeltaR2 = 9999;
    TLorentzVector Gen1, Gen2;
    const GenParticle* mc3a = 0;
    const GenParticle* mc2a = 0;
    Gen1.SetPtEtaPhiE(SelectedJetsAll.at(_closeIndex[leptonCounter])->pt(),SelectedJetsAll.at(_closeIndex[leptonCounter])->eta(),SelectedJetsAll.at(_closeIndex[leptonCounter])->phi(),SelectedJetsAll.at(_closeIndex[leptonCounter])->energy());
    
    for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
        int id = TMath::Abs(p->pdgId());
        
        if ((id > 0 && id < 6) || (id == 21) || (id == 22)) {
            if (p->status() != 2) {
                if (fabs(p->eta()) <= 10 )
                    Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
                else continue;
                //Gen2.SetPtEtaPhiM(p->pt(),0.00001,p->phi(),0);
                double deltaRcur = Gen1.DeltaR(Gen2);
                if (deltaRcur < minDeltaR3) {
                    mc3a = &*p;
                    minDeltaR3 = deltaRcur;
                }
                
            } else if (p->status() == 2) {
                if (fabs(p->eta()) <= 10 )
                    Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
                else continue;
                //Gen2.SetPtEtaPhiM(p->pt(),0.00001,p->phi(),0);
                double deltaRcur = Gen1.DeltaR(Gen2);
                if (deltaRcur < minDeltaR2) {
                    mc2a = &*p;
                    minDeltaR2 = deltaRcur;
                }
            }
        }
    }
    
    if ((minDeltaR3 < 0.5 || minDeltaR2 == 9999) && mc3a!=0) {
        _closeJetPtAllMC[leptonCounter] = mc3a->pt();
        _closeJetPtAllstatus[leptonCounter] = 3;
        _partonIdMatched[leptonCounter] = mc3a->pdgId();
    } else if (mc2a!=0) {
        _closeJetPtAllMC[leptonCounter] = mc2a->pt();
        _closeJetPtAllstatus[leptonCounter] = 2;
        _partonIdMatched[leptonCounter] = mc2a->pdgId();
    } else {
        _closeJetPtAllMC[leptonCounter] = 0;
        _closeJetPtAllstatus[leptonCounter] = 0;
        _partonIdMatched[leptonCounter] = 0;
    }
}



void HNL::matchAnyJet(const int index) {
    double minDeltaR3 = 9999;
    TLorentzVector Gen1, Gen2;
    const GenParticle* mc3a = 0;
    Gen1.SetPtEtaPhiE(_jetPt[index],_jetEta[index],_jetPhi[index],_jetE[index]);
    
    for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
        int id = TMath::Abs(p->pdgId());
        
        if ((id > 0 && id < 6) || (id == 21) || (id == 22)) {
            //if (p->status() != 2) {
            if (fabs(p->eta()) <= 3.0 )
                Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
            else continue;
            //Gen2.SetPtEtaPhiM(p->pt(),0.00001,p->phi(),0);
            double deltaRcur = Gen1.DeltaR(Gen2);
            if (deltaRcur < minDeltaR3 && fabs((p->pt()-_jetPt[index])/_jetPt[index]) < 0.5) {
                mc3a = &*p;
                minDeltaR3 = deltaRcur;
            }
            //}
        }
    }
    
    if ((minDeltaR3 < 0.5) && mc3a!=0) {
        _jetpFlav[index] = TMath::Abs(mc3a->pdgId());
        //std::cout<<"flav: "<<_jetFlav[index]<<"; Delta R="<<minDeltaR3<<" "<<(mc3a->pt()-_jetPt[index])/_jetPt[index]<<"; status = "<<mc3a->status()<<std::endl;
    } else {
        _jetpFlav[index] = -1;
    }
}

double HNL::matchCloseParticle(const int leptonCounter, const int pdgIDm) {
    double minDeltaR = 9999;
    TLorentzVector Gen2;
    
    for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
        int id = TMath::Abs(p->pdgId());
        
        if (id == pdgIDm) {
            if (fabs(p->eta()) <= 10 )
                Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
            else continue;
            //Gen2.SetPtEtaPhiM(p->pt(),0.00001,p->phi(),0);
            double deltaRcur = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( Gen2 );
            if (deltaRcur < minDeltaR) {
                minDeltaR = deltaRcur;
            }
            
        }
    }
    
    return minDeltaR;
}


void HNL::fillIsoMCVars(const int leptonCounter) {
    if( TheGenParticles.isValid() )
    {
        //if (_isolation[0] == 0) {
        //    std::cout<<"==="<<std::endl;
        //}
        _isolationMC[leptonCounter][2] = 0;
        _isolationMC[leptonCounter][3] = 0;
        for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
            //int id = TMath::Abs(p->pdgId());
            if (p->status() == 1) {
                TLorentzVector pmc; pmc.SetPtEtaPhiM( p->pt(), p->eta(), p->phi(), p->mass() );
                double ang = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pmc );
                if (ang < 0.3) {
                    _isolationMC[leptonCounter][2]+=p->pt();
                    if (fabs(p->pdgId())!= 12 && fabs(p->pdgId())!= 14 && fabs(p->pdgId())!= 16) {
                        _isolationMC[leptonCounter][3]+=p->pt();
                        //if (_isolation[0] == 0) {
                        //    std::cout<<" "<<p->pdgId()<<" "<<p->pt()<<
                        //    " "<<p->eta()<<" "<<p->phi()<<" "<<p->energy()<<" in "<<ang<<std::endl;
                        //}
                    }
                }
            }
        }
        double leptpt = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _isolationMC[leptonCounter][2]/=leptpt;
        _isolationMC[leptonCounter][3]/=leptpt;
        _isolationMC[leptonCounter][2]-=1;
        _isolationMC[leptonCounter][3]-=1;
        //std::cout<<"*************"<<std::endl;
        //for (int i=0; i!=4; ++i)
        //    std::cout<<_isolationMC[i]<<" ";
        //std::cout<<_isolation[0]<<" "<<_closeJetPtAllMC[0]/leptpt-1<<" "<<_closeJetPtAll[0]/leptpt-1<<std::endl;
    }
}

double HNL::lepMVAvalue(const int leptonCounter) {
    
    LepGood_pt = _lPt[leptonCounter];
    LepGood_eta = _lEta[leptonCounter];
    LepGood_jetNDauChargedMVASel = _trackSelectionMultiplicity[leptonCounter];
    LepGood_miniRelIsoCharged = _miniisolationCharged[leptonCounter][0];
    LepGood_miniRelIsoNeutral = _miniisolation[leptonCounter][0] - _miniisolationCharged[leptonCounter][0];
    LepGood_jetPtRelv2 = _ptrel[leptonCounter];
    LepGood_jetPtRatio = TMath::Min(_ptratio[leptonCounter],1.5);
    LepGood_jetBTagCSV = TMath::Max(_closeJetCSVAll[leptonCounter],0.);
    LepGood_sip3d = _3dIPsig[leptonCounter];
    LepGood_dxy = TMath::Log(fabs(_ipPV[leptonCounter]));
    LepGood_dz = TMath::Log(fabs(_ipZPV[leptonCounter]));
    LepGood_mvaIdSpring15 = _mvaValue[leptonCounter];
    LepGood_segmentCompatibility = _muonSegmentComp[leptonCounter];
    
    return reader[_flavors[leptonCounter]]->EvaluateMVA( "BDTG method" );  
}


void HNL::fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu) {
    
    hJet_ptRaw = (jet->correctedP4("Uncorrected")).Pt();
    hJet_genPt = genpt;
    hJet_pt = jet->pt();
    hJet_phi = jet->phi();
    hJet_eta = jet->eta();
    hJet_e = jet->energy();
    
    hJet_ptLeadTrack = 0;
    
    const reco::TrackRefVector &tracks =  jet->associatedTracks();
    for (unsigned int k=0; k!= tracks.size(); ++k) {
        if(tracks[k]->pt() > hJet_ptLeadTrack)
            hJet_ptLeadTrack = tracks[k]->pt();
    }
    
    hJet_vtx3dL = 0;
    hJet_vtx3deL = 0;
    hJet_vtxMass = 0;
    hJet_vtxPt = 0;
    
    const reco::SecondaryVertexTagInfo* scdVtx = jet->tagInfoSecondaryVertex("secondaryVertex");
    
    if (scdVtx) {
        //std::cout<<"Vertetx info: "<<scdVtx->nVertices()<<std::endl;
        if (scdVtx->nVertices()) {
            const reco::Vertex &sv1 = scdVtx->secondaryVertex(0);
            if (!sv1.isFake()) {
                Measurement1D distance1 = scdVtx->flightDistance(0, true);
                hJet_vtx3dL = distance1.value();
                hJet_vtx3deL = distance1.error();
                
                math::XYZTLorentzVectorD p4vtx = sv1.p4();
                hJet_vtxMass = p4vtx.M();
                hJet_vtxPt = p4vtx.Pt();
            }
        }
    }
    
    
    hJet_cef = jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction();
    hJet_nconstituents = jet->getPFConstituents().size();
    hJet_JECUnc = fMetCorrector->getJECUncertainty(jet->pt(),jet->eta());
    
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    
    hJet_SoftLeptptRel = pLep.Perp(pJet.Vect());
    hJet_SoftLeptPt = mu->pt();
    hJet_SoftLeptdR = pLep.DeltaR(pJet);
    
    hJet_SoftLeptIdlooseMu = 1;
    hJet_SoftLeptId95 = 1;
}

void HNL::fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* mu) {
    
    hJet_ptRaw = (jet->correctedP4("Uncorrected")).Pt();
    hJet_genPt = genpt;
    hJet_pt = jet->pt();
    hJet_phi = jet->phi();
    hJet_eta = jet->eta();
    hJet_e = jet->energy();
    
    hJet_ptLeadTrack = 0;
    
    const reco::TrackRefVector &tracks =  jet->associatedTracks();
    for (unsigned int k=0; k!= tracks.size(); ++k) {
        if(tracks[k]->pt() > hJet_ptLeadTrack)
            hJet_ptLeadTrack = tracks[k]->pt();
    }
    
    hJet_vtx3dL = 0;
    hJet_vtx3deL = 0;
    hJet_vtxMass = 0;
    hJet_vtxPt = 0;
    
    const reco::SecondaryVertexTagInfo* scdVtx = jet->tagInfoSecondaryVertex("secondaryVertex");
    
    if (scdVtx) {
        //std::cout<<"Vertetx info: "<<scdVtx->nVertices()<<std::endl;
        if (scdVtx->nVertices()) {
            const reco::Vertex &sv1 = scdVtx->secondaryVertex(0);
            if (!sv1.isFake()) {
                Measurement1D distance1 = scdVtx->flightDistance(0, true);
                hJet_vtx3dL = distance1.value();
                hJet_vtx3deL = distance1.error();
                
                math::XYZTLorentzVectorD p4vtx = sv1.p4();
                hJet_vtxMass = p4vtx.M();
                hJet_vtxPt = p4vtx.Pt();
            }
        }
    }
    
    
    hJet_cef = jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction();
    hJet_nconstituents = jet->getPFConstituents().size();
    hJet_JECUnc = fMetCorrector->getJECUncertainty(jet->pt(),jet->eta());
    
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    
    hJet_SoftLeptptRel = pLep.Perp(pJet.Vect());
    hJet_SoftLeptPt = mu->pt();
    hJet_SoftLeptdR = pLep.DeltaR(pJet);
    
    hJet_SoftLeptIdlooseMu = 1;
    hJet_SoftLeptId95 = 1;
}

int HNL::photonOrigin(const GenParticle* photon) {
    //implementation of a flag for the photon truth matching:
    //https://hypernews.cern.ch/HyperNews/CMS/get/susy-interpretations/192.html
    //-1: not a photon
    //0: direct prompt photons (prompt and delta R > 0.4)
    //1: fragmentation photons (prompt and delta R < 0.4)
    //2: non-prompt photons
    
    if (photon->pdgId() != 22) return -1;
    if (!photon->isPromptFinalState()) return 2;
    if (photon->pt() < 10) return 1;
    
    TLorentzVector photonP4; photonP4.SetPtEtaPhiE(photon->pt(),photon->eta(),photon->phi(),photon->energy());
    TLorentzVector partonP4;
    bool smallDR = false;
    for( unsigned p=0; p<TheGenParticles->size(); ++p ){
        if ((*TheGenParticles)[p].status() != 23) continue;
        if (!((*TheGenParticles)[p].pdgId() == 21 || (fabs((*TheGenParticles)[p].pdgId()) > 0 && fabs((*TheGenParticles)[p].pdgId()) <7))) continue;
        partonP4.SetPtEtaPhiE((*TheGenParticles)[p].pt(),(*TheGenParticles)[p].eta(),(*TheGenParticles)[p].phi(),(*TheGenParticles)[p].energy());
        double deltaR = photonP4.DeltaR(partonP4);
        smallDR |= (deltaR < 0.05);
    }
    if (smallDR) return 1;
    else return 0;
}

void HNL::bookTree() {
    
    //outputTree->Branch("_leptonP4", "TClonesArray", &_leptonP4, 32000, 0);
    //outputTree->Branch("_jetP4", "TClonesArray", &_jetP4, 32000, 0);
    outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
    outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
    outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");

    outputTree->Branch("_nGenLep", &_nGenLep, "_nGenLep/I");
    outputTree->Branch("_nGenE", &_nGenE, "_nGenE/I");
    outputTree->Branch("_nGenMu", &_nGenMu, "_nGenMu/I");
    outputTree->Branch("_nGenStatusNot1", &_nGenStatusNot1, "_nGenStatusNot1/I");
    outputTree->Branch("_nGenHNL",&_nGenHNL,"_nGenHNL/I");
    
    outputTree->Branch("mChi20", &mChi20, "mChi20/D");
    outputTree->Branch("mChi10", &mChi10, "mChi10/D");
    
    outputTree->Branch("_weight", &_weight, "_weight/D");
    outputTree->Branch("_genHT", &_genHT, "_genHT/D");
    outputTree->Branch("_nLeptons", &_nLeptons, "_nLeptons/I");
    

    //zhud: For the Efficiency measurement
    outputTree->Branch("_nGenHNLMu", &_nGenHNLMu,"_nGenHNLMu/I");
    outputTree->Branch("_nlHNLMu", &_nlHNLMu,"_nlHNLMu/I");
    
    outputTree->Branch("_GenHNLP",    &_GenHNLP,    "_GenHNLP[_nGenLep]/D");
    outputTree->Branch("_GenHNLEta",   &_GenHNLEta,   "_GenHNLEta[_nGenLep]/D");
    outputTree->Branch("_GenHNLPhi",   &_GenHNLPhi,   "_GenHNLPhi[_nGenLep]/D");
    outputTree->Branch("_GenHNLE",     &_GenHNLE,     "_GenHNLE[_nGenLep]/D");
    outputTree->Branch("_GenHNLGamma", &_GenHNLGamma, "_GenHNLGamma[_nGenLep]/D");
    outputTree->Branch("_GenHNLBeta ", &_GenHNLBeta , "_GenHNLBeta[_nGenLep]/D");


    outputTree->Branch("_GenHNLMuPt", &_GenHNLMuPt,"_GenHNLMuPt[_nGenLep]/D");
    outputTree->Branch("_GenHNLMuEta", &_GenHNLMuEta,"_GenHNLMuEta[_nGenLep]/D");
    outputTree->Branch("_GenHNLMuPhi", &_GenHNLMuPhi,"_GenHNLMuPhi[_nGenLep]/D");
    outputTree->Branch("_GenHNLMuE", &_GenHNLMuE,"_GenHNLMuE[_nGenLep]/D");   

    outputTree->Branch("_lHNLMuPtmc", &_lHNLMuPtmc, "_lHNLMuPtmc[_nLeptons]/D");
    
    outputTree->Branch("_lMuPt", &_lMuPt, "_lMuPt[_nLeptons]/D");
    outputTree->Branch("_lMuPhi", &_lMuPhi, "_lMuPhi[_nLeptons]/D");
    outputTree->Branch("_lMuEta", &_lMuEta, "_lMuEta[_nLeptons]/D");
    outputTree->Branch("_lMuE", &_lMuE, "_lMuE[_nLeptons]/D");

    outputTree->Branch("_GenMuPt", &_GenMuPt, "_GenMuPt[_nGenLep]/D");
    outputTree->Branch("_GenMuPhi", &_GenMuPhi, "_GenMuPhi[_nGenLep]/D");
    outputTree->Branch("_GenMuEta", &_GenMuEta, "_GenMuEta[_nGenLep]/D");
    outputTree->Branch("_GenMuE", &_GenMuE, "_GenMuE[_nGenLep]/D");

    outputTree->Branch("_GenMuVx", &_GenMuVx,"_GenMuVx[_nGenLep]/D"); 
    outputTree->Branch("_GenMuVy", &_GenMuVy,"_GenMuVy[_nGenLep]/D");
    outputTree->Branch("_GenMuVxy", &_GenMuVxy,"_GenMuVxy[_nGenLep]/D");
    outputTree->Branch("_GenMuVz", &_GenMuVz,"_GenMuVz[_nGenLep]/D");
    outputTree->Branch("_GenMuVxyz", &_GenMuVxyz,"_GenMuVxyz[_nGenLep]/D");

    outputTree->Branch("_GenHNLVxProd", &_GenHNLVxProd,"_GenHNLVxProd/D"); 
    outputTree->Branch("_GenHNLVyProd", &_GenHNLVyProd,"_GenHNLVyProd/D");
    outputTree->Branch("_GenHNLVxyProd", &_GenHNLVxyProd,"_GenHNLVxyProd/D");
    outputTree->Branch("_GenHNLVzProd", &_GenHNLVzProd,"_GenHNLVzProd/D");
    outputTree->Branch("_GenHNLVxyzProd", &_GenHNLVxyzProd,"_GenHNLVxyzProd/D");

    outputTree->Branch("_GenHNLVxDecay", &_GenHNLVxDecay,"_GenHNLVxDecay/D"); 
    outputTree->Branch("_GenHNLVyDecay", &_GenHNLVyDecay,"_GenHNLVyDecay/D");
    outputTree->Branch("_GenHNLVxyDecay", &_GenHNLVxyDecay,"_GenHNLVxyDecay/D");
    outputTree->Branch("_GenHNLVzDecay", &_GenHNLVzDecay,"_GenHNLVzDecay/D");
    outputTree->Branch("_GenHNLVxyzDecay", &_GenHNLVxyzDecay,"_GenHNLVxyzDecay/D");
    outputTree->Branch("_GenHNLVxyzDecay", &_GenHNLVxyzDecay,"_GenHNLVxyzDecay/D");
    outputTree->Branch("_GenHNLMu3DDisplacement", &_GenHNLMu3DDisplacement,"_GenHNLMu3DDisplacement/D");

    outputTree->Branch("_foundGenMuPt", &_foundGenMuPt, "_foundGenMuPt[_nGenLep]/D");
    outputTree->Branch("_isdetectedMu", &_isdetectedMu, "_isdetectedMu[_nGenLep]/O");
    outputTree->Branch("_isfromHNLMu", &_isfromHNLMu, "_isfromHNLMu[_nGenLep]/O");
    outputTree->Branch("_ispromptMu", &_ispromptMu, "_ispromptMu[_nGenLep]/O");
    outputTree->Branch("_isdetectedDSAMu", &_isdetectedDSAMu, "_isdetectedDSAMu[_nGenLep]/O");


    //zhud: end efficiency measurement KPIs


    outputTree->Branch("_nDSAMu", &_nDSAMu, "_nDSAMu/I");

    outputTree->Branch("_nEle", &_nEle, "_nEle/I");
    outputTree->Branch("_nMu", &_nMu, "_nMu/I");
    outputTree->Branch("_nTau", &_nTau, "_nTau/I");
    
    outputTree->Branch("_lPt", &_lPt, "_lPt[_nLeptons]/D");
    outputTree->Branch("_lEta", &_lEta, "_lEta[_nLeptons]/D");
    outputTree->Branch("_lPhi", &_lPhi, "_lPhi[_nLeptons]/D");
    outputTree->Branch("_lE", &_lE, "_lE[_nLeptons]/D");
    
    outputTree->Branch("_flavors", &_flavors, "_flavors[_nLeptons]/I");
    outputTree->Branch("_charges", &_charges, "_charges[_nLeptons]/I");
    outputTree->Branch("_chargesMC", &_chargesMC, "_chargesMC[_nLeptons]/I");

    outputTree->Branch("_chargeCons", &_chargeCons, "_chargeCons[_nLeptons]/O");


    
    
    outputTree->Branch("_l1HTT", &_l1HTT, "_l1HTT/D");
    outputTree->Branch("_genPhotOrigin", &_genPhotOrigin, "_genPhotOrigin/I");
    

    
    outputTree->Branch("_ptSystem", &_ptSystem, "_ptSystem/D");
    outputTree->Branch("_ptTTSystem", &_ptTTSystem, "_ptTTSystem/D");
    
    outputTree->Branch("_nLooseLeptons", &_nLooseLeptons, "_nLooseLeptons/I");
    
    outputTree->Branch("_pdgmc", &_pdgmc, "_pdgmc[_nLeptons]/I");

    
    
    
     outputTree->Branch("_lPtmc", &_lPtmc, "_lPtmc[_nLeptons]/D");
     outputTree->Branch("_lEtamc", &_lEtamc, "_lEtamc[_nLeptons]/D");
     outputTree->Branch("_lPhimc", &_lPhimc, "_lPhimc[_nLeptons]/D");
     outputTree->Branch("_lEmc", &_lEmc, "_lEmc[_nLeptons]/D");
     
     outputTree->Branch("_nuPtmc", &_nuPtmc, "_nuPtmc[_nLeptons]/D");
     outputTree->Branch("_nuEtamc", &_nuEtamc, "_nuEtamc[_nLeptons]/D");
     outputTree->Branch("_nuPhimc", &_nuPhimc, "_nuPhimc[_nLeptons]/D");
     outputTree->Branch("_nuEmc", &_nuEmc, "_nuEmc[_nLeptons]/D");
     
     outputTree->Branch("_mtmc", &_mtmc, "_mtmc[_nLeptons]/D");



    
    // outputTree->Branch("_triggers2l", &_triggers2l, "_triggers2l[2][6]/O");
    // outputTree->Branch("_triggers2lbkp", &_triggers2lbkp, "_triggers2lbkp[2][6]/O");
    // outputTree->Branch("_triggersCS", &_triggersCS, "_triggersCS[5][5]/O");
    // outputTree->Branch("_triggersCSb", &_triggersCSb, "_triggersCSb[2]/O");
    // outputTree->Branch("_triggers1l", &_triggers1l, "_triggers1l[8]/O");
    // //outputTree->Branch("_triggersTau", &_triggersTau, "_triggersTau[6]/O");
    // outputTree->Branch("_triggersPost", &_triggersPost, "_triggersPost[20]/O");

    
    // outputTree->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
    // outputTree->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
    // outputTree->Branch("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, "Flag_HBHENoiseIsoFilter/O");
    // outputTree->Branch("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, "Flag_globalTightHalo2016Filter/O");
    // outputTree->Branch("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, "Flag_CSCTightHalo2015Filter/O");
    // outputTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
    // outputTree->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
    // outputTree->Branch("filterbadChCandidate", &filterbadChCandidate, "filterbadChCandidate/O");
    // outputTree->Branch("filterbadPFMuon", &filterbadPFMuon, "filterbadPFMuon/O");
    
    outputTree->Branch("_isolation", &_isolation, "_isolation[_nLeptons][2]/D");
    outputTree->Branch("_miniisolation", &_miniisolation, "_miniisolation[_nLeptons][3]/D");
    outputTree->Branch("_multiisolation", &_multiisolation, "_multiisolation[_nLeptons][3]/O");
    outputTree->Branch("_mvaValue", &_mvaValue, "_mvaValue[_nLeptons]/D");
    outputTree->Branch("_lepMVA", &_lepMVA, "_lepMVA[_nLeptons]/D");
    outputTree->Branch("_trackSelectionMultiplicity", &_trackSelectionMultiplicity, "_trackSelectionMultiplicity[_nLeptons]/I");
    
    
    
    /*
     outputTree->Branch("_isolationMC", &_isolationMC, "_isolationMC[_nLeptons][4]/D");
     
     outputTree->Branch("_index1", &_index1, "_index1/I");
     outputTree->Branch("_index2", &_index2, "_index2/I");
     
     outputTree->Branch("_sb", &_sb, "_sb/O");
     outputTree->Branch("_doubleF", &_doubleF, "_doubleF/O");*/
    
    outputTree->Branch("_originDetailed", &_originDetailed, "_originDetailed[_nLeptons]/I");
    outputTree->Branch("_origin", &_origin, "_origin[_nLeptons]/I");
    outputTree->Branch("_originPhot", &_originPhot, "_originPhot[_nLeptons]/I");
    
    
    
    
    outputTree->Branch("_isPromptFinalState", &_isPromptFinalState, "_isPromptFinalState[_nLeptons]/O");
    outputTree->Branch("_fromHardProcessFinalState", &_fromHardProcessFinalState, "_fromHardProcessFinalState[_nLeptons]/O");
    
    outputTree->Branch("_PVchi2", &_PVchi2, "_PVchi2/D");
    outputTree->Branch("_PVerr", &_PVerr, "_PVerr[3]/D");
    
    outputTree->Branch("_ipPV", &_ipPV, "_ipPV[_nLeptons]/D");
    outputTree->Branch("_ipPVerr", &_ipPVerr, "_ipPVerr[_nLeptons]/D");
    outputTree->Branch("_ipZPV", &_ipZPV, "_ipZPV[_nLeptons]/D");
    outputTree->Branch("_ipZPVerr", &_ipZPVerr, "_ipZPVerr[_nLeptons]/D");
    
    outputTree->Branch("_ipPVmc", &_ipPVmc, "_ipPVmc[_nLeptons]/D");
    
    outputTree->Branch("_3dIP", &_3dIP, "_3dIP[_nLeptons]/D");
    outputTree->Branch("_3dIPerr", &_3dIPerr, "_3dIPerr[_nLeptons]/D");
    outputTree->Branch("_3dIPsig", &_3dIPsig, "_3dIPsig[_nLeptons]/D");

    outputTree->Branch("_missingHits", &_missingHits, "_missingHits[_nLeptons]/I");

    

    outputTree->Branch("_dptoverpt", &_dptoverpt, "_dptoverpt[_nLeptons]/D");

    
    
    // outputTree->Branch("_mt", &_mt, "_mt[_nLeptons]/D");
    // outputTree->Branch("_mllZ", &_mllZ, "_mllZ[_nLeptons]/D");
    // outputTree->Branch("_mllG", &_mllG, "_mllG[_nLeptons]/D");
    // outputTree->Branch("_mllZj", &_mllZj, "_mllZj[_nLeptons]/D");
    // outputTree->Branch("_mllGj", &_mllGj, "_mllGj[_nLeptons]/D");
    // outputTree->Branch("_mll", &_mll, "_mll[_nLeptons][10]/D");
    // outputTree->Branch("_lepDeltaRloose", &_lepDeltaRloose, "_lepDeltaRloose[_nLeptons]/D");
    // outputTree->Branch("_lepDeltaR", &_lepDeltaR, "_lepDeltaR[_nLeptons][10]/D");
    
    
    // outputTree->Branch("_tauIDold", &_tauIDold, "_tauIDold[_nLeptons][5]/O");
    // outputTree->Branch("_tauIDnew", &_tauIDnew, "_tauIDnew[_nLeptons][5]/O");
    // outputTree->Branch("_tauID03", &_tauID03, "_tauID03[_nLeptons][4]/O");

    // outputTree->Branch("_tauIDeleVeto", &_tauIDeleVeto, "_tauIDeleVeto[_nLeptons][5]/O");
    // outputTree->Branch("_tauIDmuVeto", &_tauIDmuVeto, "_tauIDmuVeto[_nLeptons][2]/O");
    // outputTree->Branch("_tauDeltaR", &_tauDeltaR, "_tauDeltaR[_nLeptons]/D");

    
    // outputTree->Branch("_isloose", &_isloose, "_isloose[_nLeptons]/O");
    // outputTree->Branch("_istight", &_istight, "_istight[_nLeptons]/O");
    // outputTree->Branch("_istightID", &_istightID, "_istightID[_nLeptons]/O");
    // outputTree->Branch("_istightIso", &_istightIso, "_istightIso[_nLeptons]/O");
    // outputTree->Branch("_triggerMatch", &_triggerMatch, "_triggerMatch[_nLeptons]/O");
    // outputTree->Branch("_triggerIsoMatch", &_triggerIsoMatch, "_triggerIsoMatch[_nLeptons]/O");
    
    outputTree->Branch("_ptrel", &_ptrel, "_ptrel[_nLeptons]/D");
    //outputTree->Branch("_ptrel2", &_ptrel2, "_ptrel2[_nLeptons]/D");
    outputTree->Branch("_ptratio", &_ptratio, "_ptratio[_nLeptons]/D");
    
    // outputTree->Branch("_closeJetPtAll", &_closeJetPtAll, "_closeJetPtAll[_nLeptons]/D");
    // outputTree->Branch("_closeJetEtaAll", &_closeJetEtaAll, "_closeJetEtaAll[_nLeptons]/D");
    // outputTree->Branch("_closeJetPhiAll", &_closeJetPhiAll, "_closeJetPhiAll[_nLeptons]/D");
    // outputTree->Branch("_closeJetCSVAll", &_closeJetCSVAll, "_closeJetCSVAll[_nLeptons]/D");
    // outputTree->Branch("_closeJetNconstAll", &_closeJetNconstAll, "_closeJetNconstAll[_nLeptons]/I");
    // outputTree->Branch("_closeJetAngAll", &_closeJetAngAll, "_closeJetAngAll[_nLeptons]/D");
    // outputTree->Branch("_ptRelAll", &_ptRelAll, "_ptRelAll[_nLeptons]/D");
    
    
    // outputTree->Branch("_closeJetPtAllMC", &_closeJetPtAllMC, "_closeJetPtAllMC[_nLeptons]/D");
    // outputTree->Branch("_closeJetPtAllstatus", &_closeJetPtAllstatus, "_closeJetPtAllstatus[_nLeptons]/D");
    // outputTree->Branch("_partonIdMatched", &_partonIdMatched, "_partonIdMatched[_nLeptons]/I");
    // outputTree->Branch("_sameParton", &_sameParton, "_sameParton[_nLeptons]/O");
    
    
    // outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");
    // outputTree->Branch("_n_Interactions", &_n_Interactions, "_n_Interactions/I");
    // outputTree->Branch("_n_trueInteractions", &_n_trueInteractions, "_n_trueInteractions/D");
    
    
    // outputTree->Branch("_met", &_met, "_met/D");
    // outputTree->Branch("_met_phi", &_met_phi, "_met_phi/D");
    // outputTree->Branch("_HT", &_HT, "_HT/D");
    
    // outputTree->Branch("_genmet", &_genmet, "_genmet/D");
    // outputTree->Branch("_genmet_phi", &_genmet_phi, "_genmet_phi/D");
    
    // outputTree->Branch("_mompt", &_mompt, "_mompt[_nLeptons]/D");
    // outputTree->Branch("_momphi", &_momphi, "_momphi[_nLeptons]/D");
    // outputTree->Branch("_mometa", &_mometa, "_mometa[_nLeptons]/D");
    // outputTree->Branch("_mompdg", &_mompdg, "_mompdg[_nLeptons]/I");
    
    // outputTree->Branch("_n_bJets", &_n_bJets, "_n_bJets/I");
    // outputTree->Branch("_n_Jets", &_n_Jets, "_n_Jets/I");
    // outputTree->Branch("_n_Jets30", &_n_Jets30, "_n_Jets30/I");
    // outputTree->Branch("_bTagged", &_bTagged, "_bTagged[_n_Jets]/O");
    // outputTree->Branch("_jetEta", &_jetEta, "_jetEta[_n_Jets]/D");
    // outputTree->Branch("_jetPhi", &_jetPhi, "_jetPhi[_n_Jets]/D");
    // outputTree->Branch("_jetPt", &_jetPt, "_jetPt[_n_Jets]/D");
    // outputTree->Branch("_jetM", &_jetM, "_jetM[_n_Jets]/D");
    // outputTree->Branch("_jetE", &_jetE, "_jetE[_n_Jets]/D");
    // outputTree->Branch("_jethFlav", &_jethFlav, "_jethFlav[_n_Jets]/I");
    // outputTree->Branch("_jetpFlav", &_jetpFlav, "_jetpFlav[_n_Jets]/I");
    
    
    // outputTree->Branch("_csv", &_csv, "_csv[_n_Jets]/D");
    // outputTree->Branch("_jetDeltaR", &_jetDeltaR, "_jetDeltaR[_n_Jets][10]/D");
    // outputTree->Branch("_jetDeltaRloose", &_jetDeltaRloose, "_jetDeltaRloose[_n_Jets]/D");
    
    
    
}
DEFINE_FWK_MODULE(HNL);
