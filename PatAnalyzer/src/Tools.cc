#include "SUSYAnalyzer/PatAnalyzer/interface/Tools.h"

double tools::getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                      const reco::Candidate* ptcl,
                      double r_iso_min, double r_iso_max, double kt_scale,
                      bool use_pfweight, bool charged_only, double rho) {
    
    if (ptcl->pt()<5.) return 99999.;
    
    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
        if (fabs(dynamic_cast<const pat::Electron *>(ptcl)->superCluster()->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;} //0.08
    } else if(ptcl->isMuon()) {
        deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
    } else {
        //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }
    
    double iso_nh(0.); double iso_ch(0.);
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
        if (abs(pfc.pdgId())<7) continue;
        
        double dr = deltaR(pfc, *ptcl);
        if (dr > r_iso) continue;
        
        //////////////////  NEUTRALS  /////////////////////////
        if (pfc.charge()==0){
            if (pfc.pt()>ptThresh) {
                double wpf(1.);
                if (use_pfweight){
                    double wpv(0.), wpu(0.);
                    for (const pat::PackedCandidate &jpfc : *pfcands) {
                        double jdr = deltaR(pfc, jpfc);
                        if (pfc.charge()!=0 || jdr<0.00001) continue;
                        double jpt = jpfc.pt();
                        if (pfc.fromPV()>1) wpv *= jpt/jdr;
                        else wpu *= jpt/jdr;
                    }
                    wpv = log(wpv);
                    wpu = log(wpu);
                    wpf = wpv/(wpv+wpu);
                }
                /////////// PHOTONS ////////////
                if (abs(pfc.pdgId())==22) {
                    if(dr < deadcone_ph) continue;
                    iso_ph += wpf*pfc.pt();
                    /////////// NEUTRAL HADRONS ////////////
                } else if (abs(pfc.pdgId())==130) {
                    if(dr < deadcone_nh) continue;
                    iso_nh += wpf*pfc.pt();
                }
            }
            //////////////////  CHARGED from PV  /////////////////////////
        } else if (pfc.fromPV()>1){
            if (abs(pfc.pdgId())==211) {
                if(dr < deadcone_ch) continue;
                iso_ch += pfc.pt();
            }
            //////////////////  CHARGED from PU  /////////////////////////
        } else {
            if (pfc.pt()>ptThresh){
                if(dr < deadcone_pu) continue;
                iso_pu += pfc.pt();
            }
        }
    }
    double iso(0.);
    
    int em = 0;
    if(ptcl->isMuon())
        em = 1;
    
    //double Aeff[2][5] = {{ 0.1013, 0.0988, 0.0572, 0.0842, 0.1530 },{ 0.0913, 0.0765, 0.0546, 0.0728, 0.1177 }};
    double Aeff[2][7] = {{ 0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687 },{ 0.0735, 0.0619, 0.0465, 0.0433, 0.0577,0,0 }};
    
    /*0 ≤ abs(eta) < 1 : 0.1752
    1 ≤ abs(eta) < 1.479 : 0.1862
    1.479 ≤ abs(eta) < 2.0 : 0.1411
    2.0 ≤ abs(eta) < 2.2 : 0.1534
    2.2 ≤ abs(eta) ≤ 2.3 : 0.1903
    2.3 ≤ abs(eta) ≤ 2.4 : 0.2243
    2.4 ≤ abs(eta) ≤ 2.5 : 0.2687
    
    0 ≤ abs(eta) < 0.8 : 0.0735
    0.8 ≤ abs(eta) < 1.3 : 0.0619
    1.3 ≤ abs(eta) < 2.0 : 0.0465
    2.0 ≤ abs(eta) < 2.2 : 0.0433
    2.2 ≤ abs(eta) ≤ 2.5 : 0.0577*/
    
    
    double CorrectedTerm=0.0;
    double riso2 = r_iso*r_iso;
    if(ptcl->isMuon()) {
        if( TMath::Abs( ptcl->eta() ) < 0.8 ) CorrectedTerm = rho * Aeff[em][ 0 ]*(riso2/0.09);
        else if( TMath::Abs( ptcl->eta() ) > 0.8 && TMath::Abs( ptcl->eta() ) < 1.3  )   CorrectedTerm = rho * Aeff[em][ 1 ]*(riso2/0.09);
        else if( TMath::Abs( ptcl->eta() ) > 1.3 && TMath::Abs( ptcl->eta() ) < 2.0  )   CorrectedTerm = rho * Aeff[em][ 2 ]*(riso2/0.09);
        else if( TMath::Abs( ptcl->eta() ) > 2.0 && TMath::Abs( ptcl->eta() ) < 2.2  )   CorrectedTerm = rho * Aeff[em][ 3 ]*(riso2/0.09);
        else if( TMath::Abs( ptcl->eta() ) > 2.2 && TMath::Abs( ptcl->eta() ) < 2.5  )   CorrectedTerm = rho * Aeff[em][ 4 ]*(riso2/0.09);
    } else {
        double elEta = fabs(dynamic_cast<const pat::Electron *>(ptcl)->superCluster()->eta());
        if( elEta < 1.0 ) CorrectedTerm = rho * Aeff[em][ 0 ]*(riso2/0.09);
        else if( elEta > 1.0 && elEta < 1.479  )   CorrectedTerm = rho * Aeff[em][ 1 ]*(riso2/0.09);
        else if( elEta > 1.479 && elEta < 2.0  )   CorrectedTerm = rho * Aeff[em][ 2 ]*(riso2/0.09);
        else if( elEta > 2.0 && elEta < 2.2  )   CorrectedTerm = rho * Aeff[em][ 3 ]*(riso2/0.09);
        else if( elEta > 2.2 && elEta < 2.3  )   CorrectedTerm = rho * Aeff[em][ 4 ]*(riso2/0.09);
        else if( elEta > 2.3 && elEta < 2.4  )   CorrectedTerm = rho * Aeff[em][ 5 ]*(riso2/0.09);
        else if( elEta > 2.4 && elEta < 2.5  )   CorrectedTerm = rho * Aeff[em][ 6 ]*(riso2/0.09);
    }
    /*std::cout<<"pt "<<ptcl->pt();

    std::cout<<"; CorrectedTerm "<<CorrectedTerm;
    std::cout<<"; iso_ph "<<iso_ph;
    std::cout<<"; iso_nh "<<iso_nh;
    std::cout<<"; iso_ch "<<iso_ch;
     */
    
    if (charged_only){
        iso = iso_ch;
    } else {
        //std::cout<<iso_ph<<" "<<iso_nh<<" "<<CorrectedTerm<<" "<<iso_ch<<std::endl;
        //std::cout<<r_iso<<" "<<ptcl->eta()<<" "<<rho * Aeff[em][ 2 ]*(riso2/0.09)<<std::endl;
        //std::cout<<iso_ph/ptcl->pt()<<" "<<iso_nh/ptcl->pt()<<" "<<CorrectedTerm/ptcl->pt()<<" "<<iso_ch/ptcl->pt()<<std::endl;
        iso = iso_ph + iso_nh;
        //if (!use_pfweight) iso -= 0.5*iso_pu;
        iso -= CorrectedTerm; //EA PU correction
        if (iso>0) iso += iso_ch;
        else iso = iso_ch;
    }
    iso = iso/ptcl->pt();
    //std::cout<<"; iso "<<iso<<std::endl;

    
    return iso;
}

std::vector<const pat::Muon* > tools::ssbLooseMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                            double v_muon_pt,
                                                            reco::Vertex::Point PV,
                                                            double v_muon_d0,
                                                           const bool usedz)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    if (!usedz) v_muon_dz = 9999999.;
    
    
    /*int v_muon_numberOfMatchedStations = 2;
     int v_muon_nPixelValidHits = 1;
     int v_muon_numberOfValidMuonHits = 1;
     int v_muon_nValidHits = 6;
     double v_muon_chi2Norm = 10.;*/
    
    //double v_muon_emVetoEt=4.;
    //double v_muon_hadVetoEt= 6.;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
        /*std::cout<<"muon"<<std::endl;
         std::cout<<mu->pt()<<std::endl;
         std::cout<<mu->eta()<<std::endl;
         std::cout<<mu->isGlobalMuon()<<std::endl;
         std::cout<<mu->isPFMuon()<<std::endl;
         std::cout<<mu->numberOfMatchedStations()<<std::endl;
         */
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        
        if ( !(mu->isGlobalMuon() || mu->isTrackerMuon() )) continue;
        if ( !(mu->isPFMuon()) ) continue;
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        /*if ( !mu->isGlobalMuon()  ) continue;
         if ( !mu->isPFMuon() ) continue;
         if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
         
         const reco::TrackRef innerTrack = mu->innerTrack();
         if( innerTrack.isNull() ) continue;
         
         if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
         if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
         
         const reco::TrackRef globalTrack = mu->globalTrack() ;
         if( globalTrack.isNull() ) continue;
         
         if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
         if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
         
         if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
         if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
         //std::cout<<"passed"<<std::endl;*/
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

std::vector<const pat::Muon* > tools::ssbMediumMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    
    
    /*int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;*/
    
    //double v_muon_emVetoEt=4.;
    //double v_muon_hadVetoEt= 6.;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
    {
        /*std::cout<<"muon"<<std::endl;
         std::cout<<mu->pt()<<std::endl;
         std::cout<<mu->eta()<<std::endl;
         std::cout<<mu->isGlobalMuon()<<std::endl;
         std::cout<<mu->isPFMuon()<<std::endl;
         std::cout<<mu->numberOfMatchedStations()<<std::endl;
         */
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        
        if ( !(mu->isGlobalMuon() || mu->isTrackerMuon())) continue;
        if ( !(mu->isPFMuon()) ) continue;

        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;

        
        bool goodGlb = mu->isGlobalMuon() && mu->globalTrack()->normalizedChi2() < 3
        && mu->combinedQuality().chi2LocalPosition < 12 && mu->combinedQuality().trkKink < 20;
        bool good = mu->innerTrack()->validFraction() >= 0.8 && mu->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
        
        if (!good) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        /*if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        //std::cout<<"passed"<<std::endl;*/
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

bool tools::triggerEmulator(const pat::Electron *iE) {
    
    bool passed = true;

    double etasc = TMath::Abs(iE->superCluster()->eta());
    
    if (TMath::Abs(iE->hadronicOverEm()) >=(0.10-0.03*(etasc>1.479))) passed = false;
    if (TMath::Abs(iE->deltaEtaSuperClusterTrackAtVtx()) >= (0.01-0.002*(etasc>1.479))) passed = false;
    if (TMath::Abs(iE->deltaPhiSuperClusterTrackAtVtx()) >= (0.04+0.03*(etasc>1.479))) passed = false;
    double eInvMinusPInv = 1.0/iE->ecalEnergy() - iE->eSuperClusterOverP()/iE->ecalEnergy();
    
    if (eInvMinusPInv<=-0.05) passed = false;
    if (eInvMinusPInv>=(0.01-0.005*(abs(etasc)>1.479))) passed = false;
    if (TMath::Abs(iE->full5x5_sigmaIetaIeta()) >= (0.011+0.019*(etasc>1.479)))passed = false;

    return passed;
}

bool tools::isoTriggerEmulator(const pat::Electron *iE) {
    
    bool passed = true;
    
    //if (iE->pfIsolationVariables().sumChargedHadronPt/iE->pt() > 0.2) passed = false;
    if (iE->dr03TkSumPt()/iE->pt() > 0.2) passed = false;
    else if (iE->hcalPFClusterIso()/iE->pt() > 0.25) passed = false;
    else if (iE->ecalPFClusterIso()/iE->pt() > 0.45) passed = false;
    
    return passed;
}
bool tools::isoTriggerEmulator(const pat::Muon *iM) {
    return true;
}


std::vector<const pat::Electron* > tools::ssbMVAElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                                 reco::BeamSpot::Point BS,
                                                                 bool tight,
                                                                 const bool usedz,
                                                                 const bool missHitsCut)
{
    double v_electron_eta = 2.5;
    double v_electron_dz = 0.1;
    if (!usedz) v_electron_dz = 9999999.;

    //bool bool_electron_ecalDriven = true;
    
    /*double looseMVA[3][2]; //{{0.35, 0.20, -0.52}, {0.73, 0.57, 0.05}};//{0.8, 1.479, };
    looseMVA[0][0] = 0.35;
    looseMVA[1][0] = 0.20;
    looseMVA[2][0] = -0.52;
    looseMVA[0][1] = 0.73;
    looseMVA[1][1] = 0.57;
    looseMVA[2][1] = 0.05;
    
    //double tightMVA[3] = {0.73, 0.57, 0.05};//{0.8, 1.479, };
    
    std::vector<std::string> myManualCatWeigths;
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml");
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml");
    myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml");
    
    vector<string> myManualCatWeigthsTrig;
    string the_path;
    for (unsigned i  = 0 ; i < myManualCatWeigths.size() ; i++){
        the_path = edm::FileInPath ( myManualCatWeigths[i] ).fullPath();
        myManualCatWeigthsTrig.push_back(the_path);
    }
    
    EGammaMvaEleEstimatorCSA14* myMVATrig = new EGammaMvaEleEstimatorCSA14();
    myMVATrig->initialize("BDT",
                          EGammaMvaEleEstimatorCSA14::kNonTrigPhys14,
                          true,
                          myManualCatWeigthsTrig);
     
    */
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
        if (missHitsCut) {
            if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
        }
        
        //comment for SS fakes
        
        /*double mvaValueE = myMVATrig->mvaValue(*el,false);
        
        bool passed = false;
        if (TMath::Abs(el->eta()) < 0.8 ) {
            passed = mvaValueE > looseMVA[0][tight];
        } else if (TMath::Abs(el->eta()) < 1.479 ) {
            passed = mvaValueE > looseMVA[1][tight];
        } else {
            passed = mvaValueE > looseMVA[2][tight];
        }

        
        //std::cout<<"passed"<<std::endl;
        if (passed)*/
            vElectrons.push_back(&*el );
    }
    return vElectrons;
}


bool tools::ssbMVAElectronSelectorPassed(const edm::Ptr<pat::Electron> el,
                                                                 double v_electron_pt,
                                                                 reco::Vertex::Point PV,
                                                                 double v_electron_d0,
                                                                 bool bool_electron_chargeConsistency,
                                                                 edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                                 reco::BeamSpot::Point BS,
                                         bool tight,
                                         const bool usedz,
                                         const bool missHitsCut)
{
    double v_electron_eta = 2.5;
    double v_electron_dz = 0.1;
    if (!usedz) v_electron_dz = 9999999.;

    bool passed = false;
    
    if( el->pt() < v_electron_pt ) return passed;
    if( TMath::Abs(el->eta()) > v_electron_eta ) return passed;
    //if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
    
    const reco::GsfTrackRef gsfTrack = el->gsfTrack();
    
    if (!gsfTrack.isNonnull()) {
        return passed;
    }
    if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  return passed;
    if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) return passed;
    
    if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  return passed;
    
    bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
    if( vtxFitConversion )  return passed;
    
    if (missHitsCut) {
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) return passed;
    }
    
    passed = true;
    return passed;
}



float tools::mass(float pt1 , float pt2, float eta1 , float eta2, float phi1, float phi2)
{
    TLorentzVector v1,v2;
    v1.SetPtEtaPhiM(pt1, eta1, phi1, 0.);
    v2.SetPtEtaPhiM(pt2,eta2, phi2, 0.);
    v1 = v1 + v2;
    return v1.Mag();
}

void tools::ERR( edm::InputTag& IT )
{
    cerr << "[ERROR] : " << IT << " is not a valid input label for this event.  SKIPPING EVENT " << endl;
}

//Muon pfRelIso
double tools::pfRelIso(const pat::Muon *mu)
{
    double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
    double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
    double photonIso = mu->pfIsolationR03().sumPhotonEt;
    double beta = mu->pfIsolationR03().sumPUPt;
    double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
    return pfRelIsoMu;
}

//Electron pfRelIso

/*double tools::pfRelIso(const pat::Electron *iE)
{
    double beta = 0;//iE->pfIsolationR03().sumPUPt;
    double pfRelIsoE  = (iE->chargedHadronIso() + TMath::Max(0.0, iE->neutralHadronIso() + iE->photonIso() - 0.5 * beta ) ) /iE->pt() ;
    return pfRelIsoE;
}*/

double tools::pfRelIso(const pat::Electron *iE)
{
    
    //double pfRelIsoE = (iE->chargedHadronIso() + TMath::Max(0.0, iE->neutralHadronIso() + iE->photonIso() - CorrectedTerm ) ) /iE->pt() ;
    double pfRelIsoE = (iE->pfIsolationVariables().sumChargedHadronPt + TMath::Max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - 0.5*iE->pfIsolationVariables().sumPUPt ) ) /iE->pt() ;
    return pfRelIsoE;
}

double tools::pfRelIso(const pat::Electron *iE, double myRho)
{
    int em = 0;
    double Aeff[2][7] = {{ 0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687 },{ 0.0735, 0.0619, 0.0465, 0.0433, 0.0577,0,0 }};
    
    /*0 ≤ abs(eta) < 1 : 0.1752
     1 ≤ abs(eta) < 1.479 : 0.1862
     1.479 ≤ abs(eta) < 2.0 : 0.1411
     2.0 ≤ abs(eta) < 2.2 : 0.1534
     2.2 ≤ abs(eta) ≤ 2.3 : 0.1903
     2.3 ≤ abs(eta) ≤ 2.4 : 0.2243
     2.4 ≤ abs(eta) ≤ 2.5 : 0.2687
     
     0 ≤ abs(eta) < 0.8 : 0.0735
     0.8 ≤ abs(eta) < 1.3 : 0.0619
     1.3 ≤ abs(eta) < 2.0 : 0.0465
     2.0 ≤ abs(eta) < 2.2 : 0.0433
     2.2 ≤ abs(eta) ≤ 2.5 : 0.0577*/
    
    
    double CorrectedTerm=0.0;
    double elEta = fabs(iE->superCluster()->eta());
    if( elEta < 1.0 ) CorrectedTerm = myRho * Aeff[em][ 0 ];
    else if( elEta > 1.0 && elEta < 1.479  )   CorrectedTerm = myRho * Aeff[em][ 1 ];
    else if( elEta > 1.479 && elEta < 2.0  )   CorrectedTerm = myRho * Aeff[em][ 2 ];
    else if( elEta > 2.0 && elEta < 2.2  )   CorrectedTerm = myRho * Aeff[em][ 3 ];
    else if( elEta > 2.2 && elEta < 2.3  )   CorrectedTerm = myRho * Aeff[em][ 4 ];
    else if( elEta > 2.3 && elEta < 2.4  )   CorrectedTerm = myRho * Aeff[em][ 5 ];
    else if( elEta > 2.4 && elEta < 2.5  )   CorrectedTerm = myRho * Aeff[em][ 6 ];
    
    double pfRelIsoE = (iE->pfIsolationVariables().sumChargedHadronPt + TMath::Max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - CorrectedTerm ) ) /iE->pt() ;

    return pfRelIsoE;
}


double tools::pfRelIso(const pat::Muon *mu, double myRho)
{
    int em = 1;
    double Aeff[2][7] = {{ 0.1752, 0.1862, 0.1411, 0.1534, 0.1903, 0.2243, 0.2687 },{ 0.0735, 0.0619, 0.0465, 0.0433, 0.0577,0,0 }};
    
    double CorrectedTerm=0.0;
    if( TMath::Abs( mu->eta() ) < 0.8 ) CorrectedTerm = myRho * Aeff[em][ 0 ];
    else if( TMath::Abs( mu->eta() ) > 0.8 && TMath::Abs( mu->eta() ) < 1.3  )   CorrectedTerm = myRho * Aeff[em][ 1 ];
    else if( TMath::Abs( mu->eta() ) > 1.3 && TMath::Abs( mu->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[em][ 2 ];
    else if( TMath::Abs( mu->eta() ) > 2.0 && TMath::Abs( mu->eta() ) < 2.2  )   CorrectedTerm = myRho * Aeff[em][ 3 ];
    else if( TMath::Abs( mu->eta() ) > 2.2 && TMath::Abs( mu->eta() ) < 2.5  )   CorrectedTerm = myRho * Aeff[em][ 4 ];

    double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
    double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
    double photonIso = mu->pfIsolationR03().sumPhotonEt;
    //double beta = mu->pfIsolationR03().sumPUPt;
    
    double pfRelIsoE = (chargedHadronIso + TMath::Max(0.0, neutralHadronIso + photonIso - CorrectedTerm ) ) /mu->pt() ;
    
    return pfRelIsoE;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool tools::isTightMuon(const pat::Muon* mu,
                 double v_muon_pt,
                 reco::Vertex::Point PV,
                 double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    //double v_muon_emVetoEt=4.;
    //double v_muon_hadVetoEt= 6.;

    if ( mu->pt()  < v_muon_pt ) return false;
    if ( TMath::Abs( mu->eta() ) > v_muon_eta ) return false;
    //if ( !mu->isTrackerMuon() ) continue;
    if ( !mu->isGlobalMuon()  ) return false;
    if ( !mu->isPFMuon() ) return false;
    if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) return false;   //we should add this to skim
    
    const reco::TrackRef innerTrack = mu->innerTrack();
    if( innerTrack.isNull() ) return false;
    
    /*std::cout<<innerTrack->hitPattern().trackerLayersWithMeasurement()<<std::endl;
     std::cout<<innerTrack->hitPattern().numberOfValidPixelHits()<<std::endl;
     */
    if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) return false;
    if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) return false;
    
    const reco::TrackRef globalTrack = mu->globalTrack() ;
    if( globalTrack.isNull() ) return false;
    
    /*std::cout<<globalTrack->normalizedChi2()<<std::endl;
     std::cout<<globalTrack->hitPattern().numberOfValidMuonHits()<<std::endl;
     
     std::cout<<TMath::Abs(innerTrack->dxy(PV))<<std::endl;
     std::cout<<TMath::Abs(innerTrack->dz(PV))<<std::endl;
     */
    if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) return false;
    if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) return false;
    
    if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) return false;
    if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) return false;
    
    return true;
    
}


std::vector<const pat::Muon* > tools::ssbMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;

    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    //double v_muon_emVetoEt=4.;
    //double v_muon_hadVetoEt= 6.;
    

    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        /*std::cout<<"muon"<<std::endl;
        std::cout<<mu->pt()<<std::endl;
        std::cout<<mu->eta()<<std::endl;
        std::cout<<mu->isGlobalMuon()<<std::endl;
        std::cout<<mu->isPFMuon()<<std::endl;
        std::cout<<mu->numberOfMatchedStations()<<std::endl;
        */
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        
        /*std::cout<<innerTrack->hitPattern().trackerLayersWithMeasurement()<<std::endl;
        std::cout<<innerTrack->hitPattern().numberOfValidPixelHits()<<std::endl;
         */
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        
        /*std::cout<<globalTrack->normalizedChi2()<<std::endl;
        std::cout<<globalTrack->hitPattern().numberOfValidMuonHits()<<std::endl;

        std::cout<<TMath::Abs(innerTrack->dxy(PV))<<std::endl;
        std::cout<<TMath::Abs(innerTrack->dz(PV))<<std::endl;
        */
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        //if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
        //if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;
        
        //std::cout<<"passed"<<std::endl;
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int > tools::ssbMuonSelectorIndex(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    double v_muon_emVetoEt=4.;
    double v_muon_hadVetoEt= 6.;
    
    std::vector<int> vMuons;
    int counter = -1;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        counter++;
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
        if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;
        
        vMuons.push_back(counter);
    }
    
    return vMuons;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Electron* > tools::ssbElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_eta = 2.5;
    double v_electron_dz = 0.1;
    //bool bool_electron_ecalDriven = true;
    

    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        
        /*std::cout<<"electron"<<std::endl;
        std::cout<<el->pt()<<std::endl;
        std::cout<<el->eta()<<std::endl;
        std::cout<<TMath::Abs(el->superCluster()->eta())<<std::endl;
        */
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        //if (!el->trackerDrivenSeed() ) continue;
        
        /*std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
        std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;
        std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;
        */
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;

        /*
        if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        
        //std::cout<<vtxFitConversion<<std::endl;

        if( vtxFitConversion )  continue;
        
        /*std::cout<<gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<<std::endl;
        std::cout<<TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx())<<std::endl;
        std::cout<<TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx())<<std::endl;
        std::cout<<TMath::Abs(el->full5x5_sigmaIetaIeta())<<std::endl;
        std::cout<<TMath::Abs(el->hadronicOverEm())<<std::endl;
        */
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
  
        if( TMath::Abs(el->superCluster()->eta()) < 1.479  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->full5x5_sigmaIetaIeta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
	    }
        else {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->full5x5_sigmaIetaIeta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;   /// at the HLT 0.075  recommended is 0.1
	    }
        
        //std::cout<<"passed"<<std::endl;
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> tools::ssbElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
    //bool bool_electron_ecalDriven = true;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        //if (!el->trackerDrivenSeed() ) continue;

        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        /*if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
      
        if( TMath::Abs( el->eta()) < 1.5  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.1  ) continue;  //recommended is 0.12 but HLT applies 0.1
	    }
        else if( TMath::Abs( el->eta() ) < 2.4 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.075 ) continue;   /// at the HLT 0.075  recommended is 0.1
	    }
        
        vElectrons.push_back(index);
	}
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Muon* > tools::ssbMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0,
                                                      std::vector<const pat::Jet* > SelectedBJets)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    double v_muon_emVetoEt=4.;
    double v_muon_hadVetoEt= 6.;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
        if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;

        TLorentzVector pLep; pLep.SetPtEtaPhiE(mu->pt(), mu->eta(), mu->phi(), mu->energy());
        bool RemoveL = false;
        for (unsigned int nj = 0; nj!=SelectedBJets.size(); ++nj) {
            TLorentzVector pJet; pJet.SetPtEtaPhiM( SelectedBJets[nj]->pt(), SelectedBJets[nj]->eta(), SelectedBJets[nj]->phi(), 0 );
            double ang = pLep.DeltaR( pJet );
            if (ang < 0.4) {
                RemoveL = true;
                break;
            }
        }
        if (!RemoveL)
            vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int > tools::ssbMuonSelectorIndex(const std::vector<pat::Muon>  & thePatMuons,
                                              double v_muon_pt,
                                              reco::Vertex::Point PV,
                                              double v_muon_d0,
                                              std::vector<const pat::Jet* > SelectedBJets)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    double v_muon_emVetoEt=4.;
    double v_muon_hadVetoEt= 6.;
    
    std::vector<int> vMuons;
    int counter = -1;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        counter++;
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
        if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;
        
        TLorentzVector pLep; pLep.SetPtEtaPhiE(mu->pt(), mu->eta(), mu->phi(), mu->energy());
        bool RemoveL = false;
        for (unsigned int nj = 0; nj!=SelectedBJets.size(); ++nj) {
            TLorentzVector pJet; pJet.SetPtEtaPhiM( SelectedBJets[nj]->pt(), SelectedBJets[nj]->eta(), SelectedBJets[nj]->phi(), 0 );
            double ang = pLep.DeltaR( pJet );
            if (ang < 0.4) {
                RemoveL = true;
                break;
            }
        }
        if (!RemoveL)
        vMuons.push_back(counter);
    }
    
    return vMuons;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Electron* > tools::ssbElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS,
                                                              std::vector<const pat::Jet* > SelectedBJets)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
    //bool bool_electron_ecalDriven = true;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        //if (!el->trackerDrivenSeed() ) continue;
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        /*if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;

        if( TMath::Abs( el->eta()) < 1.479  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
	    }
        else if( TMath::Abs( el->eta() ) < 2.4 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;   /// at the HLT 0.075  recommended is 0.1
	    }
        
        TLorentzVector pLep; pLep.SetPtEtaPhiE(el->pt(), el->eta(), el->phi(), el->energy());
        bool RemoveL = false;
        for (unsigned int nj = 0; nj!=SelectedBJets.size(); ++nj) {
            TLorentzVector pJet; pJet.SetPtEtaPhiM( SelectedBJets[nj]->pt(), SelectedBJets[nj]->eta(), SelectedBJets[nj]->phi(), 0 );
            double ang = pLep.DeltaR( pJet );
            if (ang < 0.4) {
                RemoveL = true;
                break;
            }
        }
        if (!RemoveL)
            vElectrons.push_back(&*el );
	}
    return vElectrons;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> tools::ssbElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                 double v_electron_pt,
                                                 reco::Vertex::Point PV,
                                                 double v_electron_d0,
                                                 bool bool_electron_chargeConsistency,
                                                 edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                 reco::BeamSpot::Point BS,
                                                 std::vector<const pat::Jet* > SelectedBJets)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
    //bool bool_electron_ecalDriven = true;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        if (!el->trackerDrivenSeed() ) continue;
        
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        /*if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
 
        if( TMath::Abs( el->eta()) < 1.5  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.1  ) continue;  //recommended is 0.12 but HLT applies 0.1
	    }
        else if( TMath::Abs( el->eta() ) < 2.4 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.075 ) continue;   /// at the HLT 0.075  recommended is 0.1
	    }
        
        TLorentzVector pLep; pLep.SetPtEtaPhiE(el->pt(), el->eta(), el->phi(), el->energy());
        bool RemoveL = false;
        for (unsigned int nj = 0; nj!=SelectedBJets.size(); ++nj) {
            TLorentzVector pJet; pJet.SetPtEtaPhiM( SelectedBJets[nj]->pt(), SelectedBJets[nj]->eta(), SelectedBJets[nj]->phi(), 0 );
            double ang = pLep.DeltaR( pJet );
            if (ang < 0.4) {
                RemoveL = true;
                break;
            }
        }
        if (!RemoveL)
            vElectrons.push_back(index);
	}
    return vElectrons;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Muon Selector

std::vector<const pat::Muon* > tools::fakeMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    //double v_muon_dz = 0.2;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        //if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}


std::vector<int> tools::fakeMuonSelectorIndex(const std::vector<pat::Muon>  & thePatMuons,
                                                       double v_muon_pt,
                                                       reco::Vertex::Point PV,
                                                       double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.2;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    std::vector<int > vMuons;
    int index = -1;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        index++;
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(index);
    }
    
    return vMuons;
}


std::vector<const pat::Muon* > tools::ewkMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    //double v_muon_dz = 0.2; <2012
    double v_muon_dz = 0.5;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
        
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;

        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

std::vector<int> tools::ewkMuonSelectorIndex(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.5;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    int index = -1;
    std::vector<int > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        index++;
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(index);
	}
    return vMuons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<const pat::Electron* > tools::ewkElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.2;
    //bool bool_electron_ecalDriven = true;
    double v_electron_eta=2.4;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        } 
     /*   std::cout<<TMath::Abs(el->eta())<<std::endl;
        std::cout<<TMath::Abs(el->superCluster()->eta())<<std::endl;
        std::cout<<el->pt()<<std::endl;
        std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
        std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;
        std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;
        std::cout<<"Conversion  "<<ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS)<<std::endl;
        std::cout<<gsfTrack->trackerExpectedHitsInner().numberOfHits()<<std::endl;
        
        std::cout<<"********"<<std::endl;*/
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        //if (!el->trackerDrivenSeed() ) continue;

        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;

        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
    
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
                
        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 1 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;

        if( TMath::Abs( el->eta()) < 1.5  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else if( TMath::Abs( el->eta() ) < 2.4 )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
        }
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<int> tools::ewkElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                           double v_electron_pt,
                                                           reco::Vertex::Point PV,
                                                           double v_electron_d0,
                                                           bool bool_electron_chargeConsistency,
                                                           edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                           reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.2;
    //bool bool_electron_ecalDriven = true;
    double v_electron_eta=2.4;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        //if (!el->trackerDrivenSeed() ) continue;

        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 1 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;

        if( TMath::Abs( el->eta()) < 1.5  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
	    }
        else if( TMath::Abs( el->eta() ) < 2.4 )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
	    }
        
        vElectrons.push_back(index);
	}
    return vElectrons;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<const pat::Electron* > tools::fakeElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.1;
    double v_electron_eta=2.5;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        /*   std::cout<<TMath::Abs(el->eta())<<std::endl;
         std::cout<<TMath::Abs(el->superCluster()->eta())<<std::endl;
         std::cout<<el->pt()<<std::endl;
         std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
         std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;
         std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;
         std::cout<<"Conversion  "<<ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS)<<std::endl;
         std::cout<<gsfTrack->trackerExpectedHitsInner().numberOfHits()<<std::endl;
         
         std::cout<<"********"<<std::endl;*/
        
        if( el->pt() < v_electron_pt ) continue;
        
        //std::cout<<TMath::Abs(el->eta())<<std::endl;
        
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;

        //std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
        //std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;

        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;

        //std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;

        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        //std::cout<<vtxFitConversion<<std::endl;
        if( vtxFitConversion )  continue;
        //std::cout<<gsfTrack->trackerExpectedHitsInner().numberOfHits()<<std::endl;
        
        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;

        if( TMath::Abs( el->superCluster()->eta()) < 1.479  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else //if( TMath::Abs( el->superCluster()->eta() ) < v_electron_eta )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
        }
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}


std::vector<const pat::Electron* > tools::phys14LooseElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                               double v_electron_pt,
                                                               reco::Vertex::Point PV,
                                                               double v_electron_d0,
                                                               bool bool_electron_chargeConsistency,
                                                               edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                               reco::BeamSpot::Point BS)
{
    //double v_electron_dz = 0.54342;
    double v_electron_eta=2.5;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        if( el->pt() < v_electron_pt ) continue;

        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;

        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 1 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;

        if( TMath::Abs( el->superCluster()->eta()) < 1.479  )
        {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.072624  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.012442 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.010557 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.121476  ) continue;  //recommended is 0.12 but HLT applies 0.1
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.221803 ) continue;
            if( TMath::Abs(gsfTrack->dz(PV))  > 0.173670  ) continue;
        }
        else //if( TMath::Abs( el->superCluster()->eta() ) < v_electron_eta )
        {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.145129 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.010654 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.032602 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.131862 ) continue;   /// at the HLT 0.075  recommended is 0.1
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.142283 ) continue;
            if( TMath::Abs(gsfTrack->dz(PV))  > 0.198444  ) continue;
        }
        
        vElectrons.push_back(&*el );
    }
    return vElectrons;
}

std::vector<const pat::Electron* > tools::csa14MediumElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                                     double v_electron_pt,
                                                                     reco::Vertex::Point PV,
                                                                     double v_electron_d0,
                                                                     bool bool_electron_chargeConsistency,
                                                                     edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                                     reco::BeamSpot::Point BS)
{
    //double v_electron_dz = 0.54342;
    double v_electron_eta=2.5;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
    {
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        if( el->pt() < v_electron_pt ) continue;
        
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;

        
        if( TMath::Abs( el->superCluster()->eta()) < 1.479  )
        {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.0323  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.0106 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.0107 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.067  ) continue;  //recommended is 0.12 but HLT applies 0.1
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.1043 ) continue;
            if( TMath::Abs(gsfTrack->dz(PV))  > 0.22310  ) continue;
        }
        else //if( TMath::Abs( el->superCluster()->eta() ) < v_electron_eta )
        {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.0455 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.0108 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.0318 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.097 ) continue;   /// at the HLT 0.075  recommended is 0.1
            if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.1201 ) continue;
            if( TMath::Abs(gsfTrack->dz(PV))  > 0.7523  ) continue;
        }
        
        vElectrons.push_back(&*el );
    }
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<int> tools::fakeElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                 double v_electron_pt,
                                                 reco::Vertex::Point PV,
                                                 double v_electron_d0,
                                                 bool bool_electron_chargeConsistency,
                                                 edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                 reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.1;
    double v_electron_eta=2.5;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;
        
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;

        
        if( TMath::Abs( el->superCluster()->eta()) < 1.479  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else //if( TMath::Abs( el->superCluster()->eta() ) < v_electron_eta )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
        }
        
        vElectrons.push_back(index);
	}
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<const pat::Jet* > tools::JetSelectorAll(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta)
{
    bool    bool_jet_id= false;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            //if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() + jet->photonMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        vJets.push_back( &*jet );
        
        /*unsigned int nConst = jet->getPFConstituents().size();
         std::cout<<"Number of constituents "<<nConst<<std::endl;
         for (unsigned int i=0; i!=nConst; ++i) {
         std::cout<<jet->getPFConstituent(i)->reco::LeafCandidate::vz()<<std::endl;
         }*/
        
        
    }
    return vJets;
}

std::vector<int> tools::JetSelectorIndexAll(const std::vector<pat::Jet>  & thePatJets,
                                         double  v_jet_pt,
                                         double  v_jet_eta)
{
    bool    bool_jet_id= false;
    
    std::vector<int> vJets;
    
    int i = -1;
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        i++;
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            //if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() + jet->photonMultiplicity() ) < 2 ) continue;

            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        vJets.push_back( i );
        
    }
    return vJets;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta)
{
    bool    bool_jet_id= true;
    
    //looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4)
    /*NHF = pfjet->neutralHadronEnergyFraction();
    NEMF = pfjet->neutralEmEnergyFraction();
    CHF = pfjet->chargedHadronEnergyFraction();
    MUF = pfjet->muonEnergyFraction();
    CEMF = pfjet->chargedEmEnergyFraction();
    NumConst = pfjet->chargedMultiplicity()+pfjet->neutralMultiplicity();
    CHM = pfjet->chargedMultiplicity();*/
    
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        //std::cout<<"Jet "<<jet->pt()<<" "<<jet->eta()<<std::endl;
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->chargedMultiplicity() + jet->neutralMultiplicity() ) < 2 ) continue;
            //if(  jet->muonEnergyFraction() >=0.8 ) continue;
            //if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() + jet->photonMultiplicity() ) < 2 ) continue;

            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        
        
        /*if( bool_jet_id )
        {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            //if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() + jet->photonMultiplicity() ) < 2 ) continue;
            
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
        }*/
        //std::cout<<"passed "<<std::endl;
        vJets.push_back( &*jet );
        
        /*unsigned int nConst = jet->getPFConstituents().size();
        std::cout<<"Number of constituents "<<nConst<<std::endl;
        for (unsigned int i=0; i!=nConst; ++i) {
            std::cout<<jet->getPFConstituent(i)->reco::LeafCandidate::vz()<<std::endl;
        }*/
        
        
    }
    return vJets;
}

std::vector<int> tools::JetSelectorIndex(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta)
{
    bool    bool_jet_id= true;
    
    std::vector<int> vJets;
    
    int i = -1;
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        i++;
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            //if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() + jet->photonMultiplicity() ) < 2 ) continue;

            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        vJets.push_back( i );
        
    }
    return vJets;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons)
{
    bool    bool_jet_id= true;
    double  v_jet_leptonVetoDR=0.4;
    //double  v_jet_leptonVetoDR = -1.;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_jet_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
	    
        if( vetoJet ) continue;
	    
        
        vJets.push_back( &*jet );
    }
    return vJets;
}

std::vector<int> tools::JetSelectorIndex(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons)
{
    bool    bool_jet_id= true;
    double  v_jet_leptonVetoDR=0.4;
    //double  v_jet_leptonVetoDR = -1.;
    
    std::vector<int> vJets;
    
    int index = -1;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        index++;
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_jet_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
	    
        if( vetoJet ) continue;
	    
        
        vJets.push_back( index );
    }
    return vJets;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons,
                                                 std::map<const reco::PFTau*, int> SelectedTaus)
{
    bool    bool_jet_id= true;
    double  v_jet_leptonVetoDR=0.4;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_jet_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }

        for(std::map<const reco::PFTau*, int >::iterator it = SelectedTaus.begin() ; it != SelectedTaus.end() ;it++ ){
            
            const reco::PFTau *itau = it->first;

            float dphi = TMath::ACos( TMath::Cos( itau->phi()-jet->phi() ) );
            float deta = itau->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
        
        if( vetoJet ) continue;
	    
        
        vJets.push_back( &*jet );
    }
    return vJets;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Photon* > tools::PhotonSelector(const std::vector<pat::Photon>  & thePatPhotons,
                                                 double  v_photon_pt,
                                                 double  v_photon_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons,
                                                 std::map<const reco::PFTau*, int> SelectedTaus)
{
    bool    bool_photon_id= true;
    double  v_photon_leptonVetoDR=0.2;
    
    std::vector< const pat::Photon* > vJets;
    
    for( std::vector<pat::Photon>::const_iterator jet = thePatPhotons.begin(); jet != thePatPhotons.end(); jet++ )
	{
        //std::cout<<jet->pt()<<" "<<jet->eta()<<" "<<jet->hadTowOverEm()<<" "<<jet->sigmaIetaIeta()<<std::endl;
        if( jet->pt() < v_photon_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_photon_eta) continue;
        if( bool_photon_id )
	    {
            if( jet->hadTowOverEm() >= 0.05 ) continue;
            if( TMath::Abs(jet->sigmaIetaIeta()) > 0.012 ) continue;
            if( TMath::Abs( jet->eta() ) < 1.479 )
            {
                if( jet->hadTowOverEm() >= 0.05 ) continue;
                if( TMath::Abs(jet->sigmaIetaIeta()) > 0.034 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_photon_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_photon_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
        
        for(std::map<const reco::PFTau*, int >::iterator it = SelectedTaus.begin() ; it != SelectedTaus.end() ;it++ ){
            
            const reco::PFTau *itau = it->first;
            
            float dphi = TMath::ACos( TMath::Cos( itau->phi()-jet->phi() ) );
            float deta = itau->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_photon_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
        
        if( vetoJet ) continue;
	    
        
        vJets.push_back( &*jet );
    }
    return vJets;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Isolated Muon Selector
std::vector<const pat::Muon* > tools::MuonSelector_Iso(const std::vector< const pat::Muon *>  & thePatMuons, double v_muon_reliso, bool usePFiso, int Iso)
{
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<const pat::Muon *>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        
        float pfRelIsoMu = pfRelIso(*mu);
        float detRelIso = ((*mu)->isolationR03().emEt + (*mu)->isolationR03().hadEt + (*mu)->isolationR03().sumPt) / (*mu)->pt() ;
        float RelIso = detRelIso;
        if(usePFiso) RelIso = pfRelIsoMu;
        if( RelIso  < v_muon_reliso && Iso==0) continue;  //NonIso
        if( RelIso  > v_muon_reliso && Iso==1) continue;  //Iso
        
        vMuons.push_back(*mu);
    }
    
    return vMuons;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Isolated Electron Selector

std::vector<const pat::Electron* > tools::ElectronSelector_Iso(const std::vector<const pat::Electron *>  & thePatElectrons,
                                                               double v_electron_reliso, bool usePFiso, double myRho, int Iso)
{
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<const pat::Electron *>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        
        float ecalIso = TMath::Abs((*el)->eta()) > 1.47 ? (*el)->dr03EcalRecHitSumEt() : TMath::Max((*el)->dr03EcalRecHitSumEt()-1.,0.);
        float detRelIso = ((*el)->dr03TkSumPt() + (*el)->dr03HcalTowerSumEt() + ecalIso ) / (*el)->pt() ;
        float pfRelIsoEl = pfRelIso(*el, myRho);
        float RelIso = detRelIso;
        if(usePFiso) RelIso = pfRelIsoEl;
        if( RelIso  < v_electron_reliso && Iso==0) continue;  //NonIso
        if( RelIso  > v_electron_reliso && Iso==1) continue;  //Iso
        
        vElectrons.push_back(*el );
	}
    return vElectrons;
}

// clean vs selected muons
std::vector<const pat::Electron* > tools::ElectronSelector_Iso(const std::vector<const pat::Electron *>  & thePatElectrons,
                                                               double v_electron_reliso, bool usePFiso, double myRho, int Iso,
                                                               std::vector<const pat::Muon* > & thePatMuons)
{
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<const pat::Electron *>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        
        float ecalIso = TMath::Abs((*el)->eta()) > 1.47 ? (*el)->dr03EcalRecHitSumEt() : TMath::Max((*el)->dr03EcalRecHitSumEt()-1.,0.);
        float detRelIso = ((*el)->dr03TkSumPt() + (*el)->dr03HcalTowerSumEt() + ecalIso ) / (*el)->pt() ;
        float pfRelIsoEl = pfRelIso(*el, myRho);
        float RelIso = detRelIso;
        if(usePFiso) RelIso = pfRelIsoEl;
        if( RelIso  < v_electron_reliso && Iso==0) continue;  //NonIso
        if( RelIso  > v_electron_reliso && Iso==1) continue;  //Iso
        
        bool Remove = false;
        TLorentzVector Ele; Ele.SetPtEtaPhiE( (*el)->pt(), (*el)->eta(), (*el)->phi(), (*el)->energy() );
        for( std::vector<const pat::Muon *>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
            TLorentzVector Mu; Mu.SetPtEtaPhiE( (*mu)->pt(), (*mu)->eta(), (*mu)->phi(), (*mu)->energy() );
            float dR = Mu.DeltaR( Ele );
            if( dR < 0.1 ) {
                Remove = true;
                break;
            }
        }
        if (!Remove)
            vElectrons.push_back(*el );
	}
    return vElectrons;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tools::MT_calc(TLorentzVector Vect, double MET, double MET_Phi){
    
    double MT=sqrt(2* Vect.Pt() * MET * ( 1 - (TMath::Cos(Vect.Phi() - MET_Phi )) ) );
    
    return MT;
}

double tools::Mll_calc(TLorentzVector Vect1, TLorentzVector Vect2){
    return (Vect1 + Vect2).Mag();
}

//ID=#MET+5*(#MT+3*(#Mll+3*#category))
int tools::srID(double met, double mt, double mll, double channel) {
    if ((channel > 0) && (channel!=4))
        return TMath::Min(int(met/50),4) + 5*(int(mt>120)+int(mt>160) + 3*(2*int(mll>100) + 3*channel));
    else 
        return TMath::Min(int(met/50),4) + 5*(int(mt>120)+int(mt>160) + 3*(int(mll>75) + int(mll>105) + 3*channel));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tools::RelIso_El(const pat::Electron *el, int usePFiso, double myRho){
    double reliso= 999.;
    //double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13  };
    double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14  };
    
    float ecalIso = TMath::Abs(el->eta()) > 1.47 ? el->dr03EcalRecHitSumEt() : TMath::Max(el->dr03EcalRecHitSumEt()-1.,0.);
    float detRelIso = (el->dr03TkSumPt() + el->dr03HcalTowerSumEt() + ecalIso ) / el->pt() ;
    
    // PF Isolation
    double CorrectedTerm=0.0;
    if( TMath::Abs( el->superCluster()->eta() ) < 1.0 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
    else CorrectedTerm = myRho * Aeff[ 6 ];
   
    float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
    
    if(usePFiso) reliso = pfRelIso;
    if(!usePFiso) reliso = detRelIso;
    
    
    return reliso;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double tools::RelIso_Mu(const pat::Muon *mu, int usePFiso){
    double reliso= 999.;
    
    float beta_i    = mu->pfIsolationR03().sumPUPt;
    float pfRelIso  = (mu->pfIsolationR03().sumChargedHadronPt +
                       TMath::Max ( 0.0 ,(mu->pfIsolationR03().sumNeutralHadronEt +
                                          mu->pfIsolationR03().sumPhotonEt - 0.5 *beta_i )) ) /mu->pt() ;
    
    float detRelIso = (mu->isolationR03().emEt + mu->isolationR03().hadEt + mu->isolationR03().sumPt) / mu->pt() ;
    
    
    
    if(usePFiso) reliso = pfRelIso;
    if(!usePFiso) reliso = detRelIso;
    
    
    return reliso;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int tools::DY(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons){
    TLorentzVector lep1, lep2, dilep, lep3;  // lepton 4-vectors
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
        const pat::Electron *ele1 = thePatElectrons[i];
        for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
            const pat::Electron *ele2 = thePatElectrons[j];
            if(ele2->charge()==ele1->charge()) continue;
            lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
            lep2.SetPtEtaPhiE(ele2->pt(),ele2->eta(),ele2->phi(),ele2->energy());
            dilep=lep1+lep2;
            if(dilep.M()<12) return 1;
        }
    }
	
    for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
        const pat::Muon *mu1 = thePatMuons[i];
        for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
            const pat::Muon *mu2 = thePatMuons[j];
            if(mu2->charge()==mu1->charge()) continue;
            lep1.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
            lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
            dilep=lep1+lep2;
            if(dilep.M()<12) return 1;
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int tools::OSSF(std::vector<const pat::Electron*>  & thePatElectrons, std::vector<const pat::Muon*>  & thePatMuons){
    int ossf= 0;
    for(unsigned int i = 0   ; i < thePatElectrons.size() ;i++ ) {
        for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
            const pat::Electron *ele1 = thePatElectrons[i];
            const pat::Electron *ele2 = thePatElectrons[j];
            
            if( ( ele1->charge() + ele2->charge() ) ==0 ) ossf = 1;
        }}
    
    for(unsigned int i = 0   ; i < thePatMuons.size() ;i++ ) {
        for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
            const pat::Muon *mu1 = thePatMuons[i];
            const pat::Muon *mu2 = thePatMuons[j];
            
            if( ( mu1->charge() + mu2->charge() ) ==0 ) ossf = 1;
        }}
    
    return ossf;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool tools::PassEventAcceptance(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons){
    int n20=0;
    int n10=0;
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
        const pat::Electron *ele1 = thePatElectrons[i];
        if(ele1->pt()>20){
            n20++;
            n10++;
        }
        else if(ele1->pt()>10) n10++;
    }
    for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
        const pat::Muon *mu1 = thePatMuons[i];
        if(mu1->pt()>20){
            n20++;
            n10++;
        }
        else if(mu1->pt()>10) n10++;
    }
    if(n20>0 && n10==3) return kTRUE;
    else return kFALSE;
}

//*****************************************************************************************************************************
//**** Tau Selector ***********************************************************************************************************
//*****************************************************************************************************************************
std::map<const reco::PFTau*, int > tools::TauSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
                                                    double v_tau_pt,
                                                    double v_tau_eta,
                                                    edm::Handle<reco::PFTauDiscriminator> & electron,
                                                    edm::Handle<reco::PFTauDiscriminator> & muon,
                                                    edm::Handle<reco::PFTauDiscriminator> & iso,
                                                    edm::Handle<reco::PFTauDiscriminator> & decay){
    std::map<const reco::PFTau*, int> vTaus;
    
    for( unsigned  i=0; i<PFTaus->size(); i++ ) {
        
        //std::cout<<"Tau "<<(*PFTaus)[i].pt()<<" "<< (*PFTaus)[i].eta()<<" "<< (*PFTaus)[i].phi()<<std::endl;
        //std::cout<<"pt threshold"<<std::endl;
        
        if((*PFTaus)[i].pt()<v_tau_pt) continue;
        
        //std::cout<<"eta cut"<<std::endl;
        
        if(TMath::Abs((*PFTaus)[i].eta())>v_tau_eta) continue;
        
        reco::PFTauRef tauCandidate(PFTaus, i);
        //std::cout<<"electron discr."<<std::endl;
        if( (*electron)[tauCandidate] < 0.5 ) continue;
        
        //std::cout<<"electron discr."<<std::endl;
        if( (*muon)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"muon discr."<<std::endl;
        if( (*iso)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"iso discr."<<std::endl;
        if( (*decay)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"tau pass"<<std::endl;
        vTaus.insert(std::pair<const reco::PFTau*, int >(&((*PFTaus)[i]), i));
    }
    return vTaus;
}


std::map<const reco::PFTau*, int > tools::TauSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
                                                      double v_tau_pt,
                                                      double v_tau_eta,
                                                      edm::Handle<reco::PFTauDiscriminator> & electron,
                                                      edm::Handle<reco::PFTauDiscriminator> & muon,
                                                      edm::Handle<reco::PFTauDiscriminator> & iso,
                                                      edm::Handle<reco::PFTauDiscriminator> & decay,
                                                      std::vector<const pat::Muon*> & thePatMuons,
                                                      const std::vector<const pat::Electron*>  & thePatElectrons){
    std::map<const reco::PFTau*, int> vTaus;
    
    for( unsigned  i=0; i<PFTaus->size(); i++ ) {
        
        //std::cout<<"Tau "<<(*PFTaus)[i].pt()<<" "<< (*PFTaus)[i].eta()<<" "<< (*PFTaus)[i].phi()<<std::endl;
        //std::cout<<"pt threshold"<<std::endl;
        
        if((*PFTaus)[i].pt()<v_tau_pt) continue;
        
        //std::cout<<"eta cut"<<std::endl;
        
        if(TMath::Abs((*PFTaus)[i].eta())>v_tau_eta) continue;
        
        reco::PFTauRef tauCandidate(PFTaus, i);
        //std::cout<<"electron discr."<<std::endl;
        if( (*electron)[tauCandidate] < 0.5 ) continue;
        
        //std::cout<<"electron discr."<<std::endl;
        if( (*muon)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"muon discr."<<std::endl;
        if( (*iso)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"iso discr."<<std::endl;
        if( (*decay)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"tau pass"<<std::endl;
        
        bool Remove = false;
        TLorentzVector Tau; Tau.SetPtEtaPhiE( (*PFTaus)[i].pt(), (*PFTaus)[i].eta(), (*PFTaus)[i].phi(), (*PFTaus)[i].energy() );
        for( std::vector<const pat::Muon *>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
            TLorentzVector Mu; Mu.SetPtEtaPhiE( (*mu)->pt(), (*mu)->eta(), (*mu)->phi(), (*mu)->energy() );
            float dR = Mu.DeltaR( Tau );
            if( dR < 0.1 ) {
                Remove = true;
                break;
            }
        }
        if (!Remove) {
            for( std::vector<const pat::Electron *>::const_iterator mu = thePatElectrons.begin() ; mu != thePatElectrons.end() ; mu++ ) {
                TLorentzVector Mu; Mu.SetPtEtaPhiE( (*mu)->pt(), (*mu)->eta(), (*mu)->phi(), (*mu)->energy() );
                float dR = Mu.DeltaR( Tau );
                if( dR < 0.1 ) {
                    Remove = true;
                    break;
                }
            }
        }
        if (!Remove)
            vTaus.insert(std::pair<const reco::PFTau*, int >(&((*PFTaus)[i]), i));
    }
    return vTaus;
}


std::vector<const pat::Electron* > tools::ssbElectronVetoSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                                  double Rho,
                                                                  reco::Vertex::Point PV,
                                                                  const char* cutName)
{
    double v_electron_pt = 0 ;
    TString CutName = cutName;
    if( CutName == "gammastar" ) v_electron_pt = 5;
    else if( CutName == "Z"    ) v_electron_pt = 10;
    
    double v_electron_eta = 2.4 ;
    double v_electron_d0 = 0.04;
    double v_electron_reliso = 0.2;
    double v_electron_dz = 0.2;
    //bool bool_electron_ecalDriven = true;
    
    std::vector<const pat::Electron* > vElectrons;
    
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        
        //    if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
        //    if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        if(TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
        if(TMath::Abs(el->gsfTrack()->dz(PV)) > v_electron_dz  ) continue;
        
        if( TMath::Abs( el->eta() ) < 1.4442 ){
            //  if( el->isEB() ){
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.8 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.15 ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else if( TMath::Abs( el->eta() ) < 2.4 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.7 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
        }
        
        float pfRelIsoLocal = pfRelIso(&*el, Rho) ;
        if( pfRelIsoLocal > v_electron_reliso ) continue;
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Muon* > tools::ssbMuonVetoSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                          const char* cutName)
{
    //bool bool_pfIsolation = true;
    float v_muon_pt = 0;
    float v_muon_eta = 2.4;
    //  float v_muon_dz  = 999.;
    float v_muon_iso = 0.2;
    //  bool ZMass = false;
    //double ZMass = 0;
    
    TString CutName = cutName;
    
    if( CutName == "gammastar" ) v_muon_pt = 5;
    else if( CutName == "Z"    ) v_muon_pt = 10;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        if ( mu->isTrackerMuon() || mu->isGlobalMuon() ) {
            if ( !mu->isPFMuon() ) continue;            
            // PF Isolation
            float pfRelIsoLocal = pfRelIso(&*mu) ;
            if( pfRelIsoLocal > v_muon_iso ) continue;
            vMuons.push_back(&*mu);
        }
    }
    return vMuons;
}

bool tools::cleanUp(const pat::Muon* testMu,
                   std::vector<const pat::Muon* > allMu,
                   const char* cutName)
{
    
    TLorentzVector TLVec; TLVec.SetPtEtaPhiE( testMu->pt(), testMu->eta(), testMu->phi(), testMu->energy() );
    TString CutName = cutName;

    for (unsigned int i=0; i!=allMu.size(); ++i) {
        if (testMu->charge() + allMu[i]->charge() != 0) continue;
    
        TLorentzVector P1; P1.SetPtEtaPhiE( allMu[i]->pt(), allMu[i]->eta(), allMu[i]->phi(), allMu[i]->energy() );
        
        float Mass = ( P1 + TLVec ).M();
        if( CutName == "Z" ){
            if( Mass < 106. && Mass > 76. )  {
                return true;
            }
        }
        else if( CutName == "gammastar" ){
            if( Mass < 12. )  {
                return true;
            }
        }
    }
    
    return false;
}

bool tools::cleanUp(const pat::Electron* testMu,
                   std::vector<const pat::Electron* > allMu,
                   const char* cutName)
{
    
    TLorentzVector TLVec; TLVec.SetPtEtaPhiE( testMu->pt(), testMu->eta(), testMu->phi(), testMu->energy() );
    TString CutName = cutName;

    for (unsigned int i=0; i!=allMu.size(); ++i) {
        if (testMu->charge() + allMu[i]->charge() != 0) continue;
        
        TLorentzVector P1; P1.SetPtEtaPhiE( allMu[i]->pt(), allMu[i]->eta(), allMu[i]->phi(), allMu[i]->energy() );
        
        float Mass = ( P1 + TLVec ).M();
        if( CutName == "Z" ){
            if( Mass < 106. && Mass > 76. )  {
                return true;
            }
        }
        else if( CutName == "gammastar" ){
            if( Mass < 12. )  {
                return true;
            }
        }
    }
    
    return false;
}


/*void tools::removeZCandForEE(const std::vector<pat::Electron>  & thePatElectrons,
                             edm::Handle< std::vector<reco::Conversion> > &theConversions,
                             reco::BeamSpot::Point BS,
                             reco::Vertex::Point PV,
                             double  Rho,
                             std::vector< Leptons > & LepCand,
                             bool signalLepOnly,
                             const char* cutName ) {
    
    
    TString CutName = cutName;
    
    double ZMass = 0;
    
    
    std::pair<unsigned int , TString  > myPair;
    std::vector< pair<unsigned int,TString >> VecOfPair;
 
    double v_electron_pt = 0 ;
    
    
    if( CutName == "gammastar" ) v_electron_pt = 5;
    else if( CutName == "Z"    ) v_electron_pt = 10;
    
    
    double v_electron_eta = 2.4 ;
    double v_electron_d0 = 0.04;
    double v_electron_reliso = 0.2;
    double v_electron_dz = 0.2;
    //bool bool_electron_ecalDriven = true;
    
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        
        //    if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
        //    if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        if(TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
        if(TMath::Abs(el->gsfTrack()->dz(PV)) > v_electron_dz  ) continue;
        
        if( TMath::Abs( el->eta() ) < 1.4442 ){
            //  if( el->isEB() ){
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.8 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.15 ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else if( TMath::Abs( el->eta() ) < 2.4 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.7 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
        }
        
                
        float pfRelIsoLocal = pfRelIso(el, Rho) ;
        if( pfRelIsoLocal > v_electron_reliso ) continue;
        
        std::vector< Leptons >::iterator iter;
        for( iter = LepCand.begin() ; iter < LepCand.end(); iter++ ){
            
            TString Type         = iter->type;
            TString Flavor       = iter->flavor;
            int Charge           = iter->charge;
            TLorentzVector TLVec = iter->tlvector;
            
            if( Flavor == "m" ) continue;
            if( signalLepOnly && Type == "sideband" ) continue;
            
            if( el->charge() == Charge  ) continue;
            
            TLorentzVector P1; P1.SetPtEtaPhiE( el->pt(), el->eta(), el->phi(), el->energy() );
            
            float Mass = ( P1 + TLVec ).M();
            
            if( CutName == "Z" ){
                if( Mass < 106. && Mass > 76. )  {
                    ZMass = Mass;
                    myPair.first  = iter->index;
                    myPair.second = iter->flavor;
                    VecOfPair.push_back( myPair );
                }
            }
            else if( CutName == "gammastar" ){
                if( Mass < 12. )  {
                    ZMass = Mass;
                    myPair.first  = iter->index;
                    myPair.second = iter->flavor;
                    VecOfPair.push_back( myPair );
                }
            }
        }
    }
    
    std::vector< Leptons >::iterator iter = LepCand.begin();
    while ( iter != LepCand.end() ) {
        
        TString Type    = iter->type;
        TString Flavor  = iter->flavor;
        unsigned int Index       = iter->index;
        
        //       if( Flavor == "m" ) continue;
        //       if( signalLepOnly && Type == "sideband" ) continue;
        
        bool remove = false;
        std::vector< pair<unsigned int,TString >>::iterator tmp;
        for( tmp = VecOfPair.begin(); tmp < VecOfPair.end(); tmp++ ){
            
            if( (*tmp).first == Index &&  (*tmp).second == Flavor ) remove = true;
        }
        
        if( remove ) iter = LepCand.erase( iter );
        else
            ++iter;
        
    }
}




void tools::removeZCandForMM(const std::vector<pat::Muon>  & thePatMuons,
                             reco::Vertex::Point PV,
                             std::vector< Leptons > & LepCand,
                             bool signalLepOnly,
                             const char* cutName ){
    
    
    //bool bool_pfIsolation = true;
    float v_muon_pt = 0;
    float v_muon_eta = 2.4;
    //  float v_muon_dz  = 999.;
    float v_muon_iso = 0.2;
    //  bool ZMass = false;
    double ZMass = 0;
    
    TString CutName = cutName;
    
    
    if( CutName == "gammastar" ) v_muon_pt = 5;
    else if( CutName == "Z"    ) v_muon_pt = 10;
    
    std::pair<unsigned int , TString  > myPair;
    std::vector< pair< unsigned int,TString >> VecOfPair;
    
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        if ( mu->isTrackerMuon() || mu->isGlobalMuon() ) {
            if ( !mu->isPFMuon() ) continue;
            
            //    TVector3 momPerp(0,0,0);
            //    momPerp.SetPtEtaPhi(mu->pt(),mu->eta(),mu->phi());
            //    TVector3 posPerp(mu->vx()-PV.x(), mu->vy() - PV.y(), 0);
            //    float dzcorr = mu->vz() - PV.z() - posPerp.Dot(momPerp)/mu->pt() * (mu->pz()/mu->pt());
            //   if(TMath::Abs(dzcorr ) > v_muon_dz  ) continue;
            
            // PF Isolation
            float chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
            float neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
            float photonIso = mu->pfIsolationR03().sumPhotonEt;
            float beta = mu->pfIsolationR03().sumPUPt;
            float pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
            
            // Det Isolation
            float detRelIso = (mu->isolationR03().emEt + mu->isolationR03().hadEt + mu->isolationR03().sumPt) / mu->pt() ;
            
            if( bool_pfIsolation ){
                
                if( pfRelIso > v_muon_iso ) continue;
            }
            else{
                
                if( detRelIso > v_muon_iso ) continue;
            }
            
            
            std::vector< Leptons >::iterator iter;
            for( iter = LepCand.begin() ; iter < LepCand.end(); iter++ ){
                
                TString Type         = iter->type;
                TString Flavor       = iter->flavor;
                int Charge           = iter->charge;
                TLorentzVector TLVec = iter->tlvector;
                
                if( signalLepOnly && Type == "sideband" ) continue;
                
                if( Flavor == "e" ) continue;
                if( mu->charge() == Charge  ) continue;
                
                TLorentzVector P1; P1.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
                
                float Mass = ( P1 + TLVec ).M();
                
                
                if( CutName == "Z" ){
                    if( Mass < 106. && Mass > 76. )  {
                        ZMass = Mass;
                        myPair.first  = iter->index;
                        myPair.second = iter->flavor;
                        VecOfPair.push_back( myPair );
                    }
                }
                else if( CutName == "gammastar" ){
                    if( Mass < 12. )  {
                        ZMass = Mass;
                        myPair.first  = iter->index;
                        myPair.second = iter->flavor;
                        VecOfPair.push_back( myPair );
                    }
                }
            }
        }
    }
    
    std::vector< Leptons >::iterator iter = LepCand.begin();
    while ( iter != LepCand.end() ) {
        
        TString Type    = iter->type;
        TString Flavor  = iter->flavor;
        unsigned int Index       = iter->index;
        
        //   if( Flavor == "e" ) continue;
        //   if( signalLepOnly && Type == "sideband" ) continue;
        
        bool remove = false;
        std::vector< pair< unsigned int,TString >>::iterator tmp;
        for( tmp = VecOfPair.begin(); tmp < VecOfPair.end(); tmp++ ){
            
            if( (*tmp).first == Index &&  (*tmp).second == Flavor ) remove = true;
        }
        
        if( remove ) iter = LepCand.erase( iter );
        else
            ++iter;
        
    }
    
}

*/

double tools::JER (double eta) {

    double feta = fabs(eta);
    double scf = 1;
    if (feta < 0.5)
        scf = 1.052;
    else if (feta < 1.1)
        scf = 1.057;
    else if (feta < 1.7)
        scf = 1.096;
    else if (feta < 2.3)
        scf = 1.134;
    else scf = 1.288;
    
    return scf;
}

float tools::quadsum(float a, float b) {
    return sqrt(a*a+b*b);
}

//______________________________________________________________________________
#define JERSUMMER11
float tools::smear_pt_res(float pt, float genpt, float eta)
{
    eta = fabs(eta);
    if (genpt>15. && (fabs(pt - genpt) / pt)<0.5) {  // limit the effect to the core
        double res    = 1.0;
        //double resErr = 0.0;
#ifdef JERSUMMER11
        // from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        if (eta <= 0.5) {
            res    = 1.052;
            //resErr = quadsum(0.012, 0.062);
        } else if (0.5 < eta && eta <= 1.1) {
            res    = 1.057;
            //resErr = quadsum(0.012, 0.062);
        } else if (1.1 < eta && eta <= 1.7) {
            res    = 1.096;
            //resErr = quadsum(0.017, 0.063);
        } else if (1.7 < eta && eta <= 2.3) {
            res    = 1.134;
            //resErr = quadsum(0.035, 0.087);
        } else {
            res    = 1.288;
            //resErr = quadsum(0.127, 0.155);
        }
#else
        // from VHbb analysis
        if (eta <= 1.1) {
            res    = 1.05;
            //resErr = 0.05;
        } else if (1.1 < eta && eta <= 2.5) {
            res    = 1.10;
            //resErr = 0.10;
        } else {
            res    = 1.30;
            //resErr = 0.20;
        }
#endif
        float deltapt = (pt - genpt) * res;
        return TMath::Max(float(0.), genpt + deltapt);
    }
    return pt;
}

bool  tools::srSSbID(int nJets, int nbJets, double met, double Ht, int id) {

    if ((nbJets >= 2 ) && (nJets >= 2) && (Ht > 80) && (met > 30) && (id == 1)) return true;
    if ((nbJets >= 2 ) && (nJets >= 2) && (Ht > 80) && (met > 30) && (id == 2)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 200) && (met > 120) && (id == 3)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 200) && (met > 50) && (id == 4)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 320) && (met > 50) && (id == 5)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 320) && (met > 120) && (id == 6)) return true;
    if ((nbJets >= 3 ) && (nJets >= 3) && (Ht > 200) && (met > 50) && (id == 7)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 320) && (met > 0) && (id == 8)) return true;

    return false;
}

bool  tools::srSSbID28(int nJets, int nbJets, double met, double Ht, int id) {
    
    if ((id == 0) || (id == 9) || (id == 10) || (id == 19) || (id == 20)) return false;
    
    if (nJets < 2) return false;
    if (met < 50) return false;
    if (Ht < 200) return false;
    
    if (!((id < 10) || ((nbJets == 1) && (id/10 == 1)) || ((nbJets > 1) && (id/10 == 2)))) return false;
    if (!(((met < 120) && (id%10 < 5)) || ((met > 120) && (id%10 >=5)))) return false;
    if (!(((Ht < 400) && (id%2 == 1)) || ((Ht > 400) && (id%2 != 1)))) return false;
    if (!(((nJets < 4) && ((((id-1)%10)%4)/2 == 0)) || ((nJets >= 4) && ((((id-1)%10)%4)/2 != 0)))) return false;
    
    return true;

}

//______________________________________________________________________________
#include "TLorentzVector.h"
double tools::evalEt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Et();
}

double tools::evalMt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Mt();
}

std::vector<double> tools::RegressionVars(const pat::Jet *jet, float genpt, const pat::Muon* mu) {
    std::vector<double> vars;
    
    vars.push_back(smear_pt_res((jet->correctedP4("Uncorrected")).Pt(), genpt, jet->eta()));
    vars.push_back(jet->pt());
    vars.push_back(evalEt(jet->pt(), jet->eta(), jet->phi(), jet->energy()));
    vars.push_back(evalMt(jet->pt(), jet->eta(), jet->phi(), jet->energy()));
    
    double hJet_ptLeadTrack = 0;
    const reco::TrackRefVector &tracks =  jet->associatedTracks();
    for (unsigned int k=0; k!= tracks.size(); ++k) {
        if(tracks[k]->pt() > hJet_ptLeadTrack)
            hJet_ptLeadTrack = tracks[k]->pt();
    }
    vars.push_back(hJet_ptLeadTrack);
    
    double hJet_vtx3dL = 0;
    double hJet_vtx3deL = 0;
    double hJet_vtxMass = 0;
    double hJet_vtxPt = 0;
    
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
    } //else
      //  std::cout<<"no info"<<std::endl;
    
    vars.push_back(TMath::Max(0.,hJet_vtx3dL));
    vars.push_back(TMath::Max(0.,hJet_vtx3deL));
    vars.push_back(TMath::Max(0.,hJet_vtxMass));
    vars.push_back(TMath::Max(0.,hJet_vtxPt));
    
    vars.push_back(jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction());
    //vars.push_back(jet->getPFConstituents().size());
    vars.push_back(jet->numberOfDaughters());
    
    vars.push_back(0.); //JECUnc in the main file
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    vars.push_back(pLep.Perp(pJet.Vect()));
    vars.push_back(mu->pt());
    vars.push_back(pLep.DeltaR(pJet));
    
    /*values[0] = "breg_rawptJER := smear_pt_res(hJet_ptRaw, hJet_genPt, hJet_eta)";
    values[1] = "breg_pt := hJet_pt";
    values[2] = "breg_et := evalEt(hJet_pt, hJet_eta, hJet_phi, hJet_e)";
    values[3] = "breg_mt := evalMt(hJet_pt, hJet_eta, hJet_phi, hJet_e)";
    values[4] = "breg_leadtrackpt := hJet_ptLeadTrack";
    values[5] = "breg_vtx3dL := max(0,hJet_vtx3dL)";
    values[6] = "breg_vtx3deL := max(0,hJet_vtx3deL)";
    values[7] = "breg_vtxMass := max(0,hJet_vtxMass)";
    values[8] = "breg_vtxPt := max(0,hJet_vtxPt)";
    values[9] = "breg_cef := hJet_cef";
    values[10] = "breg_ntot := hJet_nconstituents";
    values[11] = "breg_eJEC := hJet_JECUnc";
    values[12] = "breg_softlepptrel := max(0,hJet_SoftLeptptRel*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[13] = "breg_softleppt := max(0,hJet_SoftLeptPt*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[14] = "breg_softlepdR := max(0,hJet_SoftLeptdR*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[15] = "breg_evt_rho25 := rho25";*/
    
    return vars;
    
}

