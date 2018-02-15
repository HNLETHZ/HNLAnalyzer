//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 23 00:14:30 2009 by ROOT version 5.22/00
// from TTree T/UF Analyses ROOT Tree
// found on file: ../test/UFRM.root
//////////////////////////////////////////////////////////

#ifndef AnalysisTreeClass_h
#define AnalysisTreeClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
// Begin list of arrays sizes
const int nRecoBeamSpotsMax = 25;       // maximum number of types of RecoBeamSpots listed in config file
const int nPatElectronCollectionsMax = 25;  // maximum number of patElectron Collections listed in config file
const int nPatElectronMax = 100; // maximum number of patElectrons for each Collection
const int nPatMuonCollectionsMax = 25;    // maximum number of patMuon Collections listed in config file
const int nPatMuonMax = 50;    // maximum number of patMuons for each Collection
const int nPatMuonTrackTypesMax = 3;     // outerTrack, innerTrack, globalTrack
const int nPatMuonTrackMax = nPatMuonTrackTypesMax*nPatMuonMax; // maximum number of tracks for each Collection
const int nPFCandidateMax = 3000; // maximum number of PF Candidates
const int nPatJetCollectionsMax = 25;       // maximum number of types of caloJets listed in config file
const int nPatJetMax = 200;                 // maximum number of caloJets of any particular type
const int nGenMETCollectionsMax = 25;       // maximum number of types of GenMETs listed in config file
const int nHLTBitsMax = 300; // maximum number of HLT bits
const int nGenJetCollectionsMax = 25;       // maximum number of types of genJets listed in config file
const int nGenJetMax = 200;                 // maximum number of genJets of any particular type
const int nPatMETCollectionsMax = 25; // maximum number of types of METs listed in config file
const int nGenParticleMax = 3000; // maximum number of generator particles
const int nGenParticleDaughterMax = nGenParticleMax*10; // maximum number of generator particle daughters
const int nRecoTrackCollectionsMax = 25;       // maximum number of types of tracks listed in config file
const int nRecoTrackMax = 500;                 // maximum number of tracks of any particular type
const int nRecoTrackRecHitMax = 100*nRecoTrackMax; // maximum number of reconstructed hits for one type of tracks
const int nCaloTowerCollectionsMax = 25;   // maximum number of caloTower Collections listed in config file
const int nCaloTowerMax = 1000; // maximum number of caloTowers for each Collection
// End list of arrays sizes

class AnalysisTreeClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           AnalyzerInputParameters_rootMakerVersion;
   Int_t           AnalyzerInputParameters_rootMakerExternalFlag;
   Int_t           EventAuxiliarly_bunchCrossing;
   Int_t           EventAuxiliarly_IdEvent;
   Int_t           EventAuxiliarly_IdRun;
   Bool_t          EventAuxiliarly_isRealData;
   Int_t           EventAuxiliarly_luminosityBlock;
   Int_t           EventAuxiliarly_orbitNumber;
   Float_t         EventParameters_genEventScale_genEventScale;
   Float_t         EventParameters_genEventWeight_genEventWeight;
   Int_t           EventParameters_genEventProcID_genEventProcID;
   Float_t         RecoPdfInfo_genEventPdfInfo_recoPdfInfoScalePDF;
   Int_t           RecoPdfInfo_genEventPdfInfo_recoPdfInfoId1;
   Int_t           RecoPdfInfo_genEventPdfInfo_recoPdfInfoId2;
   Float_t         RecoPdfInfo_genEventPdfInfo_recoPdfInfoX1;
   Float_t         RecoPdfInfo_genEventPdfInfo_recoPdfInfoX2;
   Float_t         RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf1;
   Float_t         RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf2;
   Bool_t          RecoBeamSpot_offlineBeamSpot_isValid;
   Float_t         RecoBeamSpot_offlineBeamSpot_x0;
   Float_t         RecoBeamSpot_offlineBeamSpot_y0;
   Float_t         RecoBeamSpot_offlineBeamSpot_z0;
   Float_t         RecoBeamSpot_offlineBeamSpot_sigmaZ;
   Float_t         RecoBeamSpot_offlineBeamSpot_dxdz;
   Float_t         RecoBeamSpot_offlineBeamSpot_beamWidth;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTIsValid;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTAccept;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTError;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTWasRun;
   Int_t           TriggerResults_TriggerResultsHLT_nHLTBits;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTBitAccept[nHLTBitsMax];   //[TriggerResults_TriggerResultsHLT_nHLTBits]
   Bool_t          TriggerResults_TriggerResultsHLT_HLTBitError[nHLTBitsMax];   //[TriggerResults_TriggerResultsHLT_nHLTBits]
   Bool_t          TriggerResults_TriggerResultsHLT_HLTBitWasRun[nHLTBitsMax];   //[TriggerResults_TriggerResultsHLT_nHLTBits]
   Int_t           MCTruth_genParticles_nGenParticle;
   Int_t           MCTruth_genParticles_genParticlePdgId[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Short_t         MCTruth_genParticles_genParticleStatus[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Short_t         MCTruth_genParticles_genParticleCharge[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Float_t         MCTruth_genParticles_genParticlePx[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Float_t         MCTruth_genParticles_genParticlePy[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Float_t         MCTruth_genParticles_genParticlePz[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Float_t         MCTruth_genParticles_genParticleMass[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Float_t         MCTruth_genParticles_genParticleVx[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Float_t         MCTruth_genParticles_genParticleVy[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Float_t         MCTruth_genParticles_genParticleVz[nGenParticleMax];   //[MCTruth_genParticles_nGenParticle]
   Int_t           MCTruth_genParticles_nGenParticleDaughter;
   Short_t         MCTruth_genParticles_genParticleDaughterIndex[nGenParticleDaughterMax];   //[MCTruth_genParticles_nGenParticleDaughter]
   Short_t         MCTruth_genParticles_genParticleDaughterParentIndex[nGenParticleDaughterMax];   //[MCTruth_genParticles_nGenParticleDaughter]
   Int_t           PF_particleFlow_nPFCandidate;
   Int_t           PF_particleFlow_PFCandidateCharge[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateEnergy[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateECALEnergy[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateHCALEnergy[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateEt[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateEta[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidatePhi[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateMT[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateP[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidatePt[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidatePx[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidatePy[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidatePz[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateVx[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateVy[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Float_t         PF_particleFlow_PFCandidateVz[nPFCandidateMax];   //[PF_particleFlow_nPFCandidate]
   Int_t           PatMu_allLayer1Muons_nPatMuon;
   Float_t         PatMu_allLayer1Muons_trackIso[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_hcalIso[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_ecalIso[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_caloIso[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_trackerIsoDepositEta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_trackerIsoDepositPhi[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_trackerIsoDepositVetoEta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_trackerIsoDepositVetoPhi[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_trackerIsoDepositVetoDR[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_trackerIsoDepositCandEnergy[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_trackerIsoDepositWithin30[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_trackerIsoDepositWithin50[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_trackerIsoDepositCountWithin30[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_trackerIsoDepositCountWithin50[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_hcalIsoDepositEta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_hcalIsoDepositPhi[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_hcalIsoDepositVetoEta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_hcalIsoDepositVetoPhi[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_hcalIsoDepositVetoDR[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_hcalIsoDepositCandEnergy[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_hcalIsoDepositWithin30[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_hcalIsoDepositWithin50[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_hcalIsoDepositCountWithin30[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_hcalIsoDepositCountWithin50[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_ecalIsoDepositEta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_ecalIsoDepositPhi[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_ecalIsoDepositVetoEta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_ecalIsoDepositVetoPhi[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_ecalIsoDepositVetoDR[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_ecalIsoDepositCandEnergy[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_ecalIsoDepositWithin30[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_ecalIsoDepositWithin50[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_ecalIsoDepositCountWithin30[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_ecalIsoDepositCountWithin50[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isMuon[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isCaloMuon[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isTrackerMuon[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isStandAloneMuon[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isGlobalMuon[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isGood[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeAll[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeAllGlobalMuons[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeAllStandAloneMuons[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeAllTrackerMuons[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeTrackerMuonArbitrated[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeAllArbitrated[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeGlobalMuonPromptTight[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeTMLastStationLoose[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeTMLastStationTight[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeTM2DCompatibilityLoose[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeTM2DCompatibilityTight[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeTMOneStationLoose[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeTMOneStationTight[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtLoose[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtTight[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isIsolationValid[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_isolationR03EmEt[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_isolationR03HadEt[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_isolationR03HoEt[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_isolationR03SumPt[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Short_t         PatMu_allLayer1Muons_isolationR03NTracks[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Short_t         PatMu_allLayer1Muons_isolationR03NJets[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_isolationR05EmEt[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_isolationR05HadEt[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_isolationR05HoEt[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_isolationR05SumPt[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Short_t         PatMu_allLayer1Muons_isolationR05NTracks[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Short_t         PatMu_allLayer1Muons_isolationR05NJets[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isEnergyValid[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_calEnergyTower[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_calEnergyTowerS9[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_calEnergyEm[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_calEnergyEmS9[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_calEnergyHad[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_calEnergyHadS9[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_calEnergyHo[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_calEnergyHoS9[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isCaloCompatibilityValid[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_caloCompatibility[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Bool_t          PatMu_allLayer1Muons_isTimeValid[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Short_t         PatMu_allLayer1Muons_timeDirection[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_timeNStations[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_timeInverseBeta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_timeInverseBetaErr[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_timeFreeInverseBeta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_timeFreeInverseBetaErr[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_timeTimeAtIpInOut[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_timeTimeAtIpInOutErr[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_timeTimeAtIpOutIn[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_timeTimeAtIpOutInErr[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Short_t         PatMu_allLayer1Muons_charge[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_px[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_py[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_pz[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_pt[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_p[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_energy[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_phi[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_eta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_theta[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_vx[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_vy[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_vz[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_globalTrackNormalizedChi2[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_innerTrackNumberOfValidHits[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_innerTrackNumberOfLostHits[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_innerTrackPhi[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_innerTrackD0[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Float_t         PatMu_allLayer1Muons_innerTrackD0BeamSpotCorrected[nPatMuonMax];   //[PatMu_allLayer1Muons_nPatMuon]
   Int_t           PatMu_allLayer1Muons_nPatMuonTrack;
   Int_t           PatMu_allLayer1Muons_trackParentPatMuonIndex[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Short_t         PatMu_allLayer1Muons_trackType[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Short_t         PatMu_allLayer1Muons_trackCharge[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackPx[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackPy[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackPz[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackPt[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackPtError[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackP[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackPhi[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackPhiError[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackEta[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackEtaError[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackTheta[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackThetaError[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackVx[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackVy[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackVz[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackChi2[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackChi2Norm[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackNdof[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackD0[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackD0Error[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackDsz[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackDszError[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackDxy[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackDxyError[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackDz[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Float_t         PatMu_allLayer1Muons_trackDzError[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Short_t         PatMu_allLayer1Muons_trackNumberOfValidHits[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Short_t         PatMu_allLayer1Muons_trackNumberOfLostHits[nPatMuonTrackMax];   //[PatMu_allLayer1Muons_nPatMuonTrack]
   Int_t           PatEl_allLayer1Electrons_nPatElectron;
   Bool_t          PatEl_allLayer1Electrons_eidRobustLoose[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Bool_t          PatEl_allLayer1Electrons_eidRobustTight[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Bool_t          PatEl_allLayer1Electrons_eidRobustHighEnergy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Bool_t          PatEl_allLayer1Electrons_eidLoose[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Bool_t          PatEl_allLayer1Electrons_eidTight[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_scSigmaEtaEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_scSigmaIEtaIEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_scE1x5[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_scE2x5Max[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_scE5x5[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackIso[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hcalIso[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_ecalIso[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_caloIso[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackerIsoDepositEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackerIsoDepositPhi[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackerIsoDepositVetoEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackerIsoDepositVetoPhi[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackerIsoDepositVetoDR[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackerIsoDepositCandEnergy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackerIsoDepositWithin30[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackerIsoDepositWithin50[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Int_t           PatEl_allLayer1Electrons_trackerIsoDepositCountWithin30[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Int_t           PatEl_allLayer1Electrons_trackerIsoDepositCountWithin50[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hcalIsoDepositEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hcalIsoDepositPhi[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hcalIsoDepositVetoEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hcalIsoDepositVetoPhi[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hcalIsoDepositVetoDR[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hcalIsoDepositCandEnergy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hcalIsoDepositWithin30[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hcalIsoDepositWithin50[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Int_t           PatEl_allLayer1Electrons_hcalIsoDepositCountWithin30[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Int_t           PatEl_allLayer1Electrons_hcalIsoDepositCountWithin50[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_ecalIsoDepositEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_ecalIsoDepositPhi[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_ecalIsoDepositVetoEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_ecalIsoDepositVetoPhi[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_ecalIsoDepositVetoDR[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_ecalIsoDepositCandEnergy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_ecalIsoDepositWithin30[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_ecalIsoDepositWithin50[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Int_t           PatEl_allLayer1Electrons_ecalIsoDepositCountWithin30[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Int_t           PatEl_allLayer1Electrons_ecalIsoDepositCountWithin50[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Bool_t          PatEl_allLayer1Electrons_isElectron[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Short_t         PatEl_allLayer1Electrons_classification[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_caloEnergy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_caloEnergyError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Bool_t          PatEl_allLayer1Electrons_isEnergyScaleCorrected[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Bool_t          PatEl_allLayer1Electrons_isMomentumCorrected[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_hadronicOverEm[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_caloPositionX[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_caloPositionY[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_caloPositionZ[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumAtVtxPx[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumAtVtxPy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumAtVtxPz[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPositionAtVtxX[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPositionAtVtxY[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPositionAtVtxZ[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumAtCaloPx[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumAtCaloPy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumAtCaloPz[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPositionAtCaloX[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPositionAtCaloY[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPositionAtCaloZ[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumOutPx[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumOutPy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackMomentumOutPz[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_eSuperClusterOverP[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_eSeedClusterOverPout[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_deltaEtaSuperClusterTrackAtVtx[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_deltaEtaSeedClusterTrackAtCalo[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_deltaPhiSuperClusterTrackAtVtx[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_deltaPhiSeedClusterTrackAtCalo[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Int_t           PatEl_allLayer1Electrons_gsfTrackChargeMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackQoverpMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackQoverpModeError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackPxMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackPyMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackPzMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackPtMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackPtModeError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackPMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackPhiMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackPhiModeError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackEtaMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackEtaModeError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackThetaMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackThetaModeError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackLambdaMode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackLambdaModeError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackParameter0Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackParameter1Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackParameter2Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix00Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix01Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix02Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix10Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix11Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix12Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix20Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix21Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix22Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackError00Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackError11Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_gsfTrackError22Mode[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Short_t         PatEl_allLayer1Electrons_trackCharge[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPx[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPz[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPt[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPtError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackP[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPhi[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackPhiError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackEtaError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackTheta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackThetaError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackVx[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackVy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackVz[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackChi2[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackChi2Norm[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackNdof[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackD0[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackD0BeamSpotCorrected[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackD0Error[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackDsz[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackDszError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackDxy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackDxyError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackDz[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_trackDzError[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Short_t         PatEl_allLayer1Electrons_trackNumberOfValidHits[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Short_t         PatEl_allLayer1Electrons_trackNumberOfLostHits[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_superClusterEnergy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_superClusterX[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_superClusterY[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_superClusterZ[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_superClusterEta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_superClusterPhi[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   UInt_t          PatEl_allLayer1Electrons_superClusterCaloID[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_superClusterPreshowerEnergy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_superClusterPhiWidth[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_superClusterEtaWidth[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Short_t         PatEl_allLayer1Electrons_charge[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_px[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_py[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_pz[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_pt[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_p[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_energy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_phi[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_eta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_theta[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_vx[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_vy[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Float_t         PatEl_allLayer1Electrons_vz[nPatElectronMax];   //[PatEl_allLayer1Electrons_nPatElectron]
   Int_t           Track_generalTracks_nRecoTrack;
   Short_t         Track_generalTracks_trackCharge[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackPx[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackPy[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackPz[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackPt[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackPtError[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackP[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackPhi[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackPhiError[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackEta[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackEtaError[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackTheta[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackThetaError[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackVx[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackVy[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackVz[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackChi2[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackChi2Norm[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackNdof[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackD0[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackD0Error[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackDsz[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackDszError[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackDxy[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackDxyError[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackDz[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Float_t         Track_generalTracks_trackDzError[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Short_t         Track_generalTracks_trackNumberOfValidHits[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Short_t         Track_generalTracks_trackNumberOfLostHits[nRecoTrackMax];   //[Track_generalTracks_nRecoTrack]
   Int_t           Track_ckfInOutTracksFromConversions_nRecoTrack;
   Short_t         Track_ckfInOutTracksFromConversions_trackCharge[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackPx[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackPy[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackPz[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackPt[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackPtError[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackP[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackPhi[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackPhiError[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackEta[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackEtaError[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackTheta[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackThetaError[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackVx[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackVy[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackVz[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackChi2[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackChi2Norm[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackNdof[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackD0[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackD0Error[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackDsz[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackDszError[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackDxy[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackDxyError[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackDz[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfInOutTracksFromConversions_trackDzError[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Short_t         Track_ckfInOutTracksFromConversions_trackNumberOfValidHits[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Short_t         Track_ckfInOutTracksFromConversions_trackNumberOfLostHits[nRecoTrackMax];   //[Track_ckfInOutTracksFromConversions_nRecoTrack]
   Int_t           Track_ckfOutInTracksFromConversions_nRecoTrack;
   Short_t         Track_ckfOutInTracksFromConversions_trackCharge[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackPx[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackPy[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackPz[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackPt[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackPtError[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackP[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackPhi[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackPhiError[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackEta[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackEtaError[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackTheta[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackThetaError[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackVx[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackVy[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackVz[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackChi2[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackChi2Norm[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackNdof[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackD0[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackD0Error[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackDsz[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackDszError[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackDxy[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackDxyError[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackDz[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Float_t         Track_ckfOutInTracksFromConversions_trackDzError[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Short_t         Track_ckfOutInTracksFromConversions_trackNumberOfValidHits[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Short_t         Track_ckfOutInTracksFromConversions_trackNumberOfLostHits[nRecoTrackMax];   //[Track_ckfOutInTracksFromConversions_nRecoTrack]
   Int_t           JetGen_iterativeCone5GenJets_nGenJet;
   Short_t         JetGen_iterativeCone5GenJets_genJetCharge[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetPx[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetPy[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetPz[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetPt[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetP[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetPhi[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetEta[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetTheta[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetEmEnergy[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetHadEnergy[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetInvisibleEnergy[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Float_t         JetGen_iterativeCone5GenJets_genJetAuxiliaryEnergy[nGenJetMax];   //[JetGen_iterativeCone5GenJets_nGenJet]
   Int_t           JetGen_sisCone5GenJets_nGenJet;
   Short_t         JetGen_sisCone5GenJets_genJetCharge[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetPx[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetPy[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetPz[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetPt[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetP[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetPhi[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetEta[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetTheta[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetEmEnergy[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetHadEnergy[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetInvisibleEnergy[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Float_t         JetGen_sisCone5GenJets_genJetAuxiliaryEnergy[nGenJetMax];   //[JetGen_sisCone5GenJets_nGenJet]
   Int_t           PatJets_allLayer1JetsIC5JPT_nPatJet;
   Short_t         PatJets_allLayer1JetsIC5JPT_charge[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Bool_t          PatJets_allLayer1JetsIC5JPT_isCaloJet[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Bool_t          PatJets_allLayer1JetsIC5JPT_isPFJet[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Bool_t          PatJets_allLayer1JetsIC5JPT_isBasicJet[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_px[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_py[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_pz[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_pt[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_p[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_phi[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_eta[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_theta[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_maxEInEmTowers[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_maxEInHadTowers[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_hadEnergyInHO[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_hadEnergyInHB[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_hadEnergyInHF[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_hadEnergyInHE[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_emEnergyInEB[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_emEnergyInEE[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_emEnergyInHF[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_energyFractionHadronic[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_emEnergyFraction[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_towersArea[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Short_t         PatJets_allLayer1JetsIC5JPT_n90[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Short_t         PatJets_allLayer1JetsIC5JPT_n60[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Int_t           PatJets_allLayer1JetsIC5JPT_partonFlavour[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexMVABJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_jetBProbabilityBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_jetProbabilityBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_simpleSecondaryVertexBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_softElectronBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_softMuonBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_softMuonByPtBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_softMuonByIP3dBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_trackCountingHighEffBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5JPT_trackCountingHighPurBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5JPT_nPatJet]
   Int_t           PatJets_allLayer1JetsIC5_nPatJet;
   Short_t         PatJets_allLayer1JetsIC5_charge[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Bool_t          PatJets_allLayer1JetsIC5_isCaloJet[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Bool_t          PatJets_allLayer1JetsIC5_isPFJet[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Bool_t          PatJets_allLayer1JetsIC5_isBasicJet[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_px[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_py[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_pz[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_pt[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_p[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_phi[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_eta[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_theta[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_maxEInEmTowers[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_maxEInHadTowers[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_hadEnergyInHO[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_hadEnergyInHB[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_hadEnergyInHF[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_hadEnergyInHE[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_emEnergyInEB[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_emEnergyInEE[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_emEnergyInHF[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_energyFractionHadronic[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_emEnergyFraction[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_towersArea[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Short_t         PatJets_allLayer1JetsIC5_n90[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Short_t         PatJets_allLayer1JetsIC5_n60[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Int_t           PatJets_allLayer1JetsIC5_partonFlavour[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_combinedSecondaryVertexBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_combinedSecondaryVertexMVABJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_jetBProbabilityBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_jetProbabilityBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_simpleSecondaryVertexBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_softElectronBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_softMuonBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_softMuonByPtBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_softMuonByIP3dBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_trackCountingHighEffBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Float_t         PatJets_allLayer1JetsIC5_trackCountingHighPurBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsIC5_nPatJet]
   Int_t           PatJets_allLayer1JetsSC5_nPatJet;
   Short_t         PatJets_allLayer1JetsSC5_charge[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Bool_t          PatJets_allLayer1JetsSC5_isCaloJet[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Bool_t          PatJets_allLayer1JetsSC5_isPFJet[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Bool_t          PatJets_allLayer1JetsSC5_isBasicJet[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_px[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_py[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_pz[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_pt[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_p[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_phi[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_eta[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_theta[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_maxEInEmTowers[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_maxEInHadTowers[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_hadEnergyInHO[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_hadEnergyInHB[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_hadEnergyInHF[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_hadEnergyInHE[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_emEnergyInEB[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_emEnergyInEE[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_emEnergyInHF[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_energyFractionHadronic[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_emEnergyFraction[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_towersArea[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Short_t         PatJets_allLayer1JetsSC5_n90[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Short_t         PatJets_allLayer1JetsSC5_n60[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Int_t           PatJets_allLayer1JetsSC5_partonFlavour[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_combinedSecondaryVertexBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_combinedSecondaryVertexMVABJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_jetBProbabilityBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_jetProbabilityBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_simpleSecondaryVertexBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_softElectronBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_softMuonBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_softMuonByPtBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_softMuonByIP3dBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_trackCountingHighEffBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         PatJets_allLayer1JetsSC5_trackCountingHighPurBJetTags[nPatJetMax];   //[PatJets_allLayer1JetsSC5_nPatJet]
   Float_t         METGen_genMet_genMETPx;
   Float_t         METGen_genMet_genMETPy;
   Float_t         METGen_genMet_genMETEt;
   Float_t         METGen_genMet_genMETPhi;
   Float_t         METGen_genMet_genMETsumEt;
   Float_t         METGen_genMet_genMETmEtSig;
   Float_t         METGen_genMet_genMETE_longitudinal;
   Float_t         METGen_genMet_genMETEmEnergy;
   Float_t         METGen_genMet_genMETHadEnergy;
   Float_t         METGen_genMet_genMETInvisibleEnergy;
   Float_t         METGen_genMet_genMETAuxiliaryEnergy;
   Float_t         METGen_genMetNoNuBSM_genMETPx;
   Float_t         METGen_genMetNoNuBSM_genMETPy;
   Float_t         METGen_genMetNoNuBSM_genMETEt;
   Float_t         METGen_genMetNoNuBSM_genMETPhi;
   Float_t         METGen_genMetNoNuBSM_genMETsumEt;
   Float_t         METGen_genMetNoNuBSM_genMETmEtSig;
   Float_t         METGen_genMetNoNuBSM_genMETE_longitudinal;
   Float_t         METGen_genMetNoNuBSM_genMETEmEnergy;
   Float_t         METGen_genMetNoNuBSM_genMETHadEnergy;
   Float_t         METGen_genMetNoNuBSM_genMETInvisibleEnergy;
   Float_t         METGen_genMetNoNuBSM_genMETAuxiliaryEnergy;
   Float_t         PatMET_allLayer1METsIC5_px;
   Float_t         PatMET_allLayer1METsIC5_py;
   Float_t         PatMET_allLayer1METsIC5_et;
   Float_t         PatMET_allLayer1METsIC5_phi;
   Float_t         PatMET_allLayer1METsIC5_sumEt;
   Float_t         PatMET_allLayer1METsIC5_mEtSig;
   Float_t         PatMET_allLayer1METsIC5_e_longitudinal;
   Int_t           PatMET_allLayer1METsIC5_nCorrections;
   Float_t         PatMET_allLayer1METsIC5_corExUncorrALL;
   Float_t         PatMET_allLayer1METsIC5_corEyUncorrALL;
   Float_t         PatMET_allLayer1METsIC5_corSumEtUncorrALL;
   Float_t         PatMET_allLayer1METsIC5_uncorrectedPtUncorrALL;
   Float_t         PatMET_allLayer1METsIC5_uncorrectedPhiUncorrALL;
   Float_t         PatMET_allLayer1METsIC5_corExUncorrJES;
   Float_t         PatMET_allLayer1METsIC5_corEyUncorrJES;
   Float_t         PatMET_allLayer1METsIC5_corSumEtUncorrJES;
   Float_t         PatMET_allLayer1METsIC5_uncorrectedPtUncorrJES;
   Float_t         PatMET_allLayer1METsIC5_uncorrectedPhiUncorrJES;
   Float_t         PatMET_allLayer1METsIC5_corExUncorrMUON;
   Float_t         PatMET_allLayer1METsIC5_corEyUncorrMUON;
   Float_t         PatMET_allLayer1METsIC5_corSumEtUncorrMUON;
   Float_t         PatMET_allLayer1METsIC5_uncorrectedPtUncorrMUON;
   Float_t         PatMET_allLayer1METsIC5_uncorrectedPhiUncorrMUON;
   Bool_t          PatMET_allLayer1METsIC5_isCaloMET;
   Bool_t          PatMET_allLayer1METsIC5_isRecoMET;
   Float_t         PatMET_allLayer1METsIC5_maxEtInEmTowers;
   Float_t         PatMET_allLayer1METsIC5_maxEtInHadTowers;
   Float_t         PatMET_allLayer1METsIC5_etFractionHadronic;
   Float_t         PatMET_allLayer1METsIC5_emEtFraction;
   Float_t         PatMET_allLayer1METsIC5_hadEtInHB;
   Float_t         PatMET_allLayer1METsIC5_hadEtInHO;
   Float_t         PatMET_allLayer1METsIC5_hadEtInHE;
   Float_t         PatMET_allLayer1METsIC5_hadEtInHF;
   Float_t         PatMET_allLayer1METsIC5_emEtInEB;
   Float_t         PatMET_allLayer1METsIC5_emEtInEE;
   Float_t         PatMET_allLayer1METsIC5_emEtInHF;
   Float_t         PatMET_allLayer1METsIC5_metSignificance;
   Float_t         PatMET_allLayer1METsIC5_CaloSETInpHF;
   Float_t         PatMET_allLayer1METsIC5_CaloSETInmHF;
   Float_t         PatMET_allLayer1METsIC5_CaloMETInpHF;
   Float_t         PatMET_allLayer1METsIC5_CaloMETInmHF;
   Float_t         PatMET_allLayer1METsIC5_CaloMETPhiInpHF;
   Float_t         PatMET_allLayer1METsIC5_CaloMETPhiInmHF;
   Float_t         PatMET_allLayer1METsPF_px;
   Float_t         PatMET_allLayer1METsPF_py;
   Float_t         PatMET_allLayer1METsPF_et;
   Float_t         PatMET_allLayer1METsPF_phi;
   Float_t         PatMET_allLayer1METsPF_sumEt;
   Float_t         PatMET_allLayer1METsPF_mEtSig;
   Float_t         PatMET_allLayer1METsPF_e_longitudinal;
   Int_t           PatMET_allLayer1METsPF_nCorrections;
   Float_t         PatMET_allLayer1METsPF_corExUncorrALL;
   Float_t         PatMET_allLayer1METsPF_corEyUncorrALL;
   Float_t         PatMET_allLayer1METsPF_corSumEtUncorrALL;
   Float_t         PatMET_allLayer1METsPF_uncorrectedPtUncorrALL;
   Float_t         PatMET_allLayer1METsPF_uncorrectedPhiUncorrALL;
   Float_t         PatMET_allLayer1METsPF_corExUncorrJES;
   Float_t         PatMET_allLayer1METsPF_corEyUncorrJES;
   Float_t         PatMET_allLayer1METsPF_corSumEtUncorrJES;
   Float_t         PatMET_allLayer1METsPF_uncorrectedPtUncorrJES;
   Float_t         PatMET_allLayer1METsPF_uncorrectedPhiUncorrJES;
   Float_t         PatMET_allLayer1METsPF_corExUncorrMUON;
   Float_t         PatMET_allLayer1METsPF_corEyUncorrMUON;
   Float_t         PatMET_allLayer1METsPF_corSumEtUncorrMUON;
   Float_t         PatMET_allLayer1METsPF_uncorrectedPtUncorrMUON;
   Float_t         PatMET_allLayer1METsPF_uncorrectedPhiUncorrMUON;
   Bool_t          PatMET_allLayer1METsPF_isCaloMET;
   Bool_t          PatMET_allLayer1METsPF_isRecoMET;
   Float_t         PatMET_allLayer1METsPF_maxEtInEmTowers;
   Float_t         PatMET_allLayer1METsPF_maxEtInHadTowers;
   Float_t         PatMET_allLayer1METsPF_etFractionHadronic;
   Float_t         PatMET_allLayer1METsPF_emEtFraction;
   Float_t         PatMET_allLayer1METsPF_hadEtInHB;
   Float_t         PatMET_allLayer1METsPF_hadEtInHO;
   Float_t         PatMET_allLayer1METsPF_hadEtInHE;
   Float_t         PatMET_allLayer1METsPF_hadEtInHF;
   Float_t         PatMET_allLayer1METsPF_emEtInEB;
   Float_t         PatMET_allLayer1METsPF_emEtInEE;
   Float_t         PatMET_allLayer1METsPF_emEtInHF;
   Float_t         PatMET_allLayer1METsPF_metSignificance;
   Float_t         PatMET_allLayer1METsPF_CaloSETInpHF;
   Float_t         PatMET_allLayer1METsPF_CaloSETInmHF;
   Float_t         PatMET_allLayer1METsPF_CaloMETInpHF;
   Float_t         PatMET_allLayer1METsPF_CaloMETInmHF;
   Float_t         PatMET_allLayer1METsPF_CaloMETPhiInpHF;
   Float_t         PatMET_allLayer1METsPF_CaloMETPhiInmHF;
   Float_t         PatMET_allLayer1METsSC5_px;
   Float_t         PatMET_allLayer1METsSC5_py;
   Float_t         PatMET_allLayer1METsSC5_et;
   Float_t         PatMET_allLayer1METsSC5_phi;
   Float_t         PatMET_allLayer1METsSC5_sumEt;
   Float_t         PatMET_allLayer1METsSC5_mEtSig;
   Float_t         PatMET_allLayer1METsSC5_e_longitudinal;
   Int_t           PatMET_allLayer1METsSC5_nCorrections;
   Float_t         PatMET_allLayer1METsSC5_corExUncorrALL;
   Float_t         PatMET_allLayer1METsSC5_corEyUncorrALL;
   Float_t         PatMET_allLayer1METsSC5_corSumEtUncorrALL;
   Float_t         PatMET_allLayer1METsSC5_uncorrectedPtUncorrALL;
   Float_t         PatMET_allLayer1METsSC5_uncorrectedPhiUncorrALL;
   Float_t         PatMET_allLayer1METsSC5_corExUncorrJES;
   Float_t         PatMET_allLayer1METsSC5_corEyUncorrJES;
   Float_t         PatMET_allLayer1METsSC5_corSumEtUncorrJES;
   Float_t         PatMET_allLayer1METsSC5_uncorrectedPtUncorrJES;
   Float_t         PatMET_allLayer1METsSC5_uncorrectedPhiUncorrJES;
   Float_t         PatMET_allLayer1METsSC5_corExUncorrMUON;
   Float_t         PatMET_allLayer1METsSC5_corEyUncorrMUON;
   Float_t         PatMET_allLayer1METsSC5_corSumEtUncorrMUON;
   Float_t         PatMET_allLayer1METsSC5_uncorrectedPtUncorrMUON;
   Float_t         PatMET_allLayer1METsSC5_uncorrectedPhiUncorrMUON;
   Bool_t          PatMET_allLayer1METsSC5_isCaloMET;
   Bool_t          PatMET_allLayer1METsSC5_isRecoMET;
   Float_t         PatMET_allLayer1METsSC5_maxEtInEmTowers;
   Float_t         PatMET_allLayer1METsSC5_maxEtInHadTowers;
   Float_t         PatMET_allLayer1METsSC5_etFractionHadronic;
   Float_t         PatMET_allLayer1METsSC5_emEtFraction;
   Float_t         PatMET_allLayer1METsSC5_hadEtInHB;
   Float_t         PatMET_allLayer1METsSC5_hadEtInHO;
   Float_t         PatMET_allLayer1METsSC5_hadEtInHE;
   Float_t         PatMET_allLayer1METsSC5_hadEtInHF;
   Float_t         PatMET_allLayer1METsSC5_emEtInEB;
   Float_t         PatMET_allLayer1METsSC5_emEtInEE;
   Float_t         PatMET_allLayer1METsSC5_emEtInHF;
   Float_t         PatMET_allLayer1METsSC5_metSignificance;
   Float_t         PatMET_allLayer1METsSC5_CaloSETInpHF;
   Float_t         PatMET_allLayer1METsSC5_CaloSETInmHF;
   Float_t         PatMET_allLayer1METsSC5_CaloMETInpHF;
   Float_t         PatMET_allLayer1METsSC5_CaloMETInmHF;
   Float_t         PatMET_allLayer1METsSC5_CaloMETPhiInpHF;
   Float_t         PatMET_allLayer1METsSC5_CaloMETPhiInmHF;
   Float_t         PatMET_allLayer1METstcMET_px;
   Float_t         PatMET_allLayer1METstcMET_py;
   Float_t         PatMET_allLayer1METstcMET_et;
   Float_t         PatMET_allLayer1METstcMET_phi;
   Float_t         PatMET_allLayer1METstcMET_sumEt;
   Float_t         PatMET_allLayer1METstcMET_mEtSig;
   Float_t         PatMET_allLayer1METstcMET_e_longitudinal;
   Int_t           PatMET_allLayer1METstcMET_nCorrections;
   Float_t         PatMET_allLayer1METstcMET_corExUncorrALL;
   Float_t         PatMET_allLayer1METstcMET_corEyUncorrALL;
   Float_t         PatMET_allLayer1METstcMET_corSumEtUncorrALL;
   Float_t         PatMET_allLayer1METstcMET_uncorrectedPtUncorrALL;
   Float_t         PatMET_allLayer1METstcMET_uncorrectedPhiUncorrALL;
   Float_t         PatMET_allLayer1METstcMET_corExUncorrJES;
   Float_t         PatMET_allLayer1METstcMET_corEyUncorrJES;
   Float_t         PatMET_allLayer1METstcMET_corSumEtUncorrJES;
   Float_t         PatMET_allLayer1METstcMET_uncorrectedPtUncorrJES;
   Float_t         PatMET_allLayer1METstcMET_uncorrectedPhiUncorrJES;
   Float_t         PatMET_allLayer1METstcMET_corExUncorrMUON;
   Float_t         PatMET_allLayer1METstcMET_corEyUncorrMUON;
   Float_t         PatMET_allLayer1METstcMET_corSumEtUncorrMUON;
   Float_t         PatMET_allLayer1METstcMET_uncorrectedPtUncorrMUON;
   Float_t         PatMET_allLayer1METstcMET_uncorrectedPhiUncorrMUON;
   Bool_t          PatMET_allLayer1METstcMET_isCaloMET;
   Bool_t          PatMET_allLayer1METstcMET_isRecoMET;
   Float_t         PatMET_allLayer1METstcMET_maxEtInEmTowers;
   Float_t         PatMET_allLayer1METstcMET_maxEtInHadTowers;
   Float_t         PatMET_allLayer1METstcMET_etFractionHadronic;
   Float_t         PatMET_allLayer1METstcMET_emEtFraction;
   Float_t         PatMET_allLayer1METstcMET_hadEtInHB;
   Float_t         PatMET_allLayer1METstcMET_hadEtInHO;
   Float_t         PatMET_allLayer1METstcMET_hadEtInHE;
   Float_t         PatMET_allLayer1METstcMET_hadEtInHF;
   Float_t         PatMET_allLayer1METstcMET_emEtInEB;
   Float_t         PatMET_allLayer1METstcMET_emEtInEE;
   Float_t         PatMET_allLayer1METstcMET_emEtInHF;
   Float_t         PatMET_allLayer1METstcMET_metSignificance;
   Float_t         PatMET_allLayer1METstcMET_CaloSETInpHF;
   Float_t         PatMET_allLayer1METstcMET_CaloSETInmHF;
   Float_t         PatMET_allLayer1METstcMET_CaloMETInpHF;
   Float_t         PatMET_allLayer1METstcMET_CaloMETInmHF;
   Float_t         PatMET_allLayer1METstcMET_CaloMETPhiInpHF;
   Float_t         PatMET_allLayer1METstcMET_CaloMETPhiInmHF;
   Int_t           CaloTower_towerMaker_nCaloTower;
   Float_t         CaloTower_towerMaker_emEnergy[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_hadEnergy[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_outerEnergy[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_emEt[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_hadEt[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_outerEt[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_emPositionX[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_emPositionY[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_emPositionZ[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_hadPositionX[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_hadPositionY[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_hadPositionZ[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Int_t           CaloTower_towerMaker_emLvl1[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Int_t           CaloTower_towerMaker_hadLvl1[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_hadEnergyHeOuterLayer[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_hadEnergyHeInnerLayer[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_ecalTime[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_hcalTime[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Int_t           CaloTower_towerMaker_ieta[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Int_t           CaloTower_towerMaker_iphi[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Int_t           CaloTower_towerMaker_numCrystals[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_px[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_py[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_pz[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_pt[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_p[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_energy[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_phi[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_eta[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Float_t         CaloTower_towerMaker_theta[nCaloTowerMax];   //[CaloTower_towerMaker_nCaloTower]
   Bool_t          TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleEle10_LW_OnlyPixelM_L1R;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleEle10_LW_OnlyPixelM_L1R;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleEle10_LW_OnlyPixelM_L1R;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleEle10_LW_OnlyPixelM_L1R;
   Short_t         TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleEle10_LW_OnlyPixelM_L1R;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleMu3;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleMu3;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleMu3;
   Bool_t          TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleMu3;
   Short_t         TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleMu3;

   // List of branches
   TBranch        *b_AnalyzerInputParameters_rootMakerVersion;   //!
   TBranch        *b_AnalyzerInputParameters_rootMakerExternalFlag;   //!
   TBranch        *b_EventAuxiliarly_bunchCrossing;   //!
   TBranch        *b_EventAuxiliarly_IdEvent;   //!
   TBranch        *b_EventAuxiliarly_IdRun;   //!
   TBranch        *b_EventAuxiliarly_isRealData;   //!
   TBranch        *b_EventAuxiliarly_luminosityBlock;   //!
   TBranch        *b_EventAuxiliarly_orbitNumber;   //!
   TBranch        *b_EventParameters_genEventScale_genEventScale;   //!
   TBranch        *b_EventParameters_genEventWeight_genEventWeight;   //!
   TBranch        *b_EventParameters_genEventProcID_genEventProcID;   //!
   TBranch        *b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoScalePDF;   //!
   TBranch        *b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoId1;   //!
   TBranch        *b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoId2;   //!
   TBranch        *b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoX1;   //!
   TBranch        *b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoX2;   //!
   TBranch        *b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf1;   //!
   TBranch        *b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf2;   //!
   TBranch        *b_RecoBeamSpot_offlineBeamSpot_isValid;   //!
   TBranch        *b_RecoBeamSpot_offlineBeamSpot_x0;   //!
   TBranch        *b_RecoBeamSpot_offlineBeamSpot_y0;   //!
   TBranch        *b_RecoBeamSpot_offlineBeamSpot_z0;   //!
   TBranch        *b_RecoBeamSpot_offlineBeamSpot_sigmaZ;   //!
   TBranch        *b_RecoBeamSpot_offlineBeamSpot_dxdz;   //!
   TBranch        *b_RecoBeamSpot_offlineBeamSpot_beamWidth;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTIsValid;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTAccept;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTError;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTWasRun;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_nHLTBits;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTBitAccept;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTBitError;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTBitWasRun;   //!
   TBranch        *b_MCTruth_genParticles_nGenParticle;   //!
   TBranch        *b_MCTruth_genParticles_genParticlePdgId;   //!
   TBranch        *b_MCTruth_genParticles_genParticleStatus;   //!
   TBranch        *b_MCTruth_genParticles_genParticleCharge;   //!
   TBranch        *b_MCTruth_genParticles_genParticlePx;   //!
   TBranch        *b_MCTruth_genParticles_genParticlePy;   //!
   TBranch        *b_MCTruth_genParticles_genParticlePz;   //!
   TBranch        *b_MCTruth_genParticles_genParticleMass;   //!
   TBranch        *b_MCTruth_genParticles_genParticleVx;   //!
   TBranch        *b_MCTruth_genParticles_genParticleVy;   //!
   TBranch        *b_MCTruth_genParticles_genParticleVz;   //!
   TBranch        *b_MCTruth_genParticles_nGenParticleDaughter;   //!
   TBranch        *b_MCTruth_genParticles_genParticleDaughterIndex;   //!
   TBranch        *b_MCTruth_genParticles_genParticleDaughterParentIndex;   //!
   TBranch        *b_PF_particleFlow_nPFCandidate;   //!
   TBranch        *b_PF_particleFlow_PFCandidateCharge;   //!
   TBranch        *b_PF_particleFlow_PFCandidateEnergy;   //!
   TBranch        *b_PF_particleFlow_PFCandidateECALEnergy;   //!
   TBranch        *b_PF_particleFlow_PFCandidateHCALEnergy;   //!
   TBranch        *b_PF_particleFlow_PFCandidateEt;   //!
   TBranch        *b_PF_particleFlow_PFCandidateEta;   //!
   TBranch        *b_PF_particleFlow_PFCandidatePhi;   //!
   TBranch        *b_PF_particleFlow_PFCandidateMT;   //!
   TBranch        *b_PF_particleFlow_PFCandidateP;   //!
   TBranch        *b_PF_particleFlow_PFCandidatePt;   //!
   TBranch        *b_PF_particleFlow_PFCandidatePx;   //!
   TBranch        *b_PF_particleFlow_PFCandidatePy;   //!
   TBranch        *b_PF_particleFlow_PFCandidatePz;   //!
   TBranch        *b_PF_particleFlow_PFCandidateVx;   //!
   TBranch        *b_PF_particleFlow_PFCandidateVy;   //!
   TBranch        *b_PF_particleFlow_PFCandidateVz;   //!
   TBranch        *b_PatMu_allLayer1Muons_nPatMuon;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackIso;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIso;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIso;   //!
   TBranch        *b_PatMu_allLayer1Muons_caloIso;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositEta;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositPhi;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositVetoEta;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositVetoPhi;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositVetoDR;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositCandEnergy;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositWithin30;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositWithin50;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositCountWithin30;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackerIsoDepositCountWithin50;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositEta;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositPhi;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositVetoEta;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositVetoPhi;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositVetoDR;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositCandEnergy;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositWithin30;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositWithin50;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositCountWithin30;   //!
   TBranch        *b_PatMu_allLayer1Muons_hcalIsoDepositCountWithin50;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositEta;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositPhi;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositVetoEta;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositVetoPhi;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositVetoDR;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositCandEnergy;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositWithin30;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositWithin50;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositCountWithin30;   //!
   TBranch        *b_PatMu_allLayer1Muons_ecalIsoDepositCountWithin50;   //!
   TBranch        *b_PatMu_allLayer1Muons_isMuon;   //!
   TBranch        *b_PatMu_allLayer1Muons_isCaloMuon;   //!
   TBranch        *b_PatMu_allLayer1Muons_isTrackerMuon;   //!
   TBranch        *b_PatMu_allLayer1Muons_isStandAloneMuon;   //!
   TBranch        *b_PatMu_allLayer1Muons_isGlobalMuon;   //!
   TBranch        *b_PatMu_allLayer1Muons_isGood;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeAll;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeAllGlobalMuons;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeAllStandAloneMuons;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeAllTrackerMuons;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeTrackerMuonArbitrated;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeAllArbitrated;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeGlobalMuonPromptTight;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeTMLastStationLoose;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeTMLastStationTight;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeTM2DCompatibilityLoose;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeTM2DCompatibilityTight;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeTMOneStationLoose;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeTMOneStationTight;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtTight;   //!
   TBranch        *b_PatMu_allLayer1Muons_isIsolationValid;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR03EmEt;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR03HadEt;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR03HoEt;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR03SumPt;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR03NTracks;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR03NJets;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR05EmEt;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR05HadEt;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR05HoEt;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR05SumPt;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR05NTracks;   //!
   TBranch        *b_PatMu_allLayer1Muons_isolationR05NJets;   //!
   TBranch        *b_PatMu_allLayer1Muons_isEnergyValid;   //!
   TBranch        *b_PatMu_allLayer1Muons_calEnergyTower;   //!
   TBranch        *b_PatMu_allLayer1Muons_calEnergyTowerS9;   //!
   TBranch        *b_PatMu_allLayer1Muons_calEnergyEm;   //!
   TBranch        *b_PatMu_allLayer1Muons_calEnergyEmS9;   //!
   TBranch        *b_PatMu_allLayer1Muons_calEnergyHad;   //!
   TBranch        *b_PatMu_allLayer1Muons_calEnergyHadS9;   //!
   TBranch        *b_PatMu_allLayer1Muons_calEnergyHo;   //!
   TBranch        *b_PatMu_allLayer1Muons_calEnergyHoS9;   //!
   TBranch        *b_PatMu_allLayer1Muons_isCaloCompatibilityValid;   //!
   TBranch        *b_PatMu_allLayer1Muons_caloCompatibility;   //!
   TBranch        *b_PatMu_allLayer1Muons_isTimeValid;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeDirection;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeNStations;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeInverseBeta;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeInverseBetaErr;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeFreeInverseBeta;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeFreeInverseBetaErr;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeTimeAtIpInOut;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeTimeAtIpInOutErr;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeTimeAtIpOutIn;   //!
   TBranch        *b_PatMu_allLayer1Muons_timeTimeAtIpOutInErr;   //!
   TBranch        *b_PatMu_allLayer1Muons_charge;   //!
   TBranch        *b_PatMu_allLayer1Muons_px;   //!
   TBranch        *b_PatMu_allLayer1Muons_py;   //!
   TBranch        *b_PatMu_allLayer1Muons_pz;   //!
   TBranch        *b_PatMu_allLayer1Muons_pt;   //!
   TBranch        *b_PatMu_allLayer1Muons_p;   //!
   TBranch        *b_PatMu_allLayer1Muons_energy;   //!
   TBranch        *b_PatMu_allLayer1Muons_phi;   //!
   TBranch        *b_PatMu_allLayer1Muons_eta;   //!
   TBranch        *b_PatMu_allLayer1Muons_theta;   //!
   TBranch        *b_PatMu_allLayer1Muons_vx;   //!
   TBranch        *b_PatMu_allLayer1Muons_vy;   //!
   TBranch        *b_PatMu_allLayer1Muons_vz;   //!
   TBranch        *b_PatMu_allLayer1Muons_globalTrackNormalizedChi2;   //!
   TBranch        *b_PatMu_allLayer1Muons_innerTrackNumberOfValidHits;   //!
   TBranch        *b_PatMu_allLayer1Muons_innerTrackNumberOfLostHits;   //!
   TBranch        *b_PatMu_allLayer1Muons_innerTrackPhi;   //!
   TBranch        *b_PatMu_allLayer1Muons_innerTrackD0;   //!
   TBranch        *b_PatMu_allLayer1Muons_innerTrackD0BeamSpotCorrected;   //!
   TBranch        *b_PatMu_allLayer1Muons_nPatMuonTrack;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackParentPatMuonIndex;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackType;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackCharge;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackPx;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackPy;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackPz;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackPt;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackPtError;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackP;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackPhi;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackPhiError;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackEta;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackEtaError;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackTheta;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackThetaError;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackVx;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackVy;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackVz;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackChi2;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackChi2Norm;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackNdof;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackD0;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackD0Error;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackDsz;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackDszError;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackDxy;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackDxyError;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackDz;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackDzError;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackNumberOfValidHits;   //!
   TBranch        *b_PatMu_allLayer1Muons_trackNumberOfLostHits;   //!
   TBranch        *b_PatEl_allLayer1Electrons_nPatElectron;   //!
   TBranch        *b_PatEl_allLayer1Electrons_eidRobustLoose;   //!
   TBranch        *b_PatEl_allLayer1Electrons_eidRobustTight;   //!
   TBranch        *b_PatEl_allLayer1Electrons_eidRobustHighEnergy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_eidLoose;   //!
   TBranch        *b_PatEl_allLayer1Electrons_eidTight;   //!
   TBranch        *b_PatEl_allLayer1Electrons_scSigmaEtaEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_scSigmaIEtaIEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_scE1x5;   //!
   TBranch        *b_PatEl_allLayer1Electrons_scE2x5Max;   //!
   TBranch        *b_PatEl_allLayer1Electrons_scE5x5;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackIso;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIso;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIso;   //!
   TBranch        *b_PatEl_allLayer1Electrons_caloIso;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositPhi;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositVetoEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositVetoPhi;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositVetoDR;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositCandEnergy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositWithin30;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositWithin50;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositCountWithin30;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackerIsoDepositCountWithin50;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositPhi;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositVetoEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositVetoPhi;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositVetoDR;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositCandEnergy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositWithin30;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositWithin50;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositCountWithin30;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hcalIsoDepositCountWithin50;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositPhi;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositVetoEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositVetoPhi;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositVetoDR;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositCandEnergy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositWithin30;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositWithin50;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositCountWithin30;   //!
   TBranch        *b_PatEl_allLayer1Electrons_ecalIsoDepositCountWithin50;   //!
   TBranch        *b_PatEl_allLayer1Electrons_isElectron;   //!
   TBranch        *b_PatEl_allLayer1Electrons_classification;   //!
   TBranch        *b_PatEl_allLayer1Electrons_caloEnergy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_caloEnergyError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_isEnergyScaleCorrected;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_isMomentumCorrected;   //!
   TBranch        *b_PatEl_allLayer1Electrons_hadronicOverEm;   //!
   TBranch        *b_PatEl_allLayer1Electrons_caloPositionX;   //!
   TBranch        *b_PatEl_allLayer1Electrons_caloPositionY;   //!
   TBranch        *b_PatEl_allLayer1Electrons_caloPositionZ;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumAtVtxPx;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumAtVtxPy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumAtVtxPz;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPositionAtVtxX;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPositionAtVtxY;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPositionAtVtxZ;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumAtCaloPx;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumAtCaloPy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumAtCaloPz;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPositionAtCaloX;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPositionAtCaloY;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPositionAtCaloZ;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumOutPx;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumOutPy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackMomentumOutPz;   //!
   TBranch        *b_PatEl_allLayer1Electrons_eSuperClusterOverP;   //!
   TBranch        *b_PatEl_allLayer1Electrons_eSeedClusterOverPout;   //!
   TBranch        *b_PatEl_allLayer1Electrons_deltaEtaSuperClusterTrackAtVtx;   //!
   TBranch        *b_PatEl_allLayer1Electrons_deltaEtaSeedClusterTrackAtCalo;   //!
   TBranch        *b_PatEl_allLayer1Electrons_deltaPhiSuperClusterTrackAtVtx;   //!
   TBranch        *b_PatEl_allLayer1Electrons_deltaPhiSeedClusterTrackAtCalo;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackChargeMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackQoverpMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackQoverpModeError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackPxMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackPyMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackPzMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackPtMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackPtModeError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackPMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackPhiMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackPhiModeError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackEtaMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackEtaModeError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackThetaMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackThetaModeError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackLambdaMode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackLambdaModeError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackParameter0Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackParameter1Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackParameter2Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix00Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix01Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix02Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix10Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix11Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix12Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix20Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix21Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix22Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackError00Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackError11Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_gsfTrackError22Mode;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackCharge;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPx;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPz;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPt;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPtError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackP;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPhi;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackPhiError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackEtaError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackTheta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackThetaError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackVx;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackVy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackVz;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackChi2;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackChi2Norm;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackNdof;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackD0;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackD0BeamSpotCorrected;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackD0Error;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackDsz;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackDszError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackDxy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackDxyError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackDz;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackDzError;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackNumberOfValidHits;   //!
   TBranch        *b_PatEl_allLayer1Electrons_trackNumberOfLostHits;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterEnergy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterX;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterY;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterZ;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterEta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterPhi;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterCaloID;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterPreshowerEnergy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterPhiWidth;   //!
   TBranch        *b_PatEl_allLayer1Electrons_superClusterEtaWidth;   //!
   TBranch        *b_PatEl_allLayer1Electrons_charge;   //!
   TBranch        *b_PatEl_allLayer1Electrons_px;   //!
   TBranch        *b_PatEl_allLayer1Electrons_py;   //!
   TBranch        *b_PatEl_allLayer1Electrons_pz;   //!
   TBranch        *b_PatEl_allLayer1Electrons_pt;   //!
   TBranch        *b_PatEl_allLayer1Electrons_p;   //!
   TBranch        *b_PatEl_allLayer1Electrons_energy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_phi;   //!
   TBranch        *b_PatEl_allLayer1Electrons_eta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_theta;   //!
   TBranch        *b_PatEl_allLayer1Electrons_vx;   //!
   TBranch        *b_PatEl_allLayer1Electrons_vy;   //!
   TBranch        *b_PatEl_allLayer1Electrons_vz;   //!
   TBranch        *b_Track_generalTracks_nRecoTrack;   //!
   TBranch        *b_Track_generalTracks_trackCharge;   //!
   TBranch        *b_Track_generalTracks_trackPx;   //!
   TBranch        *b_Track_generalTracks_trackPy;   //!
   TBranch        *b_Track_generalTracks_trackPz;   //!
   TBranch        *b_Track_generalTracks_trackPt;   //!
   TBranch        *b_Track_generalTracks_trackPtError;   //!
   TBranch        *b_Track_generalTracks_trackP;   //!
   TBranch        *b_Track_generalTracks_trackPhi;   //!
   TBranch        *b_Track_generalTracks_trackPhiError;   //!
   TBranch        *b_Track_generalTracks_trackEta;   //!
   TBranch        *b_Track_generalTracks_trackEtaError;   //!
   TBranch        *b_Track_generalTracks_trackTheta;   //!
   TBranch        *b_Track_generalTracks_trackThetaError;   //!
   TBranch        *b_Track_generalTracks_trackVx;   //!
   TBranch        *b_Track_generalTracks_trackVy;   //!
   TBranch        *b_Track_generalTracks_trackVz;   //!
   TBranch        *b_Track_generalTracks_trackChi2;   //!
   TBranch        *b_Track_generalTracks_trackChi2Norm;   //!
   TBranch        *b_Track_generalTracks_trackNdof;   //!
   TBranch        *b_Track_generalTracks_trackD0;   //!
   TBranch        *b_Track_generalTracks_trackD0Error;   //!
   TBranch        *b_Track_generalTracks_trackDsz;   //!
   TBranch        *b_Track_generalTracks_trackDszError;   //!
   TBranch        *b_Track_generalTracks_trackDxy;   //!
   TBranch        *b_Track_generalTracks_trackDxyError;   //!
   TBranch        *b_Track_generalTracks_trackDz;   //!
   TBranch        *b_Track_generalTracks_trackDzError;   //!
   TBranch        *b_Track_generalTracks_trackNumberOfValidHits;   //!
   TBranch        *b_Track_generalTracks_trackNumberOfLostHits;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_nRecoTrack;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackCharge;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackPx;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackPy;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackPz;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackPt;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackPtError;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackP;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackPhi;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackPhiError;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackEta;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackEtaError;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackTheta;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackThetaError;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackVx;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackVy;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackVz;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackChi2;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackChi2Norm;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackNdof;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackD0;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackD0Error;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackDsz;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackDszError;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackDxy;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackDxyError;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackDz;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackDzError;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackNumberOfValidHits;   //!
   TBranch        *b_Track_ckfInOutTracksFromConversions_trackNumberOfLostHits;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_nRecoTrack;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackCharge;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackPx;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackPy;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackPz;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackPt;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackPtError;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackP;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackPhi;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackPhiError;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackEta;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackEtaError;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackTheta;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackThetaError;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackVx;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackVy;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackVz;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackChi2;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackChi2Norm;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackNdof;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackD0;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackD0Error;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackDsz;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackDszError;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackDxy;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackDxyError;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackDz;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackDzError;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackNumberOfValidHits;   //!
   TBranch        *b_Track_ckfOutInTracksFromConversions_trackNumberOfLostHits;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_nGenJet;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetCharge;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetPx;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetPy;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetPz;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetPt;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetP;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetPhi;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetEta;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetTheta;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetEmEnergy;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetHadEnergy;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetInvisibleEnergy;   //!
   TBranch        *b_JetGen_iterativeCone5GenJets_genJetAuxiliaryEnergy;   //!
   TBranch        *b_JetGen_sisCone5GenJets_nGenJet;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetCharge;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetPx;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetPy;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetPz;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetPt;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetP;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetPhi;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetEta;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetTheta;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetEmEnergy;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetHadEnergy;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetInvisibleEnergy;   //!
   TBranch        *b_JetGen_sisCone5GenJets_genJetAuxiliaryEnergy;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_nPatJet;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_charge;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_isCaloJet;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_isPFJet;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_isBasicJet;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_px;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_py;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_pz;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_pt;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_p;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_phi;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_eta;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_theta;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_maxEInEmTowers;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_maxEInHadTowers;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_hadEnergyInHO;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_hadEnergyInHB;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_hadEnergyInHF;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_hadEnergyInHE;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_emEnergyInEB;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_emEnergyInEE;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_emEnergyInHF;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_energyFractionHadronic;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_emEnergyFraction;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_towersArea;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_n90;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_n60;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_partonFlavour;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexMVABJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_jetBProbabilityBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_jetProbabilityBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_simpleSecondaryVertexBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_softElectronBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_softMuonBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_softMuonByPtBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_softMuonByIP3dBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_trackCountingHighEffBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5JPT_trackCountingHighPurBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_nPatJet;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_charge;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_isCaloJet;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_isPFJet;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_isBasicJet;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_px;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_py;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_pz;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_pt;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_p;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_phi;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_eta;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_theta;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_maxEInEmTowers;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_maxEInHadTowers;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_hadEnergyInHO;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_hadEnergyInHB;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_hadEnergyInHF;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_hadEnergyInHE;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_emEnergyInEB;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_emEnergyInEE;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_emEnergyInHF;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_energyFractionHadronic;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_emEnergyFraction;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_towersArea;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_n90;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_n60;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_partonFlavour;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_combinedSecondaryVertexBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_combinedSecondaryVertexMVABJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_jetBProbabilityBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_jetProbabilityBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_simpleSecondaryVertexBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_softElectronBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_softMuonBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_softMuonByPtBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_softMuonByIP3dBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_trackCountingHighEffBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsIC5_trackCountingHighPurBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_nPatJet;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_charge;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_isCaloJet;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_isPFJet;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_isBasicJet;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_px;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_py;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_pz;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_pt;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_p;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_phi;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_eta;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_theta;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_maxEInEmTowers;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_maxEInHadTowers;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_hadEnergyInHO;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_hadEnergyInHB;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_hadEnergyInHF;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_hadEnergyInHE;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_emEnergyInEB;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_emEnergyInEE;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_emEnergyInHF;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_energyFractionHadronic;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_emEnergyFraction;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_towersArea;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_n90;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_n60;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_partonFlavour;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_combinedSecondaryVertexBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_combinedSecondaryVertexMVABJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_jetBProbabilityBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_jetProbabilityBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_simpleSecondaryVertexBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_softElectronBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_softMuonBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_softMuonByPtBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_softMuonByIP3dBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_trackCountingHighEffBJetTags;   //!
   TBranch        *b_PatJets_allLayer1JetsSC5_trackCountingHighPurBJetTags;   //!
   TBranch        *b_METGen_genMet_genMETPx;   //!
   TBranch        *b_METGen_genMet_genMETPy;   //!
   TBranch        *b_METGen_genMet_genMETEt;   //!
   TBranch        *b_METGen_genMet_genMETPhi;   //!
   TBranch        *b_METGen_genMet_genMETsumEt;   //!
   TBranch        *b_METGen_genMet_genMETmEtSig;   //!
   TBranch        *b_METGen_genMet_genMETE_longitudinal;   //!
   TBranch        *b_METGen_genMet_genMETEmEnergy;   //!
   TBranch        *b_METGen_genMet_genMETHadEnergy;   //!
   TBranch        *b_METGen_genMet_genMETInvisibleEnergy;   //!
   TBranch        *b_METGen_genMet_genMETAuxiliaryEnergy;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETPx;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETPy;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETEt;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETPhi;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETsumEt;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETmEtSig;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETE_longitudinal;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETEmEnergy;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETHadEnergy;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETInvisibleEnergy;   //!
   TBranch        *b_METGen_genMetNoNuBSM_genMETAuxiliaryEnergy;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_px;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_py;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_et;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_phi;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_sumEt;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_mEtSig;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_e_longitudinal;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_nCorrections;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_corExUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_corEyUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_corSumEtUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_uncorrectedPtUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_uncorrectedPhiUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_corExUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_corEyUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_corSumEtUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_uncorrectedPtUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_uncorrectedPhiUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_corExUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_corEyUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_corSumEtUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_uncorrectedPtUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_uncorrectedPhiUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_isCaloMET;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_isRecoMET;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_maxEtInEmTowers;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_maxEtInHadTowers;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_etFractionHadronic;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_emEtFraction;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_hadEtInHB;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_hadEtInHO;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_hadEtInHE;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_hadEtInHF;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_emEtInEB;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_emEtInEE;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_emEtInHF;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_metSignificance;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_CaloSETInpHF;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_CaloSETInmHF;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_CaloMETInpHF;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_CaloMETInmHF;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_CaloMETPhiInpHF;   //!
   TBranch        *b_PatMET_allLayer1METsIC5_CaloMETPhiInmHF;   //!
   TBranch        *b_PatMET_allLayer1METsPF_px;   //!
   TBranch        *b_PatMET_allLayer1METsPF_py;   //!
   TBranch        *b_PatMET_allLayer1METsPF_et;   //!
   TBranch        *b_PatMET_allLayer1METsPF_phi;   //!
   TBranch        *b_PatMET_allLayer1METsPF_sumEt;   //!
   TBranch        *b_PatMET_allLayer1METsPF_mEtSig;   //!
   TBranch        *b_PatMET_allLayer1METsPF_e_longitudinal;   //!
   TBranch        *b_PatMET_allLayer1METsPF_nCorrections;   //!
   TBranch        *b_PatMET_allLayer1METsPF_corExUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsPF_corEyUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsPF_corSumEtUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsPF_uncorrectedPtUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsPF_uncorrectedPhiUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsPF_corExUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsPF_corEyUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsPF_corSumEtUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsPF_uncorrectedPtUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsPF_uncorrectedPhiUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsPF_corExUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsPF_corEyUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsPF_corSumEtUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsPF_uncorrectedPtUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsPF_uncorrectedPhiUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsPF_isCaloMET;   //!
   TBranch        *b_PatMET_allLayer1METsPF_isRecoMET;   //!
   TBranch        *b_PatMET_allLayer1METsPF_maxEtInEmTowers;   //!
   TBranch        *b_PatMET_allLayer1METsPF_maxEtInHadTowers;   //!
   TBranch        *b_PatMET_allLayer1METsPF_etFractionHadronic;   //!
   TBranch        *b_PatMET_allLayer1METsPF_emEtFraction;   //!
   TBranch        *b_PatMET_allLayer1METsPF_hadEtInHB;   //!
   TBranch        *b_PatMET_allLayer1METsPF_hadEtInHO;   //!
   TBranch        *b_PatMET_allLayer1METsPF_hadEtInHE;   //!
   TBranch        *b_PatMET_allLayer1METsPF_hadEtInHF;   //!
   TBranch        *b_PatMET_allLayer1METsPF_emEtInEB;   //!
   TBranch        *b_PatMET_allLayer1METsPF_emEtInEE;   //!
   TBranch        *b_PatMET_allLayer1METsPF_emEtInHF;   //!
   TBranch        *b_PatMET_allLayer1METsPF_metSignificance;   //!
   TBranch        *b_PatMET_allLayer1METsPF_CaloSETInpHF;   //!
   TBranch        *b_PatMET_allLayer1METsPF_CaloSETInmHF;   //!
   TBranch        *b_PatMET_allLayer1METsPF_CaloMETInpHF;   //!
   TBranch        *b_PatMET_allLayer1METsPF_CaloMETInmHF;   //!
   TBranch        *b_PatMET_allLayer1METsPF_CaloMETPhiInpHF;   //!
   TBranch        *b_PatMET_allLayer1METsPF_CaloMETPhiInmHF;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_px;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_py;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_et;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_phi;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_sumEt;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_mEtSig;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_e_longitudinal;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_nCorrections;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_corExUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_corEyUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_corSumEtUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_uncorrectedPtUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_uncorrectedPhiUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_corExUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_corEyUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_corSumEtUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_uncorrectedPtUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_uncorrectedPhiUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_corExUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_corEyUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_corSumEtUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_uncorrectedPtUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_uncorrectedPhiUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_isCaloMET;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_isRecoMET;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_maxEtInEmTowers;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_maxEtInHadTowers;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_etFractionHadronic;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_emEtFraction;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_hadEtInHB;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_hadEtInHO;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_hadEtInHE;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_hadEtInHF;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_emEtInEB;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_emEtInEE;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_emEtInHF;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_metSignificance;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_CaloSETInpHF;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_CaloSETInmHF;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_CaloMETInpHF;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_CaloMETInmHF;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_CaloMETPhiInpHF;   //!
   TBranch        *b_PatMET_allLayer1METsSC5_CaloMETPhiInmHF;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_px;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_py;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_et;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_phi;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_sumEt;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_mEtSig;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_e_longitudinal;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_nCorrections;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_corExUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_corEyUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_corSumEtUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_uncorrectedPtUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_uncorrectedPhiUncorrALL;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_corExUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_corEyUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_corSumEtUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_uncorrectedPtUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_uncorrectedPhiUncorrJES;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_corExUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_corEyUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_corSumEtUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_uncorrectedPtUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_uncorrectedPhiUncorrMUON;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_isCaloMET;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_isRecoMET;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_maxEtInEmTowers;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_maxEtInHadTowers;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_etFractionHadronic;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_emEtFraction;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_hadEtInHB;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_hadEtInHO;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_hadEtInHE;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_hadEtInHF;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_emEtInEB;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_emEtInEE;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_emEtInHF;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_metSignificance;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_CaloSETInpHF;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_CaloSETInmHF;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_CaloMETInpHF;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_CaloMETInmHF;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_CaloMETPhiInpHF;   //!
   TBranch        *b_PatMET_allLayer1METstcMET_CaloMETPhiInmHF;   //!
   TBranch        *b_CaloTower_towerMaker_nCaloTower;   //!
   TBranch        *b_CaloTower_towerMaker_emEnergy;   //!
   TBranch        *b_CaloTower_towerMaker_hadEnergy;   //!
   TBranch        *b_CaloTower_towerMaker_outerEnergy;   //!
   TBranch        *b_CaloTower_towerMaker_emEt;   //!
   TBranch        *b_CaloTower_towerMaker_hadEt;   //!
   TBranch        *b_CaloTower_towerMaker_outerEt;   //!
   TBranch        *b_CaloTower_towerMaker_emPositionX;   //!
   TBranch        *b_CaloTower_towerMaker_emPositionY;   //!
   TBranch        *b_CaloTower_towerMaker_emPositionZ;   //!
   TBranch        *b_CaloTower_towerMaker_hadPositionX;   //!
   TBranch        *b_CaloTower_towerMaker_hadPositionY;   //!
   TBranch        *b_CaloTower_towerMaker_hadPositionZ;   //!
   TBranch        *b_CaloTower_towerMaker_emLvl1;   //!
   TBranch        *b_CaloTower_towerMaker_hadLvl1;   //!
   TBranch        *b_CaloTower_towerMaker_hadEnergyHeOuterLayer;   //!
   TBranch        *b_CaloTower_towerMaker_hadEnergyHeInnerLayer;   //!
   TBranch        *b_CaloTower_towerMaker_ecalTime;   //!
   TBranch        *b_CaloTower_towerMaker_hcalTime;   //!
   TBranch        *b_CaloTower_towerMaker_ieta;   //!
   TBranch        *b_CaloTower_towerMaker_iphi;   //!
   TBranch        *b_CaloTower_towerMaker_numCrystals;   //!
   TBranch        *b_CaloTower_towerMaker_px;   //!
   TBranch        *b_CaloTower_towerMaker_py;   //!
   TBranch        *b_CaloTower_towerMaker_pz;   //!
   TBranch        *b_CaloTower_towerMaker_pt;   //!
   TBranch        *b_CaloTower_towerMaker_p;   //!
   TBranch        *b_CaloTower_towerMaker_energy;   //!
   TBranch        *b_CaloTower_towerMaker_phi;   //!
   TBranch        *b_CaloTower_towerMaker_eta;   //!
   TBranch        *b_CaloTower_towerMaker_theta;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleEle10_LW_OnlyPixelM_L1R;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleEle10_LW_OnlyPixelM_L1R;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleEle10_LW_OnlyPixelM_L1R;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleEle10_LW_OnlyPixelM_L1R;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleEle10_LW_OnlyPixelM_L1R;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleMu3;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleMu3;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleMu3;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleMu3;   //!
   TBranch        *b_TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleMu3;   //!

   AnalysisTreeClass(TTree *tree=0);
   virtual ~AnalysisTreeClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalysisTreeClass_cxx
AnalysisTreeClass::AnalysisTreeClass(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../test/UFRM.root");
      if (!f) {
         f = new TFile("../test/UFRM.root");
      }
      tree = (TTree*)gDirectory->Get("T");

   }
   Init(tree);
}

AnalysisTreeClass::~AnalysisTreeClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalysisTreeClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalysisTreeClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalysisTreeClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("AnalyzerInputParameters_rootMakerVersion", &AnalyzerInputParameters_rootMakerVersion, &b_AnalyzerInputParameters_rootMakerVersion);
   fChain->SetBranchAddress("AnalyzerInputParameters_rootMakerExternalFlag", &AnalyzerInputParameters_rootMakerExternalFlag, &b_AnalyzerInputParameters_rootMakerExternalFlag);
   fChain->SetBranchAddress("EventAuxiliarly_bunchCrossing", &EventAuxiliarly_bunchCrossing, &b_EventAuxiliarly_bunchCrossing);
   fChain->SetBranchAddress("EventAuxiliarly_IdEvent", &EventAuxiliarly_IdEvent, &b_EventAuxiliarly_IdEvent);
   fChain->SetBranchAddress("EventAuxiliarly_IdRun", &EventAuxiliarly_IdRun, &b_EventAuxiliarly_IdRun);
   fChain->SetBranchAddress("EventAuxiliarly_isRealData", &EventAuxiliarly_isRealData, &b_EventAuxiliarly_isRealData);
   fChain->SetBranchAddress("EventAuxiliarly_luminosityBlock", &EventAuxiliarly_luminosityBlock, &b_EventAuxiliarly_luminosityBlock);
   fChain->SetBranchAddress("EventAuxiliarly_orbitNumber", &EventAuxiliarly_orbitNumber, &b_EventAuxiliarly_orbitNumber);
   fChain->SetBranchAddress("EventParameters_genEventScale_genEventScale", &EventParameters_genEventScale_genEventScale, &b_EventParameters_genEventScale_genEventScale);
   fChain->SetBranchAddress("EventParameters_genEventWeight_genEventWeight", &EventParameters_genEventWeight_genEventWeight, &b_EventParameters_genEventWeight_genEventWeight);
   fChain->SetBranchAddress("EventParameters_genEventProcID_genEventProcID", &EventParameters_genEventProcID_genEventProcID, &b_EventParameters_genEventProcID_genEventProcID);
   fChain->SetBranchAddress("RecoPdfInfo_genEventPdfInfo_recoPdfInfoScalePDF", &RecoPdfInfo_genEventPdfInfo_recoPdfInfoScalePDF, &b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoScalePDF);
   fChain->SetBranchAddress("RecoPdfInfo_genEventPdfInfo_recoPdfInfoId1", &RecoPdfInfo_genEventPdfInfo_recoPdfInfoId1, &b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoId1);
   fChain->SetBranchAddress("RecoPdfInfo_genEventPdfInfo_recoPdfInfoId2", &RecoPdfInfo_genEventPdfInfo_recoPdfInfoId2, &b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoId2);
   fChain->SetBranchAddress("RecoPdfInfo_genEventPdfInfo_recoPdfInfoX1", &RecoPdfInfo_genEventPdfInfo_recoPdfInfoX1, &b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoX1);
   fChain->SetBranchAddress("RecoPdfInfo_genEventPdfInfo_recoPdfInfoX2", &RecoPdfInfo_genEventPdfInfo_recoPdfInfoX2, &b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoX2);
   fChain->SetBranchAddress("RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf1", &RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf1, &b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf1);
   fChain->SetBranchAddress("RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf2", &RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf2, &b_RecoPdfInfo_genEventPdfInfo_recoPdfInfoPdf2);
   fChain->SetBranchAddress("RecoBeamSpot_offlineBeamSpot_isValid", &RecoBeamSpot_offlineBeamSpot_isValid, &b_RecoBeamSpot_offlineBeamSpot_isValid);
   fChain->SetBranchAddress("RecoBeamSpot_offlineBeamSpot_x0", &RecoBeamSpot_offlineBeamSpot_x0, &b_RecoBeamSpot_offlineBeamSpot_x0);
   fChain->SetBranchAddress("RecoBeamSpot_offlineBeamSpot_y0", &RecoBeamSpot_offlineBeamSpot_y0, &b_RecoBeamSpot_offlineBeamSpot_y0);
   fChain->SetBranchAddress("RecoBeamSpot_offlineBeamSpot_z0", &RecoBeamSpot_offlineBeamSpot_z0, &b_RecoBeamSpot_offlineBeamSpot_z0);
   fChain->SetBranchAddress("RecoBeamSpot_offlineBeamSpot_sigmaZ", &RecoBeamSpot_offlineBeamSpot_sigmaZ, &b_RecoBeamSpot_offlineBeamSpot_sigmaZ);
   fChain->SetBranchAddress("RecoBeamSpot_offlineBeamSpot_dxdz", &RecoBeamSpot_offlineBeamSpot_dxdz, &b_RecoBeamSpot_offlineBeamSpot_dxdz);
   fChain->SetBranchAddress("RecoBeamSpot_offlineBeamSpot_beamWidth", &RecoBeamSpot_offlineBeamSpot_beamWidth, &b_RecoBeamSpot_offlineBeamSpot_beamWidth);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTIsValid", &TriggerResults_TriggerResultsHLT_HLTIsValid, &b_TriggerResults_TriggerResultsHLT_HLTIsValid);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTAccept", &TriggerResults_TriggerResultsHLT_HLTAccept, &b_TriggerResults_TriggerResultsHLT_HLTAccept);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTError", &TriggerResults_TriggerResultsHLT_HLTError, &b_TriggerResults_TriggerResultsHLT_HLTError);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTWasRun", &TriggerResults_TriggerResultsHLT_HLTWasRun, &b_TriggerResults_TriggerResultsHLT_HLTWasRun);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_nHLTBits", &TriggerResults_TriggerResultsHLT_nHLTBits, &b_TriggerResults_TriggerResultsHLT_nHLTBits);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTBitAccept", TriggerResults_TriggerResultsHLT_HLTBitAccept, &b_TriggerResults_TriggerResultsHLT_HLTBitAccept);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTBitError", TriggerResults_TriggerResultsHLT_HLTBitError, &b_TriggerResults_TriggerResultsHLT_HLTBitError);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTBitWasRun", TriggerResults_TriggerResultsHLT_HLTBitWasRun, &b_TriggerResults_TriggerResultsHLT_HLTBitWasRun);
   fChain->SetBranchAddress("MCTruth_genParticles_nGenParticle", &MCTruth_genParticles_nGenParticle, &b_MCTruth_genParticles_nGenParticle);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticlePdgId", MCTruth_genParticles_genParticlePdgId, &b_MCTruth_genParticles_genParticlePdgId);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticleStatus", MCTruth_genParticles_genParticleStatus, &b_MCTruth_genParticles_genParticleStatus);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticleCharge", MCTruth_genParticles_genParticleCharge, &b_MCTruth_genParticles_genParticleCharge);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticlePx", MCTruth_genParticles_genParticlePx, &b_MCTruth_genParticles_genParticlePx);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticlePy", MCTruth_genParticles_genParticlePy, &b_MCTruth_genParticles_genParticlePy);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticlePz", MCTruth_genParticles_genParticlePz, &b_MCTruth_genParticles_genParticlePz);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticleMass", MCTruth_genParticles_genParticleMass, &b_MCTruth_genParticles_genParticleMass);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticleVx", MCTruth_genParticles_genParticleVx, &b_MCTruth_genParticles_genParticleVx);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticleVy", MCTruth_genParticles_genParticleVy, &b_MCTruth_genParticles_genParticleVy);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticleVz", MCTruth_genParticles_genParticleVz, &b_MCTruth_genParticles_genParticleVz);
   fChain->SetBranchAddress("MCTruth_genParticles_nGenParticleDaughter", &MCTruth_genParticles_nGenParticleDaughter, &b_MCTruth_genParticles_nGenParticleDaughter);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticleDaughterIndex", MCTruth_genParticles_genParticleDaughterIndex, &b_MCTruth_genParticles_genParticleDaughterIndex);
   fChain->SetBranchAddress("MCTruth_genParticles_genParticleDaughterParentIndex", MCTruth_genParticles_genParticleDaughterParentIndex, &b_MCTruth_genParticles_genParticleDaughterParentIndex);
   fChain->SetBranchAddress("PF_particleFlow_nPFCandidate", &PF_particleFlow_nPFCandidate, &b_PF_particleFlow_nPFCandidate);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateCharge", PF_particleFlow_PFCandidateCharge, &b_PF_particleFlow_PFCandidateCharge);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateEnergy", PF_particleFlow_PFCandidateEnergy, &b_PF_particleFlow_PFCandidateEnergy);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateECALEnergy", PF_particleFlow_PFCandidateECALEnergy, &b_PF_particleFlow_PFCandidateECALEnergy);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateHCALEnergy", PF_particleFlow_PFCandidateHCALEnergy, &b_PF_particleFlow_PFCandidateHCALEnergy);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateEt", PF_particleFlow_PFCandidateEt, &b_PF_particleFlow_PFCandidateEt);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateEta", PF_particleFlow_PFCandidateEta, &b_PF_particleFlow_PFCandidateEta);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidatePhi", PF_particleFlow_PFCandidatePhi, &b_PF_particleFlow_PFCandidatePhi);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateMT", PF_particleFlow_PFCandidateMT, &b_PF_particleFlow_PFCandidateMT);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateP", PF_particleFlow_PFCandidateP, &b_PF_particleFlow_PFCandidateP);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidatePt", PF_particleFlow_PFCandidatePt, &b_PF_particleFlow_PFCandidatePt);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidatePx", PF_particleFlow_PFCandidatePx, &b_PF_particleFlow_PFCandidatePx);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidatePy", PF_particleFlow_PFCandidatePy, &b_PF_particleFlow_PFCandidatePy);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidatePz", PF_particleFlow_PFCandidatePz, &b_PF_particleFlow_PFCandidatePz);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateVx", PF_particleFlow_PFCandidateVx, &b_PF_particleFlow_PFCandidateVx);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateVy", PF_particleFlow_PFCandidateVy, &b_PF_particleFlow_PFCandidateVy);
   fChain->SetBranchAddress("PF_particleFlow_PFCandidateVz", PF_particleFlow_PFCandidateVz, &b_PF_particleFlow_PFCandidateVz);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_nPatMuon", &PatMu_allLayer1Muons_nPatMuon, &b_PatMu_allLayer1Muons_nPatMuon);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackIso", PatMu_allLayer1Muons_trackIso, &b_PatMu_allLayer1Muons_trackIso);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIso", PatMu_allLayer1Muons_hcalIso, &b_PatMu_allLayer1Muons_hcalIso);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIso", PatMu_allLayer1Muons_ecalIso, &b_PatMu_allLayer1Muons_ecalIso);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_caloIso", PatMu_allLayer1Muons_caloIso, &b_PatMu_allLayer1Muons_caloIso);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositEta", PatMu_allLayer1Muons_trackerIsoDepositEta, &b_PatMu_allLayer1Muons_trackerIsoDepositEta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositPhi", PatMu_allLayer1Muons_trackerIsoDepositPhi, &b_PatMu_allLayer1Muons_trackerIsoDepositPhi);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositVetoEta", PatMu_allLayer1Muons_trackerIsoDepositVetoEta, &b_PatMu_allLayer1Muons_trackerIsoDepositVetoEta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositVetoPhi", PatMu_allLayer1Muons_trackerIsoDepositVetoPhi, &b_PatMu_allLayer1Muons_trackerIsoDepositVetoPhi);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositVetoDR", PatMu_allLayer1Muons_trackerIsoDepositVetoDR, &b_PatMu_allLayer1Muons_trackerIsoDepositVetoDR);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositCandEnergy", PatMu_allLayer1Muons_trackerIsoDepositCandEnergy, &b_PatMu_allLayer1Muons_trackerIsoDepositCandEnergy);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositWithin30", PatMu_allLayer1Muons_trackerIsoDepositWithin30, &b_PatMu_allLayer1Muons_trackerIsoDepositWithin30);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositWithin50", PatMu_allLayer1Muons_trackerIsoDepositWithin50, &b_PatMu_allLayer1Muons_trackerIsoDepositWithin50);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositCountWithin30", PatMu_allLayer1Muons_trackerIsoDepositCountWithin30, &b_PatMu_allLayer1Muons_trackerIsoDepositCountWithin30);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackerIsoDepositCountWithin50", PatMu_allLayer1Muons_trackerIsoDepositCountWithin50, &b_PatMu_allLayer1Muons_trackerIsoDepositCountWithin50);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositEta", PatMu_allLayer1Muons_hcalIsoDepositEta, &b_PatMu_allLayer1Muons_hcalIsoDepositEta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositPhi", PatMu_allLayer1Muons_hcalIsoDepositPhi, &b_PatMu_allLayer1Muons_hcalIsoDepositPhi);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositVetoEta", PatMu_allLayer1Muons_hcalIsoDepositVetoEta, &b_PatMu_allLayer1Muons_hcalIsoDepositVetoEta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositVetoPhi", PatMu_allLayer1Muons_hcalIsoDepositVetoPhi, &b_PatMu_allLayer1Muons_hcalIsoDepositVetoPhi);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositVetoDR", PatMu_allLayer1Muons_hcalIsoDepositVetoDR, &b_PatMu_allLayer1Muons_hcalIsoDepositVetoDR);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositCandEnergy", PatMu_allLayer1Muons_hcalIsoDepositCandEnergy, &b_PatMu_allLayer1Muons_hcalIsoDepositCandEnergy);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositWithin30", PatMu_allLayer1Muons_hcalIsoDepositWithin30, &b_PatMu_allLayer1Muons_hcalIsoDepositWithin30);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositWithin50", PatMu_allLayer1Muons_hcalIsoDepositWithin50, &b_PatMu_allLayer1Muons_hcalIsoDepositWithin50);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositCountWithin30", PatMu_allLayer1Muons_hcalIsoDepositCountWithin30, &b_PatMu_allLayer1Muons_hcalIsoDepositCountWithin30);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_hcalIsoDepositCountWithin50", PatMu_allLayer1Muons_hcalIsoDepositCountWithin50, &b_PatMu_allLayer1Muons_hcalIsoDepositCountWithin50);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositEta", PatMu_allLayer1Muons_ecalIsoDepositEta, &b_PatMu_allLayer1Muons_ecalIsoDepositEta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositPhi", PatMu_allLayer1Muons_ecalIsoDepositPhi, &b_PatMu_allLayer1Muons_ecalIsoDepositPhi);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositVetoEta", PatMu_allLayer1Muons_ecalIsoDepositVetoEta, &b_PatMu_allLayer1Muons_ecalIsoDepositVetoEta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositVetoPhi", PatMu_allLayer1Muons_ecalIsoDepositVetoPhi, &b_PatMu_allLayer1Muons_ecalIsoDepositVetoPhi);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositVetoDR", PatMu_allLayer1Muons_ecalIsoDepositVetoDR, &b_PatMu_allLayer1Muons_ecalIsoDepositVetoDR);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositCandEnergy", PatMu_allLayer1Muons_ecalIsoDepositCandEnergy, &b_PatMu_allLayer1Muons_ecalIsoDepositCandEnergy);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositWithin30", PatMu_allLayer1Muons_ecalIsoDepositWithin30, &b_PatMu_allLayer1Muons_ecalIsoDepositWithin30);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositWithin50", PatMu_allLayer1Muons_ecalIsoDepositWithin50, &b_PatMu_allLayer1Muons_ecalIsoDepositWithin50);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositCountWithin30", PatMu_allLayer1Muons_ecalIsoDepositCountWithin30, &b_PatMu_allLayer1Muons_ecalIsoDepositCountWithin30);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_ecalIsoDepositCountWithin50", PatMu_allLayer1Muons_ecalIsoDepositCountWithin50, &b_PatMu_allLayer1Muons_ecalIsoDepositCountWithin50);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isMuon", PatMu_allLayer1Muons_isMuon, &b_PatMu_allLayer1Muons_isMuon);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isCaloMuon", PatMu_allLayer1Muons_isCaloMuon, &b_PatMu_allLayer1Muons_isCaloMuon);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isTrackerMuon", PatMu_allLayer1Muons_isTrackerMuon, &b_PatMu_allLayer1Muons_isTrackerMuon);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isStandAloneMuon", PatMu_allLayer1Muons_isStandAloneMuon, &b_PatMu_allLayer1Muons_isStandAloneMuon);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isGlobalMuon", PatMu_allLayer1Muons_isGlobalMuon, &b_PatMu_allLayer1Muons_isGlobalMuon);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isGood", PatMu_allLayer1Muons_isGood, &b_PatMu_allLayer1Muons_isGood);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeAll", PatMu_allLayer1Muons_typeAll, &b_PatMu_allLayer1Muons_typeAll);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeAllGlobalMuons", PatMu_allLayer1Muons_typeAllGlobalMuons, &b_PatMu_allLayer1Muons_typeAllGlobalMuons);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeAllStandAloneMuons", PatMu_allLayer1Muons_typeAllStandAloneMuons, &b_PatMu_allLayer1Muons_typeAllStandAloneMuons);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeAllTrackerMuons", PatMu_allLayer1Muons_typeAllTrackerMuons, &b_PatMu_allLayer1Muons_typeAllTrackerMuons);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeTrackerMuonArbitrated", PatMu_allLayer1Muons_typeTrackerMuonArbitrated, &b_PatMu_allLayer1Muons_typeTrackerMuonArbitrated);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeAllArbitrated", PatMu_allLayer1Muons_typeAllArbitrated, &b_PatMu_allLayer1Muons_typeAllArbitrated);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeGlobalMuonPromptTight", PatMu_allLayer1Muons_typeGlobalMuonPromptTight, &b_PatMu_allLayer1Muons_typeGlobalMuonPromptTight);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeTMLastStationLoose", PatMu_allLayer1Muons_typeTMLastStationLoose, &b_PatMu_allLayer1Muons_typeTMLastStationLoose);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeTMLastStationTight", PatMu_allLayer1Muons_typeTMLastStationTight, &b_PatMu_allLayer1Muons_typeTMLastStationTight);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeTM2DCompatibilityLoose", PatMu_allLayer1Muons_typeTM2DCompatibilityLoose, &b_PatMu_allLayer1Muons_typeTM2DCompatibilityLoose);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeTM2DCompatibilityTight", PatMu_allLayer1Muons_typeTM2DCompatibilityTight, &b_PatMu_allLayer1Muons_typeTM2DCompatibilityTight);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeTMOneStationLoose", PatMu_allLayer1Muons_typeTMOneStationLoose, &b_PatMu_allLayer1Muons_typeTMOneStationLoose);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeTMOneStationTight", PatMu_allLayer1Muons_typeTMOneStationTight, &b_PatMu_allLayer1Muons_typeTMOneStationTight);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtLoose", PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtLoose, &b_PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtTight", PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtTight, &b_PatMu_allLayer1Muons_typeTMLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isIsolationValid", PatMu_allLayer1Muons_isIsolationValid, &b_PatMu_allLayer1Muons_isIsolationValid);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR03EmEt", PatMu_allLayer1Muons_isolationR03EmEt, &b_PatMu_allLayer1Muons_isolationR03EmEt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR03HadEt", PatMu_allLayer1Muons_isolationR03HadEt, &b_PatMu_allLayer1Muons_isolationR03HadEt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR03HoEt", PatMu_allLayer1Muons_isolationR03HoEt, &b_PatMu_allLayer1Muons_isolationR03HoEt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR03SumPt", PatMu_allLayer1Muons_isolationR03SumPt, &b_PatMu_allLayer1Muons_isolationR03SumPt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR03NTracks", PatMu_allLayer1Muons_isolationR03NTracks, &b_PatMu_allLayer1Muons_isolationR03NTracks);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR03NJets", PatMu_allLayer1Muons_isolationR03NJets, &b_PatMu_allLayer1Muons_isolationR03NJets);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR05EmEt", PatMu_allLayer1Muons_isolationR05EmEt, &b_PatMu_allLayer1Muons_isolationR05EmEt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR05HadEt", PatMu_allLayer1Muons_isolationR05HadEt, &b_PatMu_allLayer1Muons_isolationR05HadEt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR05HoEt", PatMu_allLayer1Muons_isolationR05HoEt, &b_PatMu_allLayer1Muons_isolationR05HoEt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR05SumPt", PatMu_allLayer1Muons_isolationR05SumPt, &b_PatMu_allLayer1Muons_isolationR05SumPt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR05NTracks", PatMu_allLayer1Muons_isolationR05NTracks, &b_PatMu_allLayer1Muons_isolationR05NTracks);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isolationR05NJets", PatMu_allLayer1Muons_isolationR05NJets, &b_PatMu_allLayer1Muons_isolationR05NJets);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isEnergyValid", PatMu_allLayer1Muons_isEnergyValid, &b_PatMu_allLayer1Muons_isEnergyValid);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_calEnergyTower", PatMu_allLayer1Muons_calEnergyTower, &b_PatMu_allLayer1Muons_calEnergyTower);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_calEnergyTowerS9", PatMu_allLayer1Muons_calEnergyTowerS9, &b_PatMu_allLayer1Muons_calEnergyTowerS9);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_calEnergyEm", PatMu_allLayer1Muons_calEnergyEm, &b_PatMu_allLayer1Muons_calEnergyEm);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_calEnergyEmS9", PatMu_allLayer1Muons_calEnergyEmS9, &b_PatMu_allLayer1Muons_calEnergyEmS9);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_calEnergyHad", PatMu_allLayer1Muons_calEnergyHad, &b_PatMu_allLayer1Muons_calEnergyHad);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_calEnergyHadS9", PatMu_allLayer1Muons_calEnergyHadS9, &b_PatMu_allLayer1Muons_calEnergyHadS9);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_calEnergyHo", PatMu_allLayer1Muons_calEnergyHo, &b_PatMu_allLayer1Muons_calEnergyHo);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_calEnergyHoS9", PatMu_allLayer1Muons_calEnergyHoS9, &b_PatMu_allLayer1Muons_calEnergyHoS9);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isCaloCompatibilityValid", PatMu_allLayer1Muons_isCaloCompatibilityValid, &b_PatMu_allLayer1Muons_isCaloCompatibilityValid);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_caloCompatibility", PatMu_allLayer1Muons_caloCompatibility, &b_PatMu_allLayer1Muons_caloCompatibility);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_isTimeValid", PatMu_allLayer1Muons_isTimeValid, &b_PatMu_allLayer1Muons_isTimeValid);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeDirection", PatMu_allLayer1Muons_timeDirection, &b_PatMu_allLayer1Muons_timeDirection);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeNStations", PatMu_allLayer1Muons_timeNStations, &b_PatMu_allLayer1Muons_timeNStations);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeInverseBeta", PatMu_allLayer1Muons_timeInverseBeta, &b_PatMu_allLayer1Muons_timeInverseBeta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeInverseBetaErr", PatMu_allLayer1Muons_timeInverseBetaErr, &b_PatMu_allLayer1Muons_timeInverseBetaErr);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeFreeInverseBeta", PatMu_allLayer1Muons_timeFreeInverseBeta, &b_PatMu_allLayer1Muons_timeFreeInverseBeta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeFreeInverseBetaErr", PatMu_allLayer1Muons_timeFreeInverseBetaErr, &b_PatMu_allLayer1Muons_timeFreeInverseBetaErr);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeTimeAtIpInOut", PatMu_allLayer1Muons_timeTimeAtIpInOut, &b_PatMu_allLayer1Muons_timeTimeAtIpInOut);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeTimeAtIpInOutErr", PatMu_allLayer1Muons_timeTimeAtIpInOutErr, &b_PatMu_allLayer1Muons_timeTimeAtIpInOutErr);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeTimeAtIpOutIn", PatMu_allLayer1Muons_timeTimeAtIpOutIn, &b_PatMu_allLayer1Muons_timeTimeAtIpOutIn);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_timeTimeAtIpOutInErr", PatMu_allLayer1Muons_timeTimeAtIpOutInErr, &b_PatMu_allLayer1Muons_timeTimeAtIpOutInErr);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_charge", PatMu_allLayer1Muons_charge, &b_PatMu_allLayer1Muons_charge);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_px", PatMu_allLayer1Muons_px, &b_PatMu_allLayer1Muons_px);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_py", PatMu_allLayer1Muons_py, &b_PatMu_allLayer1Muons_py);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_pz", PatMu_allLayer1Muons_pz, &b_PatMu_allLayer1Muons_pz);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_pt", PatMu_allLayer1Muons_pt, &b_PatMu_allLayer1Muons_pt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_p", PatMu_allLayer1Muons_p, &b_PatMu_allLayer1Muons_p);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_energy", PatMu_allLayer1Muons_energy, &b_PatMu_allLayer1Muons_energy);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_phi", PatMu_allLayer1Muons_phi, &b_PatMu_allLayer1Muons_phi);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_eta", PatMu_allLayer1Muons_eta, &b_PatMu_allLayer1Muons_eta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_theta", PatMu_allLayer1Muons_theta, &b_PatMu_allLayer1Muons_theta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_vx", PatMu_allLayer1Muons_vx, &b_PatMu_allLayer1Muons_vx);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_vy", PatMu_allLayer1Muons_vy, &b_PatMu_allLayer1Muons_vy);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_vz", PatMu_allLayer1Muons_vz, &b_PatMu_allLayer1Muons_vz);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_globalTrackNormalizedChi2", PatMu_allLayer1Muons_globalTrackNormalizedChi2, &b_PatMu_allLayer1Muons_globalTrackNormalizedChi2);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_innerTrackNumberOfValidHits", PatMu_allLayer1Muons_innerTrackNumberOfValidHits, &b_PatMu_allLayer1Muons_innerTrackNumberOfValidHits);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_innerTrackNumberOfLostHits", PatMu_allLayer1Muons_innerTrackNumberOfLostHits, &b_PatMu_allLayer1Muons_innerTrackNumberOfLostHits);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_innerTrackPhi", PatMu_allLayer1Muons_innerTrackPhi, &b_PatMu_allLayer1Muons_innerTrackPhi);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_innerTrackD0", PatMu_allLayer1Muons_innerTrackD0, &b_PatMu_allLayer1Muons_innerTrackD0);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_innerTrackD0BeamSpotCorrected", PatMu_allLayer1Muons_innerTrackD0BeamSpotCorrected, &b_PatMu_allLayer1Muons_innerTrackD0BeamSpotCorrected);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_nPatMuonTrack", &PatMu_allLayer1Muons_nPatMuonTrack, &b_PatMu_allLayer1Muons_nPatMuonTrack);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackParentPatMuonIndex", PatMu_allLayer1Muons_trackParentPatMuonIndex, &b_PatMu_allLayer1Muons_trackParentPatMuonIndex);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackType", PatMu_allLayer1Muons_trackType, &b_PatMu_allLayer1Muons_trackType);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackCharge", PatMu_allLayer1Muons_trackCharge, &b_PatMu_allLayer1Muons_trackCharge);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackPx", PatMu_allLayer1Muons_trackPx, &b_PatMu_allLayer1Muons_trackPx);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackPy", PatMu_allLayer1Muons_trackPy, &b_PatMu_allLayer1Muons_trackPy);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackPz", PatMu_allLayer1Muons_trackPz, &b_PatMu_allLayer1Muons_trackPz);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackPt", PatMu_allLayer1Muons_trackPt, &b_PatMu_allLayer1Muons_trackPt);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackPtError", PatMu_allLayer1Muons_trackPtError, &b_PatMu_allLayer1Muons_trackPtError);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackP", PatMu_allLayer1Muons_trackP, &b_PatMu_allLayer1Muons_trackP);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackPhi", PatMu_allLayer1Muons_trackPhi, &b_PatMu_allLayer1Muons_trackPhi);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackPhiError", PatMu_allLayer1Muons_trackPhiError, &b_PatMu_allLayer1Muons_trackPhiError);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackEta", PatMu_allLayer1Muons_trackEta, &b_PatMu_allLayer1Muons_trackEta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackEtaError", PatMu_allLayer1Muons_trackEtaError, &b_PatMu_allLayer1Muons_trackEtaError);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackTheta", PatMu_allLayer1Muons_trackTheta, &b_PatMu_allLayer1Muons_trackTheta);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackThetaError", PatMu_allLayer1Muons_trackThetaError, &b_PatMu_allLayer1Muons_trackThetaError);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackVx", PatMu_allLayer1Muons_trackVx, &b_PatMu_allLayer1Muons_trackVx);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackVy", PatMu_allLayer1Muons_trackVy, &b_PatMu_allLayer1Muons_trackVy);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackVz", PatMu_allLayer1Muons_trackVz, &b_PatMu_allLayer1Muons_trackVz);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackChi2", PatMu_allLayer1Muons_trackChi2, &b_PatMu_allLayer1Muons_trackChi2);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackChi2Norm", PatMu_allLayer1Muons_trackChi2Norm, &b_PatMu_allLayer1Muons_trackChi2Norm);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackNdof", PatMu_allLayer1Muons_trackNdof, &b_PatMu_allLayer1Muons_trackNdof);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackD0", PatMu_allLayer1Muons_trackD0, &b_PatMu_allLayer1Muons_trackD0);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackD0Error", PatMu_allLayer1Muons_trackD0Error, &b_PatMu_allLayer1Muons_trackD0Error);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackDsz", PatMu_allLayer1Muons_trackDsz, &b_PatMu_allLayer1Muons_trackDsz);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackDszError", PatMu_allLayer1Muons_trackDszError, &b_PatMu_allLayer1Muons_trackDszError);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackDxy", PatMu_allLayer1Muons_trackDxy, &b_PatMu_allLayer1Muons_trackDxy);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackDxyError", PatMu_allLayer1Muons_trackDxyError, &b_PatMu_allLayer1Muons_trackDxyError);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackDz", PatMu_allLayer1Muons_trackDz, &b_PatMu_allLayer1Muons_trackDz);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackDzError", PatMu_allLayer1Muons_trackDzError, &b_PatMu_allLayer1Muons_trackDzError);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackNumberOfValidHits", PatMu_allLayer1Muons_trackNumberOfValidHits, &b_PatMu_allLayer1Muons_trackNumberOfValidHits);
   fChain->SetBranchAddress("PatMu_allLayer1Muons_trackNumberOfLostHits", PatMu_allLayer1Muons_trackNumberOfLostHits, &b_PatMu_allLayer1Muons_trackNumberOfLostHits);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_nPatElectron", &PatEl_allLayer1Electrons_nPatElectron, &b_PatEl_allLayer1Electrons_nPatElectron);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_eidRobustLoose", PatEl_allLayer1Electrons_eidRobustLoose, &b_PatEl_allLayer1Electrons_eidRobustLoose);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_eidRobustTight", PatEl_allLayer1Electrons_eidRobustTight, &b_PatEl_allLayer1Electrons_eidRobustTight);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_eidRobustHighEnergy", PatEl_allLayer1Electrons_eidRobustHighEnergy, &b_PatEl_allLayer1Electrons_eidRobustHighEnergy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_eidLoose", PatEl_allLayer1Electrons_eidLoose, &b_PatEl_allLayer1Electrons_eidLoose);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_eidTight", PatEl_allLayer1Electrons_eidTight, &b_PatEl_allLayer1Electrons_eidTight);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_scSigmaEtaEta", PatEl_allLayer1Electrons_scSigmaEtaEta, &b_PatEl_allLayer1Electrons_scSigmaEtaEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_scSigmaIEtaIEta", PatEl_allLayer1Electrons_scSigmaIEtaIEta, &b_PatEl_allLayer1Electrons_scSigmaIEtaIEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_scE1x5", PatEl_allLayer1Electrons_scE1x5, &b_PatEl_allLayer1Electrons_scE1x5);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_scE2x5Max", PatEl_allLayer1Electrons_scE2x5Max, &b_PatEl_allLayer1Electrons_scE2x5Max);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_scE5x5", PatEl_allLayer1Electrons_scE5x5, &b_PatEl_allLayer1Electrons_scE5x5);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackIso", PatEl_allLayer1Electrons_trackIso, &b_PatEl_allLayer1Electrons_trackIso);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIso", PatEl_allLayer1Electrons_hcalIso, &b_PatEl_allLayer1Electrons_hcalIso);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIso", PatEl_allLayer1Electrons_ecalIso, &b_PatEl_allLayer1Electrons_ecalIso);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_caloIso", PatEl_allLayer1Electrons_caloIso, &b_PatEl_allLayer1Electrons_caloIso);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositEta", PatEl_allLayer1Electrons_trackerIsoDepositEta, &b_PatEl_allLayer1Electrons_trackerIsoDepositEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositPhi", PatEl_allLayer1Electrons_trackerIsoDepositPhi, &b_PatEl_allLayer1Electrons_trackerIsoDepositPhi);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositVetoEta", PatEl_allLayer1Electrons_trackerIsoDepositVetoEta, &b_PatEl_allLayer1Electrons_trackerIsoDepositVetoEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositVetoPhi", PatEl_allLayer1Electrons_trackerIsoDepositVetoPhi, &b_PatEl_allLayer1Electrons_trackerIsoDepositVetoPhi);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositVetoDR", PatEl_allLayer1Electrons_trackerIsoDepositVetoDR, &b_PatEl_allLayer1Electrons_trackerIsoDepositVetoDR);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositCandEnergy", PatEl_allLayer1Electrons_trackerIsoDepositCandEnergy, &b_PatEl_allLayer1Electrons_trackerIsoDepositCandEnergy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositWithin30", PatEl_allLayer1Electrons_trackerIsoDepositWithin30, &b_PatEl_allLayer1Electrons_trackerIsoDepositWithin30);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositWithin50", PatEl_allLayer1Electrons_trackerIsoDepositWithin50, &b_PatEl_allLayer1Electrons_trackerIsoDepositWithin50);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositCountWithin30", PatEl_allLayer1Electrons_trackerIsoDepositCountWithin30, &b_PatEl_allLayer1Electrons_trackerIsoDepositCountWithin30);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackerIsoDepositCountWithin50", PatEl_allLayer1Electrons_trackerIsoDepositCountWithin50, &b_PatEl_allLayer1Electrons_trackerIsoDepositCountWithin50);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositEta", PatEl_allLayer1Electrons_hcalIsoDepositEta, &b_PatEl_allLayer1Electrons_hcalIsoDepositEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositPhi", PatEl_allLayer1Electrons_hcalIsoDepositPhi, &b_PatEl_allLayer1Electrons_hcalIsoDepositPhi);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositVetoEta", PatEl_allLayer1Electrons_hcalIsoDepositVetoEta, &b_PatEl_allLayer1Electrons_hcalIsoDepositVetoEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositVetoPhi", PatEl_allLayer1Electrons_hcalIsoDepositVetoPhi, &b_PatEl_allLayer1Electrons_hcalIsoDepositVetoPhi);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositVetoDR", PatEl_allLayer1Electrons_hcalIsoDepositVetoDR, &b_PatEl_allLayer1Electrons_hcalIsoDepositVetoDR);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositCandEnergy", PatEl_allLayer1Electrons_hcalIsoDepositCandEnergy, &b_PatEl_allLayer1Electrons_hcalIsoDepositCandEnergy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositWithin30", PatEl_allLayer1Electrons_hcalIsoDepositWithin30, &b_PatEl_allLayer1Electrons_hcalIsoDepositWithin30);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositWithin50", PatEl_allLayer1Electrons_hcalIsoDepositWithin50, &b_PatEl_allLayer1Electrons_hcalIsoDepositWithin50);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositCountWithin30", PatEl_allLayer1Electrons_hcalIsoDepositCountWithin30, &b_PatEl_allLayer1Electrons_hcalIsoDepositCountWithin30);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hcalIsoDepositCountWithin50", PatEl_allLayer1Electrons_hcalIsoDepositCountWithin50, &b_PatEl_allLayer1Electrons_hcalIsoDepositCountWithin50);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositEta", PatEl_allLayer1Electrons_ecalIsoDepositEta, &b_PatEl_allLayer1Electrons_ecalIsoDepositEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositPhi", PatEl_allLayer1Electrons_ecalIsoDepositPhi, &b_PatEl_allLayer1Electrons_ecalIsoDepositPhi);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositVetoEta", PatEl_allLayer1Electrons_ecalIsoDepositVetoEta, &b_PatEl_allLayer1Electrons_ecalIsoDepositVetoEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositVetoPhi", PatEl_allLayer1Electrons_ecalIsoDepositVetoPhi, &b_PatEl_allLayer1Electrons_ecalIsoDepositVetoPhi);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositVetoDR", PatEl_allLayer1Electrons_ecalIsoDepositVetoDR, &b_PatEl_allLayer1Electrons_ecalIsoDepositVetoDR);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositCandEnergy", PatEl_allLayer1Electrons_ecalIsoDepositCandEnergy, &b_PatEl_allLayer1Electrons_ecalIsoDepositCandEnergy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositWithin30", PatEl_allLayer1Electrons_ecalIsoDepositWithin30, &b_PatEl_allLayer1Electrons_ecalIsoDepositWithin30);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositWithin50", PatEl_allLayer1Electrons_ecalIsoDepositWithin50, &b_PatEl_allLayer1Electrons_ecalIsoDepositWithin50);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositCountWithin30", PatEl_allLayer1Electrons_ecalIsoDepositCountWithin30, &b_PatEl_allLayer1Electrons_ecalIsoDepositCountWithin30);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_ecalIsoDepositCountWithin50", PatEl_allLayer1Electrons_ecalIsoDepositCountWithin50, &b_PatEl_allLayer1Electrons_ecalIsoDepositCountWithin50);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_isElectron", PatEl_allLayer1Electrons_isElectron, &b_PatEl_allLayer1Electrons_isElectron);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_classification", PatEl_allLayer1Electrons_classification, &b_PatEl_allLayer1Electrons_classification);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_caloEnergy", PatEl_allLayer1Electrons_caloEnergy, &b_PatEl_allLayer1Electrons_caloEnergy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_caloEnergyError", PatEl_allLayer1Electrons_caloEnergyError, &b_PatEl_allLayer1Electrons_caloEnergyError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_isEnergyScaleCorrected", PatEl_allLayer1Electrons_isEnergyScaleCorrected, &b_PatEl_allLayer1Electrons_isEnergyScaleCorrected);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumError", PatEl_allLayer1Electrons_trackMomentumError, &b_PatEl_allLayer1Electrons_trackMomentumError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_isMomentumCorrected", PatEl_allLayer1Electrons_isMomentumCorrected, &b_PatEl_allLayer1Electrons_isMomentumCorrected);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_hadronicOverEm", PatEl_allLayer1Electrons_hadronicOverEm, &b_PatEl_allLayer1Electrons_hadronicOverEm);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_caloPositionX", PatEl_allLayer1Electrons_caloPositionX, &b_PatEl_allLayer1Electrons_caloPositionX);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_caloPositionY", PatEl_allLayer1Electrons_caloPositionY, &b_PatEl_allLayer1Electrons_caloPositionY);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_caloPositionZ", PatEl_allLayer1Electrons_caloPositionZ, &b_PatEl_allLayer1Electrons_caloPositionZ);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumAtVtxPx", PatEl_allLayer1Electrons_trackMomentumAtVtxPx, &b_PatEl_allLayer1Electrons_trackMomentumAtVtxPx);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumAtVtxPy", PatEl_allLayer1Electrons_trackMomentumAtVtxPy, &b_PatEl_allLayer1Electrons_trackMomentumAtVtxPy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumAtVtxPz", PatEl_allLayer1Electrons_trackMomentumAtVtxPz, &b_PatEl_allLayer1Electrons_trackMomentumAtVtxPz);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPositionAtVtxX", PatEl_allLayer1Electrons_trackPositionAtVtxX, &b_PatEl_allLayer1Electrons_trackPositionAtVtxX);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPositionAtVtxY", PatEl_allLayer1Electrons_trackPositionAtVtxY, &b_PatEl_allLayer1Electrons_trackPositionAtVtxY);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPositionAtVtxZ", PatEl_allLayer1Electrons_trackPositionAtVtxZ, &b_PatEl_allLayer1Electrons_trackPositionAtVtxZ);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumAtCaloPx", PatEl_allLayer1Electrons_trackMomentumAtCaloPx, &b_PatEl_allLayer1Electrons_trackMomentumAtCaloPx);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumAtCaloPy", PatEl_allLayer1Electrons_trackMomentumAtCaloPy, &b_PatEl_allLayer1Electrons_trackMomentumAtCaloPy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumAtCaloPz", PatEl_allLayer1Electrons_trackMomentumAtCaloPz, &b_PatEl_allLayer1Electrons_trackMomentumAtCaloPz);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPositionAtCaloX", PatEl_allLayer1Electrons_trackPositionAtCaloX, &b_PatEl_allLayer1Electrons_trackPositionAtCaloX);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPositionAtCaloY", PatEl_allLayer1Electrons_trackPositionAtCaloY, &b_PatEl_allLayer1Electrons_trackPositionAtCaloY);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPositionAtCaloZ", PatEl_allLayer1Electrons_trackPositionAtCaloZ, &b_PatEl_allLayer1Electrons_trackPositionAtCaloZ);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumOutPx", PatEl_allLayer1Electrons_trackMomentumOutPx, &b_PatEl_allLayer1Electrons_trackMomentumOutPx);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumOutPy", PatEl_allLayer1Electrons_trackMomentumOutPy, &b_PatEl_allLayer1Electrons_trackMomentumOutPy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackMomentumOutPz", PatEl_allLayer1Electrons_trackMomentumOutPz, &b_PatEl_allLayer1Electrons_trackMomentumOutPz);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_eSuperClusterOverP", PatEl_allLayer1Electrons_eSuperClusterOverP, &b_PatEl_allLayer1Electrons_eSuperClusterOverP);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_eSeedClusterOverPout", PatEl_allLayer1Electrons_eSeedClusterOverPout, &b_PatEl_allLayer1Electrons_eSeedClusterOverPout);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_deltaEtaSuperClusterTrackAtVtx", PatEl_allLayer1Electrons_deltaEtaSuperClusterTrackAtVtx, &b_PatEl_allLayer1Electrons_deltaEtaSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_deltaEtaSeedClusterTrackAtCalo", PatEl_allLayer1Electrons_deltaEtaSeedClusterTrackAtCalo, &b_PatEl_allLayer1Electrons_deltaEtaSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_deltaPhiSuperClusterTrackAtVtx", PatEl_allLayer1Electrons_deltaPhiSuperClusterTrackAtVtx, &b_PatEl_allLayer1Electrons_deltaPhiSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_deltaPhiSeedClusterTrackAtCalo", PatEl_allLayer1Electrons_deltaPhiSeedClusterTrackAtCalo, &b_PatEl_allLayer1Electrons_deltaPhiSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackChargeMode", PatEl_allLayer1Electrons_gsfTrackChargeMode, &b_PatEl_allLayer1Electrons_gsfTrackChargeMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackQoverpMode", PatEl_allLayer1Electrons_gsfTrackQoverpMode, &b_PatEl_allLayer1Electrons_gsfTrackQoverpMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackQoverpModeError", PatEl_allLayer1Electrons_gsfTrackQoverpModeError, &b_PatEl_allLayer1Electrons_gsfTrackQoverpModeError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackPxMode", PatEl_allLayer1Electrons_gsfTrackPxMode, &b_PatEl_allLayer1Electrons_gsfTrackPxMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackPyMode", PatEl_allLayer1Electrons_gsfTrackPyMode, &b_PatEl_allLayer1Electrons_gsfTrackPyMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackPzMode", PatEl_allLayer1Electrons_gsfTrackPzMode, &b_PatEl_allLayer1Electrons_gsfTrackPzMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackPtMode", PatEl_allLayer1Electrons_gsfTrackPtMode, &b_PatEl_allLayer1Electrons_gsfTrackPtMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackPtModeError", PatEl_allLayer1Electrons_gsfTrackPtModeError, &b_PatEl_allLayer1Electrons_gsfTrackPtModeError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackPMode", PatEl_allLayer1Electrons_gsfTrackPMode, &b_PatEl_allLayer1Electrons_gsfTrackPMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackPhiMode", PatEl_allLayer1Electrons_gsfTrackPhiMode, &b_PatEl_allLayer1Electrons_gsfTrackPhiMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackPhiModeError", PatEl_allLayer1Electrons_gsfTrackPhiModeError, &b_PatEl_allLayer1Electrons_gsfTrackPhiModeError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackEtaMode", PatEl_allLayer1Electrons_gsfTrackEtaMode, &b_PatEl_allLayer1Electrons_gsfTrackEtaMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackEtaModeError", PatEl_allLayer1Electrons_gsfTrackEtaModeError, &b_PatEl_allLayer1Electrons_gsfTrackEtaModeError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackThetaMode", PatEl_allLayer1Electrons_gsfTrackThetaMode, &b_PatEl_allLayer1Electrons_gsfTrackThetaMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackThetaModeError", PatEl_allLayer1Electrons_gsfTrackThetaModeError, &b_PatEl_allLayer1Electrons_gsfTrackThetaModeError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackLambdaMode", PatEl_allLayer1Electrons_gsfTrackLambdaMode, &b_PatEl_allLayer1Electrons_gsfTrackLambdaMode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackLambdaModeError", PatEl_allLayer1Electrons_gsfTrackLambdaModeError, &b_PatEl_allLayer1Electrons_gsfTrackLambdaModeError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackParameter0Mode", PatEl_allLayer1Electrons_gsfTrackParameter0Mode, &b_PatEl_allLayer1Electrons_gsfTrackParameter0Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackParameter1Mode", PatEl_allLayer1Electrons_gsfTrackParameter1Mode, &b_PatEl_allLayer1Electrons_gsfTrackParameter1Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackParameter2Mode", PatEl_allLayer1Electrons_gsfTrackParameter2Mode, &b_PatEl_allLayer1Electrons_gsfTrackParameter2Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix00Mode", PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix00Mode, &b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix00Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix01Mode", PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix01Mode, &b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix01Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix02Mode", PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix02Mode, &b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix02Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix10Mode", PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix10Mode, &b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix10Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix11Mode", PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix11Mode, &b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix11Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix12Mode", PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix12Mode, &b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix12Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix20Mode", PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix20Mode, &b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix20Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix21Mode", PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix21Mode, &b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix21Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix22Mode", PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix22Mode, &b_PatEl_allLayer1Electrons_gsfTrackCovarianceMatrix22Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackError00Mode", PatEl_allLayer1Electrons_gsfTrackError00Mode, &b_PatEl_allLayer1Electrons_gsfTrackError00Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackError11Mode", PatEl_allLayer1Electrons_gsfTrackError11Mode, &b_PatEl_allLayer1Electrons_gsfTrackError11Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_gsfTrackError22Mode", PatEl_allLayer1Electrons_gsfTrackError22Mode, &b_PatEl_allLayer1Electrons_gsfTrackError22Mode);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackCharge", PatEl_allLayer1Electrons_trackCharge, &b_PatEl_allLayer1Electrons_trackCharge);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPx", PatEl_allLayer1Electrons_trackPx, &b_PatEl_allLayer1Electrons_trackPx);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPy", PatEl_allLayer1Electrons_trackPy, &b_PatEl_allLayer1Electrons_trackPy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPz", PatEl_allLayer1Electrons_trackPz, &b_PatEl_allLayer1Electrons_trackPz);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPt", PatEl_allLayer1Electrons_trackPt, &b_PatEl_allLayer1Electrons_trackPt);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPtError", PatEl_allLayer1Electrons_trackPtError, &b_PatEl_allLayer1Electrons_trackPtError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackP", PatEl_allLayer1Electrons_trackP, &b_PatEl_allLayer1Electrons_trackP);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPhi", PatEl_allLayer1Electrons_trackPhi, &b_PatEl_allLayer1Electrons_trackPhi);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackPhiError", PatEl_allLayer1Electrons_trackPhiError, &b_PatEl_allLayer1Electrons_trackPhiError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackEta", PatEl_allLayer1Electrons_trackEta, &b_PatEl_allLayer1Electrons_trackEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackEtaError", PatEl_allLayer1Electrons_trackEtaError, &b_PatEl_allLayer1Electrons_trackEtaError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackTheta", PatEl_allLayer1Electrons_trackTheta, &b_PatEl_allLayer1Electrons_trackTheta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackThetaError", PatEl_allLayer1Electrons_trackThetaError, &b_PatEl_allLayer1Electrons_trackThetaError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackVx", PatEl_allLayer1Electrons_trackVx, &b_PatEl_allLayer1Electrons_trackVx);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackVy", PatEl_allLayer1Electrons_trackVy, &b_PatEl_allLayer1Electrons_trackVy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackVz", PatEl_allLayer1Electrons_trackVz, &b_PatEl_allLayer1Electrons_trackVz);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackChi2", PatEl_allLayer1Electrons_trackChi2, &b_PatEl_allLayer1Electrons_trackChi2);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackChi2Norm", PatEl_allLayer1Electrons_trackChi2Norm, &b_PatEl_allLayer1Electrons_trackChi2Norm);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackNdof", PatEl_allLayer1Electrons_trackNdof, &b_PatEl_allLayer1Electrons_trackNdof);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackD0", PatEl_allLayer1Electrons_trackD0, &b_PatEl_allLayer1Electrons_trackD0);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackD0BeamSpotCorrected", PatEl_allLayer1Electrons_trackD0BeamSpotCorrected, &b_PatEl_allLayer1Electrons_trackD0BeamSpotCorrected);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackD0Error", PatEl_allLayer1Electrons_trackD0Error, &b_PatEl_allLayer1Electrons_trackD0Error);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackDsz", PatEl_allLayer1Electrons_trackDsz, &b_PatEl_allLayer1Electrons_trackDsz);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackDszError", PatEl_allLayer1Electrons_trackDszError, &b_PatEl_allLayer1Electrons_trackDszError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackDxy", PatEl_allLayer1Electrons_trackDxy, &b_PatEl_allLayer1Electrons_trackDxy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackDxyError", PatEl_allLayer1Electrons_trackDxyError, &b_PatEl_allLayer1Electrons_trackDxyError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackDz", PatEl_allLayer1Electrons_trackDz, &b_PatEl_allLayer1Electrons_trackDz);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackDzError", PatEl_allLayer1Electrons_trackDzError, &b_PatEl_allLayer1Electrons_trackDzError);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackNumberOfValidHits", PatEl_allLayer1Electrons_trackNumberOfValidHits, &b_PatEl_allLayer1Electrons_trackNumberOfValidHits);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_trackNumberOfLostHits", PatEl_allLayer1Electrons_trackNumberOfLostHits, &b_PatEl_allLayer1Electrons_trackNumberOfLostHits);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterEnergy", PatEl_allLayer1Electrons_superClusterEnergy, &b_PatEl_allLayer1Electrons_superClusterEnergy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterX", PatEl_allLayer1Electrons_superClusterX, &b_PatEl_allLayer1Electrons_superClusterX);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterY", PatEl_allLayer1Electrons_superClusterY, &b_PatEl_allLayer1Electrons_superClusterY);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterZ", PatEl_allLayer1Electrons_superClusterZ, &b_PatEl_allLayer1Electrons_superClusterZ);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterEta", PatEl_allLayer1Electrons_superClusterEta, &b_PatEl_allLayer1Electrons_superClusterEta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterPhi", PatEl_allLayer1Electrons_superClusterPhi, &b_PatEl_allLayer1Electrons_superClusterPhi);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterCaloID", PatEl_allLayer1Electrons_superClusterCaloID, &b_PatEl_allLayer1Electrons_superClusterCaloID);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterPreshowerEnergy", PatEl_allLayer1Electrons_superClusterPreshowerEnergy, &b_PatEl_allLayer1Electrons_superClusterPreshowerEnergy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterPhiWidth", PatEl_allLayer1Electrons_superClusterPhiWidth, &b_PatEl_allLayer1Electrons_superClusterPhiWidth);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_superClusterEtaWidth", PatEl_allLayer1Electrons_superClusterEtaWidth, &b_PatEl_allLayer1Electrons_superClusterEtaWidth);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_charge", PatEl_allLayer1Electrons_charge, &b_PatEl_allLayer1Electrons_charge);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_px", PatEl_allLayer1Electrons_px, &b_PatEl_allLayer1Electrons_px);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_py", PatEl_allLayer1Electrons_py, &b_PatEl_allLayer1Electrons_py);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_pz", PatEl_allLayer1Electrons_pz, &b_PatEl_allLayer1Electrons_pz);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_pt", PatEl_allLayer1Electrons_pt, &b_PatEl_allLayer1Electrons_pt);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_p", PatEl_allLayer1Electrons_p, &b_PatEl_allLayer1Electrons_p);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_energy", PatEl_allLayer1Electrons_energy, &b_PatEl_allLayer1Electrons_energy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_phi", PatEl_allLayer1Electrons_phi, &b_PatEl_allLayer1Electrons_phi);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_eta", PatEl_allLayer1Electrons_eta, &b_PatEl_allLayer1Electrons_eta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_theta", PatEl_allLayer1Electrons_theta, &b_PatEl_allLayer1Electrons_theta);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_vx", PatEl_allLayer1Electrons_vx, &b_PatEl_allLayer1Electrons_vx);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_vy", PatEl_allLayer1Electrons_vy, &b_PatEl_allLayer1Electrons_vy);
   fChain->SetBranchAddress("PatEl_allLayer1Electrons_vz", PatEl_allLayer1Electrons_vz, &b_PatEl_allLayer1Electrons_vz);
   fChain->SetBranchAddress("Track_generalTracks_nRecoTrack", &Track_generalTracks_nRecoTrack, &b_Track_generalTracks_nRecoTrack);
   fChain->SetBranchAddress("Track_generalTracks_trackCharge", Track_generalTracks_trackCharge, &b_Track_generalTracks_trackCharge);
   fChain->SetBranchAddress("Track_generalTracks_trackPx", Track_generalTracks_trackPx, &b_Track_generalTracks_trackPx);
   fChain->SetBranchAddress("Track_generalTracks_trackPy", Track_generalTracks_trackPy, &b_Track_generalTracks_trackPy);
   fChain->SetBranchAddress("Track_generalTracks_trackPz", Track_generalTracks_trackPz, &b_Track_generalTracks_trackPz);
   fChain->SetBranchAddress("Track_generalTracks_trackPt", Track_generalTracks_trackPt, &b_Track_generalTracks_trackPt);
   fChain->SetBranchAddress("Track_generalTracks_trackPtError", Track_generalTracks_trackPtError, &b_Track_generalTracks_trackPtError);
   fChain->SetBranchAddress("Track_generalTracks_trackP", Track_generalTracks_trackP, &b_Track_generalTracks_trackP);
   fChain->SetBranchAddress("Track_generalTracks_trackPhi", Track_generalTracks_trackPhi, &b_Track_generalTracks_trackPhi);
   fChain->SetBranchAddress("Track_generalTracks_trackPhiError", Track_generalTracks_trackPhiError, &b_Track_generalTracks_trackPhiError);
   fChain->SetBranchAddress("Track_generalTracks_trackEta", Track_generalTracks_trackEta, &b_Track_generalTracks_trackEta);
   fChain->SetBranchAddress("Track_generalTracks_trackEtaError", Track_generalTracks_trackEtaError, &b_Track_generalTracks_trackEtaError);
   fChain->SetBranchAddress("Track_generalTracks_trackTheta", Track_generalTracks_trackTheta, &b_Track_generalTracks_trackTheta);
   fChain->SetBranchAddress("Track_generalTracks_trackThetaError", Track_generalTracks_trackThetaError, &b_Track_generalTracks_trackThetaError);
   fChain->SetBranchAddress("Track_generalTracks_trackVx", Track_generalTracks_trackVx, &b_Track_generalTracks_trackVx);
   fChain->SetBranchAddress("Track_generalTracks_trackVy", Track_generalTracks_trackVy, &b_Track_generalTracks_trackVy);
   fChain->SetBranchAddress("Track_generalTracks_trackVz", Track_generalTracks_trackVz, &b_Track_generalTracks_trackVz);
   fChain->SetBranchAddress("Track_generalTracks_trackChi2", Track_generalTracks_trackChi2, &b_Track_generalTracks_trackChi2);
   fChain->SetBranchAddress("Track_generalTracks_trackChi2Norm", Track_generalTracks_trackChi2Norm, &b_Track_generalTracks_trackChi2Norm);
   fChain->SetBranchAddress("Track_generalTracks_trackNdof", Track_generalTracks_trackNdof, &b_Track_generalTracks_trackNdof);
   fChain->SetBranchAddress("Track_generalTracks_trackD0", Track_generalTracks_trackD0, &b_Track_generalTracks_trackD0);
   fChain->SetBranchAddress("Track_generalTracks_trackD0Error", Track_generalTracks_trackD0Error, &b_Track_generalTracks_trackD0Error);
   fChain->SetBranchAddress("Track_generalTracks_trackDsz", Track_generalTracks_trackDsz, &b_Track_generalTracks_trackDsz);
   fChain->SetBranchAddress("Track_generalTracks_trackDszError", Track_generalTracks_trackDszError, &b_Track_generalTracks_trackDszError);
   fChain->SetBranchAddress("Track_generalTracks_trackDxy", Track_generalTracks_trackDxy, &b_Track_generalTracks_trackDxy);
   fChain->SetBranchAddress("Track_generalTracks_trackDxyError", Track_generalTracks_trackDxyError, &b_Track_generalTracks_trackDxyError);
   fChain->SetBranchAddress("Track_generalTracks_trackDz", Track_generalTracks_trackDz, &b_Track_generalTracks_trackDz);
   fChain->SetBranchAddress("Track_generalTracks_trackDzError", Track_generalTracks_trackDzError, &b_Track_generalTracks_trackDzError);
   fChain->SetBranchAddress("Track_generalTracks_trackNumberOfValidHits", Track_generalTracks_trackNumberOfValidHits, &b_Track_generalTracks_trackNumberOfValidHits);
   fChain->SetBranchAddress("Track_generalTracks_trackNumberOfLostHits", Track_generalTracks_trackNumberOfLostHits, &b_Track_generalTracks_trackNumberOfLostHits);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_nRecoTrack", &Track_ckfInOutTracksFromConversions_nRecoTrack, &b_Track_ckfInOutTracksFromConversions_nRecoTrack);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackCharge", Track_ckfInOutTracksFromConversions_trackCharge, &b_Track_ckfInOutTracksFromConversions_trackCharge);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackPx", Track_ckfInOutTracksFromConversions_trackPx, &b_Track_ckfInOutTracksFromConversions_trackPx);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackPy", Track_ckfInOutTracksFromConversions_trackPy, &b_Track_ckfInOutTracksFromConversions_trackPy);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackPz", Track_ckfInOutTracksFromConversions_trackPz, &b_Track_ckfInOutTracksFromConversions_trackPz);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackPt", Track_ckfInOutTracksFromConversions_trackPt, &b_Track_ckfInOutTracksFromConversions_trackPt);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackPtError", Track_ckfInOutTracksFromConversions_trackPtError, &b_Track_ckfInOutTracksFromConversions_trackPtError);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackP", Track_ckfInOutTracksFromConversions_trackP, &b_Track_ckfInOutTracksFromConversions_trackP);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackPhi", Track_ckfInOutTracksFromConversions_trackPhi, &b_Track_ckfInOutTracksFromConversions_trackPhi);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackPhiError", Track_ckfInOutTracksFromConversions_trackPhiError, &b_Track_ckfInOutTracksFromConversions_trackPhiError);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackEta", Track_ckfInOutTracksFromConversions_trackEta, &b_Track_ckfInOutTracksFromConversions_trackEta);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackEtaError", Track_ckfInOutTracksFromConversions_trackEtaError, &b_Track_ckfInOutTracksFromConversions_trackEtaError);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackTheta", Track_ckfInOutTracksFromConversions_trackTheta, &b_Track_ckfInOutTracksFromConversions_trackTheta);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackThetaError", Track_ckfInOutTracksFromConversions_trackThetaError, &b_Track_ckfInOutTracksFromConversions_trackThetaError);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackVx", Track_ckfInOutTracksFromConversions_trackVx, &b_Track_ckfInOutTracksFromConversions_trackVx);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackVy", Track_ckfInOutTracksFromConversions_trackVy, &b_Track_ckfInOutTracksFromConversions_trackVy);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackVz", Track_ckfInOutTracksFromConversions_trackVz, &b_Track_ckfInOutTracksFromConversions_trackVz);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackChi2", Track_ckfInOutTracksFromConversions_trackChi2, &b_Track_ckfInOutTracksFromConversions_trackChi2);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackChi2Norm", Track_ckfInOutTracksFromConversions_trackChi2Norm, &b_Track_ckfInOutTracksFromConversions_trackChi2Norm);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackNdof", Track_ckfInOutTracksFromConversions_trackNdof, &b_Track_ckfInOutTracksFromConversions_trackNdof);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackD0", Track_ckfInOutTracksFromConversions_trackD0, &b_Track_ckfInOutTracksFromConversions_trackD0);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackD0Error", Track_ckfInOutTracksFromConversions_trackD0Error, &b_Track_ckfInOutTracksFromConversions_trackD0Error);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackDsz", Track_ckfInOutTracksFromConversions_trackDsz, &b_Track_ckfInOutTracksFromConversions_trackDsz);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackDszError", Track_ckfInOutTracksFromConversions_trackDszError, &b_Track_ckfInOutTracksFromConversions_trackDszError);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackDxy", Track_ckfInOutTracksFromConversions_trackDxy, &b_Track_ckfInOutTracksFromConversions_trackDxy);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackDxyError", Track_ckfInOutTracksFromConversions_trackDxyError, &b_Track_ckfInOutTracksFromConversions_trackDxyError);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackDz", Track_ckfInOutTracksFromConversions_trackDz, &b_Track_ckfInOutTracksFromConversions_trackDz);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackDzError", Track_ckfInOutTracksFromConversions_trackDzError, &b_Track_ckfInOutTracksFromConversions_trackDzError);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackNumberOfValidHits", Track_ckfInOutTracksFromConversions_trackNumberOfValidHits, &b_Track_ckfInOutTracksFromConversions_trackNumberOfValidHits);
   fChain->SetBranchAddress("Track_ckfInOutTracksFromConversions_trackNumberOfLostHits", Track_ckfInOutTracksFromConversions_trackNumberOfLostHits, &b_Track_ckfInOutTracksFromConversions_trackNumberOfLostHits);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_nRecoTrack", &Track_ckfOutInTracksFromConversions_nRecoTrack, &b_Track_ckfOutInTracksFromConversions_nRecoTrack);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackCharge", Track_ckfOutInTracksFromConversions_trackCharge, &b_Track_ckfOutInTracksFromConversions_trackCharge);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackPx", Track_ckfOutInTracksFromConversions_trackPx, &b_Track_ckfOutInTracksFromConversions_trackPx);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackPy", Track_ckfOutInTracksFromConversions_trackPy, &b_Track_ckfOutInTracksFromConversions_trackPy);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackPz", Track_ckfOutInTracksFromConversions_trackPz, &b_Track_ckfOutInTracksFromConversions_trackPz);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackPt", Track_ckfOutInTracksFromConversions_trackPt, &b_Track_ckfOutInTracksFromConversions_trackPt);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackPtError", Track_ckfOutInTracksFromConversions_trackPtError, &b_Track_ckfOutInTracksFromConversions_trackPtError);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackP", Track_ckfOutInTracksFromConversions_trackP, &b_Track_ckfOutInTracksFromConversions_trackP);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackPhi", Track_ckfOutInTracksFromConversions_trackPhi, &b_Track_ckfOutInTracksFromConversions_trackPhi);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackPhiError", Track_ckfOutInTracksFromConversions_trackPhiError, &b_Track_ckfOutInTracksFromConversions_trackPhiError);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackEta", Track_ckfOutInTracksFromConversions_trackEta, &b_Track_ckfOutInTracksFromConversions_trackEta);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackEtaError", Track_ckfOutInTracksFromConversions_trackEtaError, &b_Track_ckfOutInTracksFromConversions_trackEtaError);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackTheta", Track_ckfOutInTracksFromConversions_trackTheta, &b_Track_ckfOutInTracksFromConversions_trackTheta);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackThetaError", Track_ckfOutInTracksFromConversions_trackThetaError, &b_Track_ckfOutInTracksFromConversions_trackThetaError);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackVx", Track_ckfOutInTracksFromConversions_trackVx, &b_Track_ckfOutInTracksFromConversions_trackVx);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackVy", Track_ckfOutInTracksFromConversions_trackVy, &b_Track_ckfOutInTracksFromConversions_trackVy);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackVz", Track_ckfOutInTracksFromConversions_trackVz, &b_Track_ckfOutInTracksFromConversions_trackVz);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackChi2", Track_ckfOutInTracksFromConversions_trackChi2, &b_Track_ckfOutInTracksFromConversions_trackChi2);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackChi2Norm", Track_ckfOutInTracksFromConversions_trackChi2Norm, &b_Track_ckfOutInTracksFromConversions_trackChi2Norm);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackNdof", Track_ckfOutInTracksFromConversions_trackNdof, &b_Track_ckfOutInTracksFromConversions_trackNdof);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackD0", Track_ckfOutInTracksFromConversions_trackD0, &b_Track_ckfOutInTracksFromConversions_trackD0);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackD0Error", Track_ckfOutInTracksFromConversions_trackD0Error, &b_Track_ckfOutInTracksFromConversions_trackD0Error);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackDsz", Track_ckfOutInTracksFromConversions_trackDsz, &b_Track_ckfOutInTracksFromConversions_trackDsz);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackDszError", Track_ckfOutInTracksFromConversions_trackDszError, &b_Track_ckfOutInTracksFromConversions_trackDszError);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackDxy", Track_ckfOutInTracksFromConversions_trackDxy, &b_Track_ckfOutInTracksFromConversions_trackDxy);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackDxyError", Track_ckfOutInTracksFromConversions_trackDxyError, &b_Track_ckfOutInTracksFromConversions_trackDxyError);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackDz", Track_ckfOutInTracksFromConversions_trackDz, &b_Track_ckfOutInTracksFromConversions_trackDz);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackDzError", Track_ckfOutInTracksFromConversions_trackDzError, &b_Track_ckfOutInTracksFromConversions_trackDzError);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackNumberOfValidHits", Track_ckfOutInTracksFromConversions_trackNumberOfValidHits, &b_Track_ckfOutInTracksFromConversions_trackNumberOfValidHits);
   fChain->SetBranchAddress("Track_ckfOutInTracksFromConversions_trackNumberOfLostHits", Track_ckfOutInTracksFromConversions_trackNumberOfLostHits, &b_Track_ckfOutInTracksFromConversions_trackNumberOfLostHits);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_nGenJet", &JetGen_iterativeCone5GenJets_nGenJet, &b_JetGen_iterativeCone5GenJets_nGenJet);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetCharge", JetGen_iterativeCone5GenJets_genJetCharge, &b_JetGen_iterativeCone5GenJets_genJetCharge);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetPx", JetGen_iterativeCone5GenJets_genJetPx, &b_JetGen_iterativeCone5GenJets_genJetPx);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetPy", JetGen_iterativeCone5GenJets_genJetPy, &b_JetGen_iterativeCone5GenJets_genJetPy);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetPz", JetGen_iterativeCone5GenJets_genJetPz, &b_JetGen_iterativeCone5GenJets_genJetPz);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetPt", JetGen_iterativeCone5GenJets_genJetPt, &b_JetGen_iterativeCone5GenJets_genJetPt);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetP", JetGen_iterativeCone5GenJets_genJetP, &b_JetGen_iterativeCone5GenJets_genJetP);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetPhi", JetGen_iterativeCone5GenJets_genJetPhi, &b_JetGen_iterativeCone5GenJets_genJetPhi);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetEta", JetGen_iterativeCone5GenJets_genJetEta, &b_JetGen_iterativeCone5GenJets_genJetEta);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetTheta", JetGen_iterativeCone5GenJets_genJetTheta, &b_JetGen_iterativeCone5GenJets_genJetTheta);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetEmEnergy", JetGen_iterativeCone5GenJets_genJetEmEnergy, &b_JetGen_iterativeCone5GenJets_genJetEmEnergy);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetHadEnergy", JetGen_iterativeCone5GenJets_genJetHadEnergy, &b_JetGen_iterativeCone5GenJets_genJetHadEnergy);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetInvisibleEnergy", JetGen_iterativeCone5GenJets_genJetInvisibleEnergy, &b_JetGen_iterativeCone5GenJets_genJetInvisibleEnergy);
   fChain->SetBranchAddress("JetGen_iterativeCone5GenJets_genJetAuxiliaryEnergy", JetGen_iterativeCone5GenJets_genJetAuxiliaryEnergy, &b_JetGen_iterativeCone5GenJets_genJetAuxiliaryEnergy);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_nGenJet", &JetGen_sisCone5GenJets_nGenJet, &b_JetGen_sisCone5GenJets_nGenJet);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetCharge", JetGen_sisCone5GenJets_genJetCharge, &b_JetGen_sisCone5GenJets_genJetCharge);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetPx", JetGen_sisCone5GenJets_genJetPx, &b_JetGen_sisCone5GenJets_genJetPx);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetPy", JetGen_sisCone5GenJets_genJetPy, &b_JetGen_sisCone5GenJets_genJetPy);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetPz", JetGen_sisCone5GenJets_genJetPz, &b_JetGen_sisCone5GenJets_genJetPz);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetPt", JetGen_sisCone5GenJets_genJetPt, &b_JetGen_sisCone5GenJets_genJetPt);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetP", JetGen_sisCone5GenJets_genJetP, &b_JetGen_sisCone5GenJets_genJetP);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetPhi", JetGen_sisCone5GenJets_genJetPhi, &b_JetGen_sisCone5GenJets_genJetPhi);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetEta", JetGen_sisCone5GenJets_genJetEta, &b_JetGen_sisCone5GenJets_genJetEta);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetTheta", JetGen_sisCone5GenJets_genJetTheta, &b_JetGen_sisCone5GenJets_genJetTheta);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetEmEnergy", JetGen_sisCone5GenJets_genJetEmEnergy, &b_JetGen_sisCone5GenJets_genJetEmEnergy);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetHadEnergy", JetGen_sisCone5GenJets_genJetHadEnergy, &b_JetGen_sisCone5GenJets_genJetHadEnergy);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetInvisibleEnergy", JetGen_sisCone5GenJets_genJetInvisibleEnergy, &b_JetGen_sisCone5GenJets_genJetInvisibleEnergy);
   fChain->SetBranchAddress("JetGen_sisCone5GenJets_genJetAuxiliaryEnergy", JetGen_sisCone5GenJets_genJetAuxiliaryEnergy, &b_JetGen_sisCone5GenJets_genJetAuxiliaryEnergy);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_nPatJet", &PatJets_allLayer1JetsIC5JPT_nPatJet, &b_PatJets_allLayer1JetsIC5JPT_nPatJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_charge", PatJets_allLayer1JetsIC5JPT_charge, &b_PatJets_allLayer1JetsIC5JPT_charge);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_isCaloJet", PatJets_allLayer1JetsIC5JPT_isCaloJet, &b_PatJets_allLayer1JetsIC5JPT_isCaloJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_isPFJet", PatJets_allLayer1JetsIC5JPT_isPFJet, &b_PatJets_allLayer1JetsIC5JPT_isPFJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_isBasicJet", PatJets_allLayer1JetsIC5JPT_isBasicJet, &b_PatJets_allLayer1JetsIC5JPT_isBasicJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_px", PatJets_allLayer1JetsIC5JPT_px, &b_PatJets_allLayer1JetsIC5JPT_px);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_py", PatJets_allLayer1JetsIC5JPT_py, &b_PatJets_allLayer1JetsIC5JPT_py);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_pz", PatJets_allLayer1JetsIC5JPT_pz, &b_PatJets_allLayer1JetsIC5JPT_pz);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_pt", PatJets_allLayer1JetsIC5JPT_pt, &b_PatJets_allLayer1JetsIC5JPT_pt);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_p", PatJets_allLayer1JetsIC5JPT_p, &b_PatJets_allLayer1JetsIC5JPT_p);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_phi", PatJets_allLayer1JetsIC5JPT_phi, &b_PatJets_allLayer1JetsIC5JPT_phi);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_eta", PatJets_allLayer1JetsIC5JPT_eta, &b_PatJets_allLayer1JetsIC5JPT_eta);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_theta", PatJets_allLayer1JetsIC5JPT_theta, &b_PatJets_allLayer1JetsIC5JPT_theta);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_maxEInEmTowers", PatJets_allLayer1JetsIC5JPT_maxEInEmTowers, &b_PatJets_allLayer1JetsIC5JPT_maxEInEmTowers);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_maxEInHadTowers", PatJets_allLayer1JetsIC5JPT_maxEInHadTowers, &b_PatJets_allLayer1JetsIC5JPT_maxEInHadTowers);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_hadEnergyInHO", PatJets_allLayer1JetsIC5JPT_hadEnergyInHO, &b_PatJets_allLayer1JetsIC5JPT_hadEnergyInHO);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_hadEnergyInHB", PatJets_allLayer1JetsIC5JPT_hadEnergyInHB, &b_PatJets_allLayer1JetsIC5JPT_hadEnergyInHB);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_hadEnergyInHF", PatJets_allLayer1JetsIC5JPT_hadEnergyInHF, &b_PatJets_allLayer1JetsIC5JPT_hadEnergyInHF);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_hadEnergyInHE", PatJets_allLayer1JetsIC5JPT_hadEnergyInHE, &b_PatJets_allLayer1JetsIC5JPT_hadEnergyInHE);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_emEnergyInEB", PatJets_allLayer1JetsIC5JPT_emEnergyInEB, &b_PatJets_allLayer1JetsIC5JPT_emEnergyInEB);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_emEnergyInEE", PatJets_allLayer1JetsIC5JPT_emEnergyInEE, &b_PatJets_allLayer1JetsIC5JPT_emEnergyInEE);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_emEnergyInHF", PatJets_allLayer1JetsIC5JPT_emEnergyInHF, &b_PatJets_allLayer1JetsIC5JPT_emEnergyInHF);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_energyFractionHadronic", PatJets_allLayer1JetsIC5JPT_energyFractionHadronic, &b_PatJets_allLayer1JetsIC5JPT_energyFractionHadronic);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_emEnergyFraction", PatJets_allLayer1JetsIC5JPT_emEnergyFraction, &b_PatJets_allLayer1JetsIC5JPT_emEnergyFraction);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_towersArea", PatJets_allLayer1JetsIC5JPT_towersArea, &b_PatJets_allLayer1JetsIC5JPT_towersArea);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_n90", PatJets_allLayer1JetsIC5JPT_n90, &b_PatJets_allLayer1JetsIC5JPT_n90);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_n60", PatJets_allLayer1JetsIC5JPT_n60, &b_PatJets_allLayer1JetsIC5JPT_n60);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_partonFlavour", PatJets_allLayer1JetsIC5JPT_partonFlavour, &b_PatJets_allLayer1JetsIC5JPT_partonFlavour);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexBJetTags", PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexBJetTags, &b_PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexMVABJetTags", PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexMVABJetTags, &b_PatJets_allLayer1JetsIC5JPT_combinedSecondaryVertexMVABJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_jetBProbabilityBJetTags", PatJets_allLayer1JetsIC5JPT_jetBProbabilityBJetTags, &b_PatJets_allLayer1JetsIC5JPT_jetBProbabilityBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_jetProbabilityBJetTags", PatJets_allLayer1JetsIC5JPT_jetProbabilityBJetTags, &b_PatJets_allLayer1JetsIC5JPT_jetProbabilityBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_simpleSecondaryVertexBJetTags", PatJets_allLayer1JetsIC5JPT_simpleSecondaryVertexBJetTags, &b_PatJets_allLayer1JetsIC5JPT_simpleSecondaryVertexBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_softElectronBJetTags", PatJets_allLayer1JetsIC5JPT_softElectronBJetTags, &b_PatJets_allLayer1JetsIC5JPT_softElectronBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_softMuonBJetTags", PatJets_allLayer1JetsIC5JPT_softMuonBJetTags, &b_PatJets_allLayer1JetsIC5JPT_softMuonBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_softMuonByPtBJetTags", PatJets_allLayer1JetsIC5JPT_softMuonByPtBJetTags, &b_PatJets_allLayer1JetsIC5JPT_softMuonByPtBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_softMuonByIP3dBJetTags", PatJets_allLayer1JetsIC5JPT_softMuonByIP3dBJetTags, &b_PatJets_allLayer1JetsIC5JPT_softMuonByIP3dBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_trackCountingHighEffBJetTags", PatJets_allLayer1JetsIC5JPT_trackCountingHighEffBJetTags, &b_PatJets_allLayer1JetsIC5JPT_trackCountingHighEffBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5JPT_trackCountingHighPurBJetTags", PatJets_allLayer1JetsIC5JPT_trackCountingHighPurBJetTags, &b_PatJets_allLayer1JetsIC5JPT_trackCountingHighPurBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_nPatJet", &PatJets_allLayer1JetsIC5_nPatJet, &b_PatJets_allLayer1JetsIC5_nPatJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_charge", PatJets_allLayer1JetsIC5_charge, &b_PatJets_allLayer1JetsIC5_charge);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_isCaloJet", PatJets_allLayer1JetsIC5_isCaloJet, &b_PatJets_allLayer1JetsIC5_isCaloJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_isPFJet", PatJets_allLayer1JetsIC5_isPFJet, &b_PatJets_allLayer1JetsIC5_isPFJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_isBasicJet", PatJets_allLayer1JetsIC5_isBasicJet, &b_PatJets_allLayer1JetsIC5_isBasicJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_px", PatJets_allLayer1JetsIC5_px, &b_PatJets_allLayer1JetsIC5_px);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_py", PatJets_allLayer1JetsIC5_py, &b_PatJets_allLayer1JetsIC5_py);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_pz", PatJets_allLayer1JetsIC5_pz, &b_PatJets_allLayer1JetsIC5_pz);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_pt", PatJets_allLayer1JetsIC5_pt, &b_PatJets_allLayer1JetsIC5_pt);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_p", PatJets_allLayer1JetsIC5_p, &b_PatJets_allLayer1JetsIC5_p);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_phi", PatJets_allLayer1JetsIC5_phi, &b_PatJets_allLayer1JetsIC5_phi);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_eta", PatJets_allLayer1JetsIC5_eta, &b_PatJets_allLayer1JetsIC5_eta);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_theta", PatJets_allLayer1JetsIC5_theta, &b_PatJets_allLayer1JetsIC5_theta);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_maxEInEmTowers", PatJets_allLayer1JetsIC5_maxEInEmTowers, &b_PatJets_allLayer1JetsIC5_maxEInEmTowers);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_maxEInHadTowers", PatJets_allLayer1JetsIC5_maxEInHadTowers, &b_PatJets_allLayer1JetsIC5_maxEInHadTowers);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_hadEnergyInHO", PatJets_allLayer1JetsIC5_hadEnergyInHO, &b_PatJets_allLayer1JetsIC5_hadEnergyInHO);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_hadEnergyInHB", PatJets_allLayer1JetsIC5_hadEnergyInHB, &b_PatJets_allLayer1JetsIC5_hadEnergyInHB);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_hadEnergyInHF", PatJets_allLayer1JetsIC5_hadEnergyInHF, &b_PatJets_allLayer1JetsIC5_hadEnergyInHF);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_hadEnergyInHE", PatJets_allLayer1JetsIC5_hadEnergyInHE, &b_PatJets_allLayer1JetsIC5_hadEnergyInHE);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_emEnergyInEB", PatJets_allLayer1JetsIC5_emEnergyInEB, &b_PatJets_allLayer1JetsIC5_emEnergyInEB);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_emEnergyInEE", PatJets_allLayer1JetsIC5_emEnergyInEE, &b_PatJets_allLayer1JetsIC5_emEnergyInEE);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_emEnergyInHF", PatJets_allLayer1JetsIC5_emEnergyInHF, &b_PatJets_allLayer1JetsIC5_emEnergyInHF);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_energyFractionHadronic", PatJets_allLayer1JetsIC5_energyFractionHadronic, &b_PatJets_allLayer1JetsIC5_energyFractionHadronic);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_emEnergyFraction", PatJets_allLayer1JetsIC5_emEnergyFraction, &b_PatJets_allLayer1JetsIC5_emEnergyFraction);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_towersArea", PatJets_allLayer1JetsIC5_towersArea, &b_PatJets_allLayer1JetsIC5_towersArea);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_n90", PatJets_allLayer1JetsIC5_n90, &b_PatJets_allLayer1JetsIC5_n90);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_n60", PatJets_allLayer1JetsIC5_n60, &b_PatJets_allLayer1JetsIC5_n60);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_partonFlavour", PatJets_allLayer1JetsIC5_partonFlavour, &b_PatJets_allLayer1JetsIC5_partonFlavour);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_combinedSecondaryVertexBJetTags", PatJets_allLayer1JetsIC5_combinedSecondaryVertexBJetTags, &b_PatJets_allLayer1JetsIC5_combinedSecondaryVertexBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_combinedSecondaryVertexMVABJetTags", PatJets_allLayer1JetsIC5_combinedSecondaryVertexMVABJetTags, &b_PatJets_allLayer1JetsIC5_combinedSecondaryVertexMVABJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_jetBProbabilityBJetTags", PatJets_allLayer1JetsIC5_jetBProbabilityBJetTags, &b_PatJets_allLayer1JetsIC5_jetBProbabilityBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_jetProbabilityBJetTags", PatJets_allLayer1JetsIC5_jetProbabilityBJetTags, &b_PatJets_allLayer1JetsIC5_jetProbabilityBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_simpleSecondaryVertexBJetTags", PatJets_allLayer1JetsIC5_simpleSecondaryVertexBJetTags, &b_PatJets_allLayer1JetsIC5_simpleSecondaryVertexBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_softElectronBJetTags", PatJets_allLayer1JetsIC5_softElectronBJetTags, &b_PatJets_allLayer1JetsIC5_softElectronBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_softMuonBJetTags", PatJets_allLayer1JetsIC5_softMuonBJetTags, &b_PatJets_allLayer1JetsIC5_softMuonBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_softMuonByPtBJetTags", PatJets_allLayer1JetsIC5_softMuonByPtBJetTags, &b_PatJets_allLayer1JetsIC5_softMuonByPtBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_softMuonByIP3dBJetTags", PatJets_allLayer1JetsIC5_softMuonByIP3dBJetTags, &b_PatJets_allLayer1JetsIC5_softMuonByIP3dBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_trackCountingHighEffBJetTags", PatJets_allLayer1JetsIC5_trackCountingHighEffBJetTags, &b_PatJets_allLayer1JetsIC5_trackCountingHighEffBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsIC5_trackCountingHighPurBJetTags", PatJets_allLayer1JetsIC5_trackCountingHighPurBJetTags, &b_PatJets_allLayer1JetsIC5_trackCountingHighPurBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_nPatJet", &PatJets_allLayer1JetsSC5_nPatJet, &b_PatJets_allLayer1JetsSC5_nPatJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_charge", PatJets_allLayer1JetsSC5_charge, &b_PatJets_allLayer1JetsSC5_charge);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_isCaloJet", PatJets_allLayer1JetsSC5_isCaloJet, &b_PatJets_allLayer1JetsSC5_isCaloJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_isPFJet", PatJets_allLayer1JetsSC5_isPFJet, &b_PatJets_allLayer1JetsSC5_isPFJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_isBasicJet", PatJets_allLayer1JetsSC5_isBasicJet, &b_PatJets_allLayer1JetsSC5_isBasicJet);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_px", PatJets_allLayer1JetsSC5_px, &b_PatJets_allLayer1JetsSC5_px);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_py", PatJets_allLayer1JetsSC5_py, &b_PatJets_allLayer1JetsSC5_py);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_pz", PatJets_allLayer1JetsSC5_pz, &b_PatJets_allLayer1JetsSC5_pz);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_pt", PatJets_allLayer1JetsSC5_pt, &b_PatJets_allLayer1JetsSC5_pt);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_p", PatJets_allLayer1JetsSC5_p, &b_PatJets_allLayer1JetsSC5_p);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_phi", PatJets_allLayer1JetsSC5_phi, &b_PatJets_allLayer1JetsSC5_phi);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_eta", PatJets_allLayer1JetsSC5_eta, &b_PatJets_allLayer1JetsSC5_eta);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_theta", PatJets_allLayer1JetsSC5_theta, &b_PatJets_allLayer1JetsSC5_theta);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_maxEInEmTowers", PatJets_allLayer1JetsSC5_maxEInEmTowers, &b_PatJets_allLayer1JetsSC5_maxEInEmTowers);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_maxEInHadTowers", PatJets_allLayer1JetsSC5_maxEInHadTowers, &b_PatJets_allLayer1JetsSC5_maxEInHadTowers);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_hadEnergyInHO", PatJets_allLayer1JetsSC5_hadEnergyInHO, &b_PatJets_allLayer1JetsSC5_hadEnergyInHO);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_hadEnergyInHB", PatJets_allLayer1JetsSC5_hadEnergyInHB, &b_PatJets_allLayer1JetsSC5_hadEnergyInHB);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_hadEnergyInHF", PatJets_allLayer1JetsSC5_hadEnergyInHF, &b_PatJets_allLayer1JetsSC5_hadEnergyInHF);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_hadEnergyInHE", PatJets_allLayer1JetsSC5_hadEnergyInHE, &b_PatJets_allLayer1JetsSC5_hadEnergyInHE);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_emEnergyInEB", PatJets_allLayer1JetsSC5_emEnergyInEB, &b_PatJets_allLayer1JetsSC5_emEnergyInEB);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_emEnergyInEE", PatJets_allLayer1JetsSC5_emEnergyInEE, &b_PatJets_allLayer1JetsSC5_emEnergyInEE);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_emEnergyInHF", PatJets_allLayer1JetsSC5_emEnergyInHF, &b_PatJets_allLayer1JetsSC5_emEnergyInHF);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_energyFractionHadronic", PatJets_allLayer1JetsSC5_energyFractionHadronic, &b_PatJets_allLayer1JetsSC5_energyFractionHadronic);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_emEnergyFraction", PatJets_allLayer1JetsSC5_emEnergyFraction, &b_PatJets_allLayer1JetsSC5_emEnergyFraction);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_towersArea", PatJets_allLayer1JetsSC5_towersArea, &b_PatJets_allLayer1JetsSC5_towersArea);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_n90", PatJets_allLayer1JetsSC5_n90, &b_PatJets_allLayer1JetsSC5_n90);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_n60", PatJets_allLayer1JetsSC5_n60, &b_PatJets_allLayer1JetsSC5_n60);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_partonFlavour", PatJets_allLayer1JetsSC5_partonFlavour, &b_PatJets_allLayer1JetsSC5_partonFlavour);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_combinedSecondaryVertexBJetTags", PatJets_allLayer1JetsSC5_combinedSecondaryVertexBJetTags, &b_PatJets_allLayer1JetsSC5_combinedSecondaryVertexBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_combinedSecondaryVertexMVABJetTags", PatJets_allLayer1JetsSC5_combinedSecondaryVertexMVABJetTags, &b_PatJets_allLayer1JetsSC5_combinedSecondaryVertexMVABJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_jetBProbabilityBJetTags", PatJets_allLayer1JetsSC5_jetBProbabilityBJetTags, &b_PatJets_allLayer1JetsSC5_jetBProbabilityBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_jetProbabilityBJetTags", PatJets_allLayer1JetsSC5_jetProbabilityBJetTags, &b_PatJets_allLayer1JetsSC5_jetProbabilityBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_simpleSecondaryVertexBJetTags", PatJets_allLayer1JetsSC5_simpleSecondaryVertexBJetTags, &b_PatJets_allLayer1JetsSC5_simpleSecondaryVertexBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_softElectronBJetTags", PatJets_allLayer1JetsSC5_softElectronBJetTags, &b_PatJets_allLayer1JetsSC5_softElectronBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_softMuonBJetTags", PatJets_allLayer1JetsSC5_softMuonBJetTags, &b_PatJets_allLayer1JetsSC5_softMuonBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_softMuonByPtBJetTags", PatJets_allLayer1JetsSC5_softMuonByPtBJetTags, &b_PatJets_allLayer1JetsSC5_softMuonByPtBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_softMuonByIP3dBJetTags", PatJets_allLayer1JetsSC5_softMuonByIP3dBJetTags, &b_PatJets_allLayer1JetsSC5_softMuonByIP3dBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_trackCountingHighEffBJetTags", PatJets_allLayer1JetsSC5_trackCountingHighEffBJetTags, &b_PatJets_allLayer1JetsSC5_trackCountingHighEffBJetTags);
   fChain->SetBranchAddress("PatJets_allLayer1JetsSC5_trackCountingHighPurBJetTags", PatJets_allLayer1JetsSC5_trackCountingHighPurBJetTags, &b_PatJets_allLayer1JetsSC5_trackCountingHighPurBJetTags);
   fChain->SetBranchAddress("METGen_genMet_genMETPx", &METGen_genMet_genMETPx, &b_METGen_genMet_genMETPx);
   fChain->SetBranchAddress("METGen_genMet_genMETPy", &METGen_genMet_genMETPy, &b_METGen_genMet_genMETPy);
   fChain->SetBranchAddress("METGen_genMet_genMETEt", &METGen_genMet_genMETEt, &b_METGen_genMet_genMETEt);
   fChain->SetBranchAddress("METGen_genMet_genMETPhi", &METGen_genMet_genMETPhi, &b_METGen_genMet_genMETPhi);
   fChain->SetBranchAddress("METGen_genMet_genMETsumEt", &METGen_genMet_genMETsumEt, &b_METGen_genMet_genMETsumEt);
   fChain->SetBranchAddress("METGen_genMet_genMETmEtSig", &METGen_genMet_genMETmEtSig, &b_METGen_genMet_genMETmEtSig);
   fChain->SetBranchAddress("METGen_genMet_genMETE_longitudinal", &METGen_genMet_genMETE_longitudinal, &b_METGen_genMet_genMETE_longitudinal);
   fChain->SetBranchAddress("METGen_genMet_genMETEmEnergy", &METGen_genMet_genMETEmEnergy, &b_METGen_genMet_genMETEmEnergy);
   fChain->SetBranchAddress("METGen_genMet_genMETHadEnergy", &METGen_genMet_genMETHadEnergy, &b_METGen_genMet_genMETHadEnergy);
   fChain->SetBranchAddress("METGen_genMet_genMETInvisibleEnergy", &METGen_genMet_genMETInvisibleEnergy, &b_METGen_genMet_genMETInvisibleEnergy);
   fChain->SetBranchAddress("METGen_genMet_genMETAuxiliaryEnergy", &METGen_genMet_genMETAuxiliaryEnergy, &b_METGen_genMet_genMETAuxiliaryEnergy);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETPx", &METGen_genMetNoNuBSM_genMETPx, &b_METGen_genMetNoNuBSM_genMETPx);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETPy", &METGen_genMetNoNuBSM_genMETPy, &b_METGen_genMetNoNuBSM_genMETPy);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETEt", &METGen_genMetNoNuBSM_genMETEt, &b_METGen_genMetNoNuBSM_genMETEt);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETPhi", &METGen_genMetNoNuBSM_genMETPhi, &b_METGen_genMetNoNuBSM_genMETPhi);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETsumEt", &METGen_genMetNoNuBSM_genMETsumEt, &b_METGen_genMetNoNuBSM_genMETsumEt);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETmEtSig", &METGen_genMetNoNuBSM_genMETmEtSig, &b_METGen_genMetNoNuBSM_genMETmEtSig);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETE_longitudinal", &METGen_genMetNoNuBSM_genMETE_longitudinal, &b_METGen_genMetNoNuBSM_genMETE_longitudinal);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETEmEnergy", &METGen_genMetNoNuBSM_genMETEmEnergy, &b_METGen_genMetNoNuBSM_genMETEmEnergy);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETHadEnergy", &METGen_genMetNoNuBSM_genMETHadEnergy, &b_METGen_genMetNoNuBSM_genMETHadEnergy);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETInvisibleEnergy", &METGen_genMetNoNuBSM_genMETInvisibleEnergy, &b_METGen_genMetNoNuBSM_genMETInvisibleEnergy);
   fChain->SetBranchAddress("METGen_genMetNoNuBSM_genMETAuxiliaryEnergy", &METGen_genMetNoNuBSM_genMETAuxiliaryEnergy, &b_METGen_genMetNoNuBSM_genMETAuxiliaryEnergy);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_px", &PatMET_allLayer1METsIC5_px, &b_PatMET_allLayer1METsIC5_px);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_py", &PatMET_allLayer1METsIC5_py, &b_PatMET_allLayer1METsIC5_py);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_et", &PatMET_allLayer1METsIC5_et, &b_PatMET_allLayer1METsIC5_et);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_phi", &PatMET_allLayer1METsIC5_phi, &b_PatMET_allLayer1METsIC5_phi);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_sumEt", &PatMET_allLayer1METsIC5_sumEt, &b_PatMET_allLayer1METsIC5_sumEt);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_mEtSig", &PatMET_allLayer1METsIC5_mEtSig, &b_PatMET_allLayer1METsIC5_mEtSig);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_e_longitudinal", &PatMET_allLayer1METsIC5_e_longitudinal, &b_PatMET_allLayer1METsIC5_e_longitudinal);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_nCorrections", &PatMET_allLayer1METsIC5_nCorrections, &b_PatMET_allLayer1METsIC5_nCorrections);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_corExUncorrALL", &PatMET_allLayer1METsIC5_corExUncorrALL, &b_PatMET_allLayer1METsIC5_corExUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_corEyUncorrALL", &PatMET_allLayer1METsIC5_corEyUncorrALL, &b_PatMET_allLayer1METsIC5_corEyUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_corSumEtUncorrALL", &PatMET_allLayer1METsIC5_corSumEtUncorrALL, &b_PatMET_allLayer1METsIC5_corSumEtUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_uncorrectedPtUncorrALL", &PatMET_allLayer1METsIC5_uncorrectedPtUncorrALL, &b_PatMET_allLayer1METsIC5_uncorrectedPtUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_uncorrectedPhiUncorrALL", &PatMET_allLayer1METsIC5_uncorrectedPhiUncorrALL, &b_PatMET_allLayer1METsIC5_uncorrectedPhiUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_corExUncorrJES", &PatMET_allLayer1METsIC5_corExUncorrJES, &b_PatMET_allLayer1METsIC5_corExUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_corEyUncorrJES", &PatMET_allLayer1METsIC5_corEyUncorrJES, &b_PatMET_allLayer1METsIC5_corEyUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_corSumEtUncorrJES", &PatMET_allLayer1METsIC5_corSumEtUncorrJES, &b_PatMET_allLayer1METsIC5_corSumEtUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_uncorrectedPtUncorrJES", &PatMET_allLayer1METsIC5_uncorrectedPtUncorrJES, &b_PatMET_allLayer1METsIC5_uncorrectedPtUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_uncorrectedPhiUncorrJES", &PatMET_allLayer1METsIC5_uncorrectedPhiUncorrJES, &b_PatMET_allLayer1METsIC5_uncorrectedPhiUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_corExUncorrMUON", &PatMET_allLayer1METsIC5_corExUncorrMUON, &b_PatMET_allLayer1METsIC5_corExUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_corEyUncorrMUON", &PatMET_allLayer1METsIC5_corEyUncorrMUON, &b_PatMET_allLayer1METsIC5_corEyUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_corSumEtUncorrMUON", &PatMET_allLayer1METsIC5_corSumEtUncorrMUON, &b_PatMET_allLayer1METsIC5_corSumEtUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_uncorrectedPtUncorrMUON", &PatMET_allLayer1METsIC5_uncorrectedPtUncorrMUON, &b_PatMET_allLayer1METsIC5_uncorrectedPtUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_uncorrectedPhiUncorrMUON", &PatMET_allLayer1METsIC5_uncorrectedPhiUncorrMUON, &b_PatMET_allLayer1METsIC5_uncorrectedPhiUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_isCaloMET", &PatMET_allLayer1METsIC5_isCaloMET, &b_PatMET_allLayer1METsIC5_isCaloMET);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_isRecoMET", &PatMET_allLayer1METsIC5_isRecoMET, &b_PatMET_allLayer1METsIC5_isRecoMET);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_maxEtInEmTowers", &PatMET_allLayer1METsIC5_maxEtInEmTowers, &b_PatMET_allLayer1METsIC5_maxEtInEmTowers);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_maxEtInHadTowers", &PatMET_allLayer1METsIC5_maxEtInHadTowers, &b_PatMET_allLayer1METsIC5_maxEtInHadTowers);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_etFractionHadronic", &PatMET_allLayer1METsIC5_etFractionHadronic, &b_PatMET_allLayer1METsIC5_etFractionHadronic);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_emEtFraction", &PatMET_allLayer1METsIC5_emEtFraction, &b_PatMET_allLayer1METsIC5_emEtFraction);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_hadEtInHB", &PatMET_allLayer1METsIC5_hadEtInHB, &b_PatMET_allLayer1METsIC5_hadEtInHB);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_hadEtInHO", &PatMET_allLayer1METsIC5_hadEtInHO, &b_PatMET_allLayer1METsIC5_hadEtInHO);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_hadEtInHE", &PatMET_allLayer1METsIC5_hadEtInHE, &b_PatMET_allLayer1METsIC5_hadEtInHE);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_hadEtInHF", &PatMET_allLayer1METsIC5_hadEtInHF, &b_PatMET_allLayer1METsIC5_hadEtInHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_emEtInEB", &PatMET_allLayer1METsIC5_emEtInEB, &b_PatMET_allLayer1METsIC5_emEtInEB);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_emEtInEE", &PatMET_allLayer1METsIC5_emEtInEE, &b_PatMET_allLayer1METsIC5_emEtInEE);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_emEtInHF", &PatMET_allLayer1METsIC5_emEtInHF, &b_PatMET_allLayer1METsIC5_emEtInHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_metSignificance", &PatMET_allLayer1METsIC5_metSignificance, &b_PatMET_allLayer1METsIC5_metSignificance);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_CaloSETInpHF", &PatMET_allLayer1METsIC5_CaloSETInpHF, &b_PatMET_allLayer1METsIC5_CaloSETInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_CaloSETInmHF", &PatMET_allLayer1METsIC5_CaloSETInmHF, &b_PatMET_allLayer1METsIC5_CaloSETInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_CaloMETInpHF", &PatMET_allLayer1METsIC5_CaloMETInpHF, &b_PatMET_allLayer1METsIC5_CaloMETInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_CaloMETInmHF", &PatMET_allLayer1METsIC5_CaloMETInmHF, &b_PatMET_allLayer1METsIC5_CaloMETInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_CaloMETPhiInpHF", &PatMET_allLayer1METsIC5_CaloMETPhiInpHF, &b_PatMET_allLayer1METsIC5_CaloMETPhiInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsIC5_CaloMETPhiInmHF", &PatMET_allLayer1METsIC5_CaloMETPhiInmHF, &b_PatMET_allLayer1METsIC5_CaloMETPhiInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_px", &PatMET_allLayer1METsPF_px, &b_PatMET_allLayer1METsPF_px);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_py", &PatMET_allLayer1METsPF_py, &b_PatMET_allLayer1METsPF_py);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_et", &PatMET_allLayer1METsPF_et, &b_PatMET_allLayer1METsPF_et);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_phi", &PatMET_allLayer1METsPF_phi, &b_PatMET_allLayer1METsPF_phi);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_sumEt", &PatMET_allLayer1METsPF_sumEt, &b_PatMET_allLayer1METsPF_sumEt);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_mEtSig", &PatMET_allLayer1METsPF_mEtSig, &b_PatMET_allLayer1METsPF_mEtSig);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_e_longitudinal", &PatMET_allLayer1METsPF_e_longitudinal, &b_PatMET_allLayer1METsPF_e_longitudinal);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_nCorrections", &PatMET_allLayer1METsPF_nCorrections, &b_PatMET_allLayer1METsPF_nCorrections);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_corExUncorrALL", &PatMET_allLayer1METsPF_corExUncorrALL, &b_PatMET_allLayer1METsPF_corExUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_corEyUncorrALL", &PatMET_allLayer1METsPF_corEyUncorrALL, &b_PatMET_allLayer1METsPF_corEyUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_corSumEtUncorrALL", &PatMET_allLayer1METsPF_corSumEtUncorrALL, &b_PatMET_allLayer1METsPF_corSumEtUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_uncorrectedPtUncorrALL", &PatMET_allLayer1METsPF_uncorrectedPtUncorrALL, &b_PatMET_allLayer1METsPF_uncorrectedPtUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_uncorrectedPhiUncorrALL", &PatMET_allLayer1METsPF_uncorrectedPhiUncorrALL, &b_PatMET_allLayer1METsPF_uncorrectedPhiUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_corExUncorrJES", &PatMET_allLayer1METsPF_corExUncorrJES, &b_PatMET_allLayer1METsPF_corExUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_corEyUncorrJES", &PatMET_allLayer1METsPF_corEyUncorrJES, &b_PatMET_allLayer1METsPF_corEyUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_corSumEtUncorrJES", &PatMET_allLayer1METsPF_corSumEtUncorrJES, &b_PatMET_allLayer1METsPF_corSumEtUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_uncorrectedPtUncorrJES", &PatMET_allLayer1METsPF_uncorrectedPtUncorrJES, &b_PatMET_allLayer1METsPF_uncorrectedPtUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_uncorrectedPhiUncorrJES", &PatMET_allLayer1METsPF_uncorrectedPhiUncorrJES, &b_PatMET_allLayer1METsPF_uncorrectedPhiUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_corExUncorrMUON", &PatMET_allLayer1METsPF_corExUncorrMUON, &b_PatMET_allLayer1METsPF_corExUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_corEyUncorrMUON", &PatMET_allLayer1METsPF_corEyUncorrMUON, &b_PatMET_allLayer1METsPF_corEyUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_corSumEtUncorrMUON", &PatMET_allLayer1METsPF_corSumEtUncorrMUON, &b_PatMET_allLayer1METsPF_corSumEtUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_uncorrectedPtUncorrMUON", &PatMET_allLayer1METsPF_uncorrectedPtUncorrMUON, &b_PatMET_allLayer1METsPF_uncorrectedPtUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_uncorrectedPhiUncorrMUON", &PatMET_allLayer1METsPF_uncorrectedPhiUncorrMUON, &b_PatMET_allLayer1METsPF_uncorrectedPhiUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_isCaloMET", &PatMET_allLayer1METsPF_isCaloMET, &b_PatMET_allLayer1METsPF_isCaloMET);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_isRecoMET", &PatMET_allLayer1METsPF_isRecoMET, &b_PatMET_allLayer1METsPF_isRecoMET);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_maxEtInEmTowers", &PatMET_allLayer1METsPF_maxEtInEmTowers, &b_PatMET_allLayer1METsPF_maxEtInEmTowers);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_maxEtInHadTowers", &PatMET_allLayer1METsPF_maxEtInHadTowers, &b_PatMET_allLayer1METsPF_maxEtInHadTowers);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_etFractionHadronic", &PatMET_allLayer1METsPF_etFractionHadronic, &b_PatMET_allLayer1METsPF_etFractionHadronic);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_emEtFraction", &PatMET_allLayer1METsPF_emEtFraction, &b_PatMET_allLayer1METsPF_emEtFraction);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_hadEtInHB", &PatMET_allLayer1METsPF_hadEtInHB, &b_PatMET_allLayer1METsPF_hadEtInHB);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_hadEtInHO", &PatMET_allLayer1METsPF_hadEtInHO, &b_PatMET_allLayer1METsPF_hadEtInHO);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_hadEtInHE", &PatMET_allLayer1METsPF_hadEtInHE, &b_PatMET_allLayer1METsPF_hadEtInHE);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_hadEtInHF", &PatMET_allLayer1METsPF_hadEtInHF, &b_PatMET_allLayer1METsPF_hadEtInHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_emEtInEB", &PatMET_allLayer1METsPF_emEtInEB, &b_PatMET_allLayer1METsPF_emEtInEB);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_emEtInEE", &PatMET_allLayer1METsPF_emEtInEE, &b_PatMET_allLayer1METsPF_emEtInEE);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_emEtInHF", &PatMET_allLayer1METsPF_emEtInHF, &b_PatMET_allLayer1METsPF_emEtInHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_metSignificance", &PatMET_allLayer1METsPF_metSignificance, &b_PatMET_allLayer1METsPF_metSignificance);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_CaloSETInpHF", &PatMET_allLayer1METsPF_CaloSETInpHF, &b_PatMET_allLayer1METsPF_CaloSETInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_CaloSETInmHF", &PatMET_allLayer1METsPF_CaloSETInmHF, &b_PatMET_allLayer1METsPF_CaloSETInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_CaloMETInpHF", &PatMET_allLayer1METsPF_CaloMETInpHF, &b_PatMET_allLayer1METsPF_CaloMETInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_CaloMETInmHF", &PatMET_allLayer1METsPF_CaloMETInmHF, &b_PatMET_allLayer1METsPF_CaloMETInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_CaloMETPhiInpHF", &PatMET_allLayer1METsPF_CaloMETPhiInpHF, &b_PatMET_allLayer1METsPF_CaloMETPhiInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsPF_CaloMETPhiInmHF", &PatMET_allLayer1METsPF_CaloMETPhiInmHF, &b_PatMET_allLayer1METsPF_CaloMETPhiInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_px", &PatMET_allLayer1METsSC5_px, &b_PatMET_allLayer1METsSC5_px);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_py", &PatMET_allLayer1METsSC5_py, &b_PatMET_allLayer1METsSC5_py);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_et", &PatMET_allLayer1METsSC5_et, &b_PatMET_allLayer1METsSC5_et);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_phi", &PatMET_allLayer1METsSC5_phi, &b_PatMET_allLayer1METsSC5_phi);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_sumEt", &PatMET_allLayer1METsSC5_sumEt, &b_PatMET_allLayer1METsSC5_sumEt);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_mEtSig", &PatMET_allLayer1METsSC5_mEtSig, &b_PatMET_allLayer1METsSC5_mEtSig);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_e_longitudinal", &PatMET_allLayer1METsSC5_e_longitudinal, &b_PatMET_allLayer1METsSC5_e_longitudinal);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_nCorrections", &PatMET_allLayer1METsSC5_nCorrections, &b_PatMET_allLayer1METsSC5_nCorrections);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_corExUncorrALL", &PatMET_allLayer1METsSC5_corExUncorrALL, &b_PatMET_allLayer1METsSC5_corExUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_corEyUncorrALL", &PatMET_allLayer1METsSC5_corEyUncorrALL, &b_PatMET_allLayer1METsSC5_corEyUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_corSumEtUncorrALL", &PatMET_allLayer1METsSC5_corSumEtUncorrALL, &b_PatMET_allLayer1METsSC5_corSumEtUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_uncorrectedPtUncorrALL", &PatMET_allLayer1METsSC5_uncorrectedPtUncorrALL, &b_PatMET_allLayer1METsSC5_uncorrectedPtUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_uncorrectedPhiUncorrALL", &PatMET_allLayer1METsSC5_uncorrectedPhiUncorrALL, &b_PatMET_allLayer1METsSC5_uncorrectedPhiUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_corExUncorrJES", &PatMET_allLayer1METsSC5_corExUncorrJES, &b_PatMET_allLayer1METsSC5_corExUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_corEyUncorrJES", &PatMET_allLayer1METsSC5_corEyUncorrJES, &b_PatMET_allLayer1METsSC5_corEyUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_corSumEtUncorrJES", &PatMET_allLayer1METsSC5_corSumEtUncorrJES, &b_PatMET_allLayer1METsSC5_corSumEtUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_uncorrectedPtUncorrJES", &PatMET_allLayer1METsSC5_uncorrectedPtUncorrJES, &b_PatMET_allLayer1METsSC5_uncorrectedPtUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_uncorrectedPhiUncorrJES", &PatMET_allLayer1METsSC5_uncorrectedPhiUncorrJES, &b_PatMET_allLayer1METsSC5_uncorrectedPhiUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_corExUncorrMUON", &PatMET_allLayer1METsSC5_corExUncorrMUON, &b_PatMET_allLayer1METsSC5_corExUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_corEyUncorrMUON", &PatMET_allLayer1METsSC5_corEyUncorrMUON, &b_PatMET_allLayer1METsSC5_corEyUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_corSumEtUncorrMUON", &PatMET_allLayer1METsSC5_corSumEtUncorrMUON, &b_PatMET_allLayer1METsSC5_corSumEtUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_uncorrectedPtUncorrMUON", &PatMET_allLayer1METsSC5_uncorrectedPtUncorrMUON, &b_PatMET_allLayer1METsSC5_uncorrectedPtUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_uncorrectedPhiUncorrMUON", &PatMET_allLayer1METsSC5_uncorrectedPhiUncorrMUON, &b_PatMET_allLayer1METsSC5_uncorrectedPhiUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_isCaloMET", &PatMET_allLayer1METsSC5_isCaloMET, &b_PatMET_allLayer1METsSC5_isCaloMET);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_isRecoMET", &PatMET_allLayer1METsSC5_isRecoMET, &b_PatMET_allLayer1METsSC5_isRecoMET);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_maxEtInEmTowers", &PatMET_allLayer1METsSC5_maxEtInEmTowers, &b_PatMET_allLayer1METsSC5_maxEtInEmTowers);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_maxEtInHadTowers", &PatMET_allLayer1METsSC5_maxEtInHadTowers, &b_PatMET_allLayer1METsSC5_maxEtInHadTowers);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_etFractionHadronic", &PatMET_allLayer1METsSC5_etFractionHadronic, &b_PatMET_allLayer1METsSC5_etFractionHadronic);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_emEtFraction", &PatMET_allLayer1METsSC5_emEtFraction, &b_PatMET_allLayer1METsSC5_emEtFraction);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_hadEtInHB", &PatMET_allLayer1METsSC5_hadEtInHB, &b_PatMET_allLayer1METsSC5_hadEtInHB);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_hadEtInHO", &PatMET_allLayer1METsSC5_hadEtInHO, &b_PatMET_allLayer1METsSC5_hadEtInHO);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_hadEtInHE", &PatMET_allLayer1METsSC5_hadEtInHE, &b_PatMET_allLayer1METsSC5_hadEtInHE);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_hadEtInHF", &PatMET_allLayer1METsSC5_hadEtInHF, &b_PatMET_allLayer1METsSC5_hadEtInHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_emEtInEB", &PatMET_allLayer1METsSC5_emEtInEB, &b_PatMET_allLayer1METsSC5_emEtInEB);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_emEtInEE", &PatMET_allLayer1METsSC5_emEtInEE, &b_PatMET_allLayer1METsSC5_emEtInEE);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_emEtInHF", &PatMET_allLayer1METsSC5_emEtInHF, &b_PatMET_allLayer1METsSC5_emEtInHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_metSignificance", &PatMET_allLayer1METsSC5_metSignificance, &b_PatMET_allLayer1METsSC5_metSignificance);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_CaloSETInpHF", &PatMET_allLayer1METsSC5_CaloSETInpHF, &b_PatMET_allLayer1METsSC5_CaloSETInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_CaloSETInmHF", &PatMET_allLayer1METsSC5_CaloSETInmHF, &b_PatMET_allLayer1METsSC5_CaloSETInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_CaloMETInpHF", &PatMET_allLayer1METsSC5_CaloMETInpHF, &b_PatMET_allLayer1METsSC5_CaloMETInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_CaloMETInmHF", &PatMET_allLayer1METsSC5_CaloMETInmHF, &b_PatMET_allLayer1METsSC5_CaloMETInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_CaloMETPhiInpHF", &PatMET_allLayer1METsSC5_CaloMETPhiInpHF, &b_PatMET_allLayer1METsSC5_CaloMETPhiInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METsSC5_CaloMETPhiInmHF", &PatMET_allLayer1METsSC5_CaloMETPhiInmHF, &b_PatMET_allLayer1METsSC5_CaloMETPhiInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_px", &PatMET_allLayer1METstcMET_px, &b_PatMET_allLayer1METstcMET_px);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_py", &PatMET_allLayer1METstcMET_py, &b_PatMET_allLayer1METstcMET_py);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_et", &PatMET_allLayer1METstcMET_et, &b_PatMET_allLayer1METstcMET_et);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_phi", &PatMET_allLayer1METstcMET_phi, &b_PatMET_allLayer1METstcMET_phi);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_sumEt", &PatMET_allLayer1METstcMET_sumEt, &b_PatMET_allLayer1METstcMET_sumEt);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_mEtSig", &PatMET_allLayer1METstcMET_mEtSig, &b_PatMET_allLayer1METstcMET_mEtSig);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_e_longitudinal", &PatMET_allLayer1METstcMET_e_longitudinal, &b_PatMET_allLayer1METstcMET_e_longitudinal);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_nCorrections", &PatMET_allLayer1METstcMET_nCorrections, &b_PatMET_allLayer1METstcMET_nCorrections);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_corExUncorrALL", &PatMET_allLayer1METstcMET_corExUncorrALL, &b_PatMET_allLayer1METstcMET_corExUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_corEyUncorrALL", &PatMET_allLayer1METstcMET_corEyUncorrALL, &b_PatMET_allLayer1METstcMET_corEyUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_corSumEtUncorrALL", &PatMET_allLayer1METstcMET_corSumEtUncorrALL, &b_PatMET_allLayer1METstcMET_corSumEtUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_uncorrectedPtUncorrALL", &PatMET_allLayer1METstcMET_uncorrectedPtUncorrALL, &b_PatMET_allLayer1METstcMET_uncorrectedPtUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_uncorrectedPhiUncorrALL", &PatMET_allLayer1METstcMET_uncorrectedPhiUncorrALL, &b_PatMET_allLayer1METstcMET_uncorrectedPhiUncorrALL);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_corExUncorrJES", &PatMET_allLayer1METstcMET_corExUncorrJES, &b_PatMET_allLayer1METstcMET_corExUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_corEyUncorrJES", &PatMET_allLayer1METstcMET_corEyUncorrJES, &b_PatMET_allLayer1METstcMET_corEyUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_corSumEtUncorrJES", &PatMET_allLayer1METstcMET_corSumEtUncorrJES, &b_PatMET_allLayer1METstcMET_corSumEtUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_uncorrectedPtUncorrJES", &PatMET_allLayer1METstcMET_uncorrectedPtUncorrJES, &b_PatMET_allLayer1METstcMET_uncorrectedPtUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_uncorrectedPhiUncorrJES", &PatMET_allLayer1METstcMET_uncorrectedPhiUncorrJES, &b_PatMET_allLayer1METstcMET_uncorrectedPhiUncorrJES);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_corExUncorrMUON", &PatMET_allLayer1METstcMET_corExUncorrMUON, &b_PatMET_allLayer1METstcMET_corExUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_corEyUncorrMUON", &PatMET_allLayer1METstcMET_corEyUncorrMUON, &b_PatMET_allLayer1METstcMET_corEyUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_corSumEtUncorrMUON", &PatMET_allLayer1METstcMET_corSumEtUncorrMUON, &b_PatMET_allLayer1METstcMET_corSumEtUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_uncorrectedPtUncorrMUON", &PatMET_allLayer1METstcMET_uncorrectedPtUncorrMUON, &b_PatMET_allLayer1METstcMET_uncorrectedPtUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_uncorrectedPhiUncorrMUON", &PatMET_allLayer1METstcMET_uncorrectedPhiUncorrMUON, &b_PatMET_allLayer1METstcMET_uncorrectedPhiUncorrMUON);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_isCaloMET", &PatMET_allLayer1METstcMET_isCaloMET, &b_PatMET_allLayer1METstcMET_isCaloMET);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_isRecoMET", &PatMET_allLayer1METstcMET_isRecoMET, &b_PatMET_allLayer1METstcMET_isRecoMET);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_maxEtInEmTowers", &PatMET_allLayer1METstcMET_maxEtInEmTowers, &b_PatMET_allLayer1METstcMET_maxEtInEmTowers);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_maxEtInHadTowers", &PatMET_allLayer1METstcMET_maxEtInHadTowers, &b_PatMET_allLayer1METstcMET_maxEtInHadTowers);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_etFractionHadronic", &PatMET_allLayer1METstcMET_etFractionHadronic, &b_PatMET_allLayer1METstcMET_etFractionHadronic);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_emEtFraction", &PatMET_allLayer1METstcMET_emEtFraction, &b_PatMET_allLayer1METstcMET_emEtFraction);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_hadEtInHB", &PatMET_allLayer1METstcMET_hadEtInHB, &b_PatMET_allLayer1METstcMET_hadEtInHB);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_hadEtInHO", &PatMET_allLayer1METstcMET_hadEtInHO, &b_PatMET_allLayer1METstcMET_hadEtInHO);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_hadEtInHE", &PatMET_allLayer1METstcMET_hadEtInHE, &b_PatMET_allLayer1METstcMET_hadEtInHE);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_hadEtInHF", &PatMET_allLayer1METstcMET_hadEtInHF, &b_PatMET_allLayer1METstcMET_hadEtInHF);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_emEtInEB", &PatMET_allLayer1METstcMET_emEtInEB, &b_PatMET_allLayer1METstcMET_emEtInEB);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_emEtInEE", &PatMET_allLayer1METstcMET_emEtInEE, &b_PatMET_allLayer1METstcMET_emEtInEE);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_emEtInHF", &PatMET_allLayer1METstcMET_emEtInHF, &b_PatMET_allLayer1METstcMET_emEtInHF);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_metSignificance", &PatMET_allLayer1METstcMET_metSignificance, &b_PatMET_allLayer1METstcMET_metSignificance);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_CaloSETInpHF", &PatMET_allLayer1METstcMET_CaloSETInpHF, &b_PatMET_allLayer1METstcMET_CaloSETInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_CaloSETInmHF", &PatMET_allLayer1METstcMET_CaloSETInmHF, &b_PatMET_allLayer1METstcMET_CaloSETInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_CaloMETInpHF", &PatMET_allLayer1METstcMET_CaloMETInpHF, &b_PatMET_allLayer1METstcMET_CaloMETInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_CaloMETInmHF", &PatMET_allLayer1METstcMET_CaloMETInmHF, &b_PatMET_allLayer1METstcMET_CaloMETInmHF);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_CaloMETPhiInpHF", &PatMET_allLayer1METstcMET_CaloMETPhiInpHF, &b_PatMET_allLayer1METstcMET_CaloMETPhiInpHF);
   fChain->SetBranchAddress("PatMET_allLayer1METstcMET_CaloMETPhiInmHF", &PatMET_allLayer1METstcMET_CaloMETPhiInmHF, &b_PatMET_allLayer1METstcMET_CaloMETPhiInmHF);
   fChain->SetBranchAddress("CaloTower_towerMaker_nCaloTower", &CaloTower_towerMaker_nCaloTower, &b_CaloTower_towerMaker_nCaloTower);
   fChain->SetBranchAddress("CaloTower_towerMaker_emEnergy", CaloTower_towerMaker_emEnergy, &b_CaloTower_towerMaker_emEnergy);
   fChain->SetBranchAddress("CaloTower_towerMaker_hadEnergy", CaloTower_towerMaker_hadEnergy, &b_CaloTower_towerMaker_hadEnergy);
   fChain->SetBranchAddress("CaloTower_towerMaker_outerEnergy", CaloTower_towerMaker_outerEnergy, &b_CaloTower_towerMaker_outerEnergy);
   fChain->SetBranchAddress("CaloTower_towerMaker_emEt", CaloTower_towerMaker_emEt, &b_CaloTower_towerMaker_emEt);
   fChain->SetBranchAddress("CaloTower_towerMaker_hadEt", CaloTower_towerMaker_hadEt, &b_CaloTower_towerMaker_hadEt);
   fChain->SetBranchAddress("CaloTower_towerMaker_outerEt", CaloTower_towerMaker_outerEt, &b_CaloTower_towerMaker_outerEt);
   fChain->SetBranchAddress("CaloTower_towerMaker_emPositionX", CaloTower_towerMaker_emPositionX, &b_CaloTower_towerMaker_emPositionX);
   fChain->SetBranchAddress("CaloTower_towerMaker_emPositionY", CaloTower_towerMaker_emPositionY, &b_CaloTower_towerMaker_emPositionY);
   fChain->SetBranchAddress("CaloTower_towerMaker_emPositionZ", CaloTower_towerMaker_emPositionZ, &b_CaloTower_towerMaker_emPositionZ);
   fChain->SetBranchAddress("CaloTower_towerMaker_hadPositionX", CaloTower_towerMaker_hadPositionX, &b_CaloTower_towerMaker_hadPositionX);
   fChain->SetBranchAddress("CaloTower_towerMaker_hadPositionY", CaloTower_towerMaker_hadPositionY, &b_CaloTower_towerMaker_hadPositionY);
   fChain->SetBranchAddress("CaloTower_towerMaker_hadPositionZ", CaloTower_towerMaker_hadPositionZ, &b_CaloTower_towerMaker_hadPositionZ);
   fChain->SetBranchAddress("CaloTower_towerMaker_emLvl1", CaloTower_towerMaker_emLvl1, &b_CaloTower_towerMaker_emLvl1);
   fChain->SetBranchAddress("CaloTower_towerMaker_hadLvl1", CaloTower_towerMaker_hadLvl1, &b_CaloTower_towerMaker_hadLvl1);
   fChain->SetBranchAddress("CaloTower_towerMaker_hadEnergyHeOuterLayer", CaloTower_towerMaker_hadEnergyHeOuterLayer, &b_CaloTower_towerMaker_hadEnergyHeOuterLayer);
   fChain->SetBranchAddress("CaloTower_towerMaker_hadEnergyHeInnerLayer", CaloTower_towerMaker_hadEnergyHeInnerLayer, &b_CaloTower_towerMaker_hadEnergyHeInnerLayer);
   fChain->SetBranchAddress("CaloTower_towerMaker_ecalTime", CaloTower_towerMaker_ecalTime, &b_CaloTower_towerMaker_ecalTime);
   fChain->SetBranchAddress("CaloTower_towerMaker_hcalTime", CaloTower_towerMaker_hcalTime, &b_CaloTower_towerMaker_hcalTime);
   fChain->SetBranchAddress("CaloTower_towerMaker_ieta", CaloTower_towerMaker_ieta, &b_CaloTower_towerMaker_ieta);
   fChain->SetBranchAddress("CaloTower_towerMaker_iphi", CaloTower_towerMaker_iphi, &b_CaloTower_towerMaker_iphi);
   fChain->SetBranchAddress("CaloTower_towerMaker_numCrystals", CaloTower_towerMaker_numCrystals, &b_CaloTower_towerMaker_numCrystals);
   fChain->SetBranchAddress("CaloTower_towerMaker_px", CaloTower_towerMaker_px, &b_CaloTower_towerMaker_px);
   fChain->SetBranchAddress("CaloTower_towerMaker_py", CaloTower_towerMaker_py, &b_CaloTower_towerMaker_py);
   fChain->SetBranchAddress("CaloTower_towerMaker_pz", CaloTower_towerMaker_pz, &b_CaloTower_towerMaker_pz);
   fChain->SetBranchAddress("CaloTower_towerMaker_pt", CaloTower_towerMaker_pt, &b_CaloTower_towerMaker_pt);
   fChain->SetBranchAddress("CaloTower_towerMaker_p", CaloTower_towerMaker_p, &b_CaloTower_towerMaker_p);
   fChain->SetBranchAddress("CaloTower_towerMaker_energy", CaloTower_towerMaker_energy, &b_CaloTower_towerMaker_energy);
   fChain->SetBranchAddress("CaloTower_towerMaker_phi", CaloTower_towerMaker_phi, &b_CaloTower_towerMaker_phi);
   fChain->SetBranchAddress("CaloTower_towerMaker_eta", CaloTower_towerMaker_eta, &b_CaloTower_towerMaker_eta);
   fChain->SetBranchAddress("CaloTower_towerMaker_theta", CaloTower_towerMaker_theta, &b_CaloTower_towerMaker_theta);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleEle10_LW_OnlyPixelM_L1R", &TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleEle10_LW_OnlyPixelM_L1R, &b_TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleEle10_LW_OnlyPixelM_L1R);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleEle10_LW_OnlyPixelM_L1R", &TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleEle10_LW_OnlyPixelM_L1R, &b_TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleEle10_LW_OnlyPixelM_L1R);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleEle10_LW_OnlyPixelM_L1R", &TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleEle10_LW_OnlyPixelM_L1R, &b_TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleEle10_LW_OnlyPixelM_L1R);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleEle10_LW_OnlyPixelM_L1R", &TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleEle10_LW_OnlyPixelM_L1R, &b_TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleEle10_LW_OnlyPixelM_L1R);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleEle10_LW_OnlyPixelM_L1R", &TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleEle10_LW_OnlyPixelM_L1R, &b_TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleEle10_LW_OnlyPixelM_L1R);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleMu3", &TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleMu3, &b_TriggerResults_TriggerResultsHLT_HLTPathIsValid_HLT_DoubleMu3);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleMu3", &TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleMu3, &b_TriggerResults_TriggerResultsHLT_HLTPathAccept_HLT_DoubleMu3);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleMu3", &TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleMu3, &b_TriggerResults_TriggerResultsHLT_HLTPathError_HLT_DoubleMu3);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleMu3", &TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleMu3, &b_TriggerResults_TriggerResultsHLT_HLTPathWasRun_HLT_DoubleMu3);
   fChain->SetBranchAddress("TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleMu3", &TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleMu3, &b_TriggerResults_TriggerResultsHLT_HLTPathIndex_HLT_DoubleMu3);
   Notify();
}

Bool_t AnalysisTreeClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalysisTreeClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalysisTreeClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalysisTreeClass_cxx
