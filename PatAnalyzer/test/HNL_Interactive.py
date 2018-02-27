import sys
import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts


skim     ="NOSkim"
isMC     = bool(True)
isMSUGRA = bool(False)
isSMS    = bool(False)
useCrab  = bool(False)
doSusyTopProjection = bool(False)
inputFile=""
outputFile="/afs/cern.ch/work/d/dezhu/CMSSW_8_0_30/src/0_Results/test.root"
def getVal(arg):
    i=0
    while i < len(arg) and arg[i] != "=": i+=1
    return arg[i+1:]

## loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i] 
    if "isMC" in sys.argv[i] :
        isMC=False if getVal(sys.argv[i]) == "False" else True
    elif "isMSUGRA" in sys.argv[i]:
        isMSUGRA=True
    elif "isSMS" in sys.argv[i]:
        isSMS=True
    elif "skim" in sys.argv[i]:
        skim=getVal(sys.argv[i])
    elif "output" in sys.argv[i]:
        outputFile=getVal(sys.argv[i])
    elif "input" in sys.argv[i] :
        inputFile=getVal(sys.argv[i])

        
if skim=="" : print "WARNING: No Skim Conditions have been provided \n"

#process = cms.Process("FakeLeptons")
process = cms.Process("pippo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'WARNING' # Options: INFO, WARNING, ERROR
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# Load the producer for MVA IDs. Make sure it is also added to the sequence!
#process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
#process.load("RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff")
#process.load("RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff")

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff', 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



if isMC:
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    from Configuration.AlCa.GlobalTag import GlobalTag
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '') #80X_mcRun2_asymptotic_2016_TrancheIV_v8
else:
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
    process.GlobalTag.globaltag="80X_dataRun2_2016SeptRepro_v7" # 

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),                    
                             fileNames = cms.untracked.vstring(),      
                             )

#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

# Input MiniAOD destination
if isMC:
    process.source.fileNames = cms.untracked.vstring(               
       'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_104.root',
        #'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_16.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_17.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_21.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_26.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_43.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_66.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_7.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_83.root',

        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_mu_onshell_pre2017_NLO/heavyNeutrino_104.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_mu_onshell_pre2017_NLO/heavyNeutrino_16.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_mu_onshell_pre2017_NLO/heavyNeutrino_21.root',
        # 'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00316227766017_mu_onshell_pre2017_nlo/heavyneutrino_2.root',
        # 'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00316227766017_mu_onshell_pre2017_nlo/heavyneutrino_32.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00316227766017_mu_onshell_pre2017_nlo/heavyneutrino_43.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00316227766017_mu_onshell_pre2017_nlo/heavyneutrino_66.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00316227766017_mu_onshell_pre2017_nlo/heavyneutrino_7.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00316227766017_mu_onshell_pre2017_nlo/heavyneutrino_83.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00316227766017_mu_onshell_pre2017_nlo/heavyneutrino_87.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00316227766017_mu_onshell_pre2017_nlo/heavyneutrino_92.root',

        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_104.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_16.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_17.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_26.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_2.root',
        #'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_32.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_37.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_46.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_66.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_7.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_83.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_87.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_e_onshell_pre2017_nlo/heavyneutrino_92.root',

        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_mu_onshell_pre2017_nlo/heavyneutrino_104.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_mu_onshell_pre2017_nlo/heavyneutrino_14.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_mu_onshell_pre2017_nlo/heavyneutrino_16.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_mu_onshell_pre2017_nlo/heavyneutrino_17.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_mu_onshell_pre2017_nlo/heavyneutrino_26.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyneutrinominiaod/moriond17/displaced/heavyneutrino_trilepton_m-2.1_v-0.00244948974278_mu_onshell_pre2017_nlo/heavyneutrino_2.root',
        # 'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00244948974278_mu_onshell_pre2017_NLO/heavyNeutrino_32.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00244948974278_mu_onshell_pre2017_NLO/heavyNeutrino_43.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00244948974278_mu_onshell_pre2017_NLO/heavyNeutrino_66.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00244948974278_mu_onshell_pre2017_NLO/heavyNeutrino_7.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00244948974278_mu_onshell_pre2017_NLO/heavyNeutrino_83.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00244948974278_mu_onshell_pre2017_NLO/heavyNeutrino_87.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00244948974278_mu_onshell_pre2017_NLO/heavyNeutrino_92.root',

        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_e_onshell_pre2017_NLO/heavyNeutrino_14.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_e_onshell_pre2017_NLO/heavyNeutrino_32.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_e_onshell_pre2017_NLO/heavyNeutrino_43.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_e_onshell_pre2017_NLO/heavyNeutrino_7.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_e_onshell_pre2017_NLO/heavyNeutrino_83.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_e_onshell_pre2017_NLO/heavyNeutrino_87.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_e_onshell_pre2017_NLO/heavyNeutrino_92.root',

        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_104.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_16.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_21.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_2.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_43.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_66.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_7.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_83.root',
        'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_87.root',
        # 'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.01_mu_onshell_pre2017_NLO/heavyNeutrino_92.root',

        #'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_e_onshell_pre2017_NLO/heavyNeutrino_104.root',
        # 'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_trilepton_M-2.1_V-0.00316227766017_mu_onshell_pre2017_NLO/heavyNeutrino_104.root'
        #'file:/afs/cern.ch/user/d/dezhu/workspace/public/10A42D2B-2274-E711-AA76-047D7B2E84EC.root' #10/0.0001, ctau = 176.5 mm
       #'file:/afs/cern.ch/user/d/dezhu/workspace/public/20989B66-CA73-E711-ABD7-7845C4FC3C86.root' #5/0.001, ctau = 56.8 mm
#             'file:/cms/data/store/mc/RunIISpring16MiniAODv2/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/8E8FD1DC-C526-E611-A48A-002590E7DFFC.root'
    )

else:
   process.source.fileNames = cms.untracked.vstring( 
#         'file:/cms/data/store/data/Run2016B/DoubleEG/MINIAOD/PromptReco-v2/000/273/158/00000/32EBF0E5-E619-E611-9903-02163E014118.root'
          'file:pickevents.root'
    )
#    runOnData(process)
   

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )

## input files
#for i in inputFile.split(","):
#    print "Adding: ", i
#    process.source.fileNames.append(i)

#FakeLeptons
process.FakeLeptons = cms.EDAnalyzer("HNL",
                                       MuonLabel = cms.InputTag("slimmedMuons"),
                                       displacedGlobalMuonsLabel = cms.InputTag("displacedGlobalMuons"),
                                       displacedStandAloneMuonsLabel = cms.InputTag("displacedStandAloneMuons"),
                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
                                       TauLabel = cms.InputTag("slimmedTaus"),
                                       convLabel = cms.InputTag("reducedEgamma:reducedConversions"),
                                       #TauDiscriminatorLabel = cms.InputTag("recoPFTauDiscriminator"),
                                       L1httLabel = cms.InputTag("l1extraParticles:MHT:RECO"),
                                       JetLabel = cms.InputTag("slimmedJets"),
                                       BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                       vtxLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                       rhoLabel = cms.InputTag("fixedGridRhoFastjetAll"),
                                       rhoLabelCN = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                                       PUInfoLabel = cms.InputTag("slimmedAddPileupInfo"),
                                       generatorLabel = cms.InputTag("generator"),
                                       lheevent = cms.InputTag("externalLHEProducer"),
                                       pfcLabel = cms.InputTag("packedPFCandidates"),
                                       HLTResultsLabel = cms.InputTag("TriggerResults::HLT"),
                                       filterResultsLabel = cms.InputTag("TriggerResults::RECO"),
                                       #prescales = cms.InputTag("patTrigger"),
                                       prescales = cms.InputTag("patTrigger"),
                                       METLabel = cms.InputTag("slimmedMETs"),
                                       genPartsLabel = cms.InputTag("prunedGenParticles"),
                                       #METFilter = cms.InputTag("TriggerResults::PAT"),
                                       qualityCuts = PFTauQualityCuts,
                                       SampleLabel = cms.untracked.string("ElectronsMC"), # defines a piece of code to run; helps to avoid code recompilation
                                       mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                                       mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
                                       mvaValuesMapHZZ = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
                                       mvaCategoriesMapHZZ = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Categories"),
                                       #mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                       #mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
                                       BadChCandFilter = cms.InputTag("BadChargedCandidateFilter"),
                                       BadPFMuon = cms.InputTag("BadPFMuonFilter")
                                       )

process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter" # checks for fake PVs automatically
                                                  , filterParams =cms.PSet(
                                                                           minNdof = cms.double( 4. )
                                                                           , maxZ    = cms.double( 24. )
                                                                           , maxRho  = cms.double( 2. ) )		
                                                  , filter       = cms.bool( False ) # use only as producer
                                                  , src          = cms.InputTag( 'offlineSlimmedPrimaryVertices' )
)


##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
   reverseDecision = cms.bool(False)
)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
## for miniAOD running
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')

## for miniAOD running
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")


if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
       *process.BadChargedCandidateFilter
       *process.BadPFMuonFilter
#    *process.electronMVAValueMapProducer
    *process.egmGsfElectronIDSequence
    *process.FakeLeptons
        )
else:
   import FWCore.PythonUtilities.LumiList as LumiList
   process.source.lumisToProcess = LumiList.LumiList(filename = '/cms/data/store/user/t2/users/lesya/CMSSW_8_0_20_patch1/src/SUSYAnalyzer/PatAnalyzer/test/JSON/Cert_271036-282037_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt').getVLuminosityBlockRange()
   process.p = cms.Path(
	process.goodOfflinePrimaryVertices
#       *process.eeBadScFilter
       *process.HBHENoiseFilterResultProducer
       *process.ApplyBaselineHBHENoiseFilter
       *process.ApplyBaselineHBHEIsoNoiseFilter
       *process.BadChargedCandidateFilter
       *process.BadPFMuonFilter
#       *process.electronMVAValueMapProducer
       *process.egmGsfElectronIDSequence
       *process.FakeLeptons
   )
    

 
