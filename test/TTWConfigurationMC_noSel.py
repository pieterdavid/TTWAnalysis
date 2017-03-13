import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from cp3_llbb.Framework import Framework

from cp3_llbb.TTWAnalysis.Configuration import addTTWAnalyzer, addTTWCandidatesAnalyzer, customizeProducers

globalTag_ = "80X_mcRun2_asymptotic_2016_TrancheIV_v8"
processName_ = 'PAT'

framework = Framework.Framework(False, eras.Run2_2016, globalTag=globalTag_, processName=processName_)

## ANALYZERS
addTTWAnalyzer          (framework, applyFilter=False) ## make candidates
addTTWCandidatesAnalyzer(framework)                    ## fill most branches

from cp3_llbb.TTWAnalysis.Configuration import lepton_WPs
framework.addAnalyzer('ttWTruth', cms.PSet(
        type = cms.string('ttwtruth_analyzer'),
        prefix = cms.string(''),
        enable = cms.bool(True),
        parameters = cms.PSet(
            TTWAnalyzer=cms.string("ttW"),
            JetWP=cms.vstring(*("ID{0}_Iso{1}".format(lID,lIso) for lID,lIso in lepton_WPs.iterkeys()))
            )
        ))

## PRODUCERS
customizeProducers(framework)

framework.redoJEC()
# framework.smearJets()
# framework.doSystematics(['jec', 'jer'])
framework.doSystematics([])

process = framework.create()

# process.source.fileNames = cms.untracked.vstring(
#     '/store/mc/RunIIFall15MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/14C51DB0-D6B8-E511-8D9B-8CDCD4A9A484.root'
#     )
process.source.fileNames = cms.untracked.vstring(
    "/store/mc/RunIISummer16MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/0015BB42-9BAA-E611-8C7F-0CC47A7E0196.root" ## ttH-multilepton sync sample
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring(
        "detailedInfo",
        "critical"
        ),
    detailedInfo = cms.untracked.PSet(
        threshold = cms.untracked.string("DEBUG")
        ),
    debugModules = cms.untracked.vstring("framework"),
    categories=cms.untracked.vstring("ttW-electronID")
    )
