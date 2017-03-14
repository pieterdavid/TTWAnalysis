from cp3_llbb.Framework.CmdLine import CmdLine

options = CmdLine()

if options.runOnData:
    globalTag_ = "80X_dataRun2_2016SeptRepro_v7"
    processName_ = "RECO"
else:
    globalTag_ = "80X_mcRun2_asymptotic_2016_TrancheIV_v8"
    processName_ = "PAT"

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

from cp3_llbb.Framework import Framework

from cp3_llbb.TTWAnalysis.Configuration import addTTWAnalyzer, addTTWCandidatesAnalyzer, customizeProducers

framework = Framework.Framework(options.runOnData, eras.Run2_2016, globalTag=globalTag_, processName=processName_, hltProcessName=("HLT2" if not options.runOnData else "RECO"), verbose=True)

## ANALYZERS
addTTWAnalyzer          (framework, applyFilter=True) ## make candidates
addTTWCandidatesAnalyzer(framework) ## fill most branches

if not options.runOnData: ## MC truth analyzer
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

if options.runOnData:
    process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2015D/DoubleMuon/MINIAOD/16Dec2015-v1/10000/00039A2E-D7A7-E511-98EE-3417EBE64696.root'
        )
    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

else: ## MC

    process.source.fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/04D51FB4-B2B8-E511-A399-047D7B881D6A.root'
        )
    ## process.source.fileNames = cms.untracked.vstring(
    ##     '/store/mc/RunIIFall15MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/60000/14C51DB0-D6B8-E511-8D9B-8CDCD4A9A484.root'
    ##     )
    
    ## process.source.fileNames = cms.untracked.vstring(
    ##     '/store/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0ED0E1CB-90BF-E511-B379-0025905C4432.root'
    ##     )
    
    ## Tricky gen event from /store/mc/RunIISpring15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/00000/0014DC94-DC5C-E511-82FB-7845C4FC39F5.root
    ## First one is g g -> t tbar with one W -> bbar c
    ## Second is b bar -> t tbar semi-leptonic
    #process.source.eventsToProcess = cms.untracked.VEventRange(
    #        '1:52386:13083444',
    #        '1:34020:8496854'
    #        )
    
    ## Other tricky gen events, with lots of ISR
    ## From file:/nfs/scratch/fynu/swertz/CMSSW_7_4_15/src/cp3_llbb/TTAnalysis/test/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_miniAODv2_oneFile.root
    #process.source.eventsToProcess = cms.untracked.VEventRange(
    #        '1:321521:80300260',
    #        '1:357590:89308562',
    #        '1:387992:96901374'
    #        )
    
    #process.MessageLogger.cerr.FwkReport.reportEvery = 1

    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

    ## ttH objects sync file & debug options
    # process.source.fileNames = cms.untracked.vstring(
    #     "/store/mc/RunIISummer16MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/0015BB42-9BAA-E611-8C7F-0CC47A7E0196.root" ## ttH-multilepton sync sample
    #     )
    # 
    # process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
    # 
    # process.MessageLogger = cms.Service(
    #     "MessageLogger",
    #     destinations = cms.untracked.vstring(
    #         "detailedInfo",
    #         "critical"
    #         ),
    #     detailedInfo = cms.untracked.PSet(
    #         threshold = cms.untracked.string("DEBUG")
    #         ),
    #     debugModules = cms.untracked.vstring("framework"),
    #     categories=cms.untracked.vstring("ttW-electronID")
    #     )
