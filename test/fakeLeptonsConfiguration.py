from cp3_llbb.Framework.CmdLine import CmdLine

class TTWFakeCmdLine(CmdLine):
    def __init__(self, override=None, defaults=None):
        super(TTWFakeCmdLine, self).__init__(override=override, defaults=defaults)
    def registerOptions(self):
        super(TTWFakeCmdLine, self).registerOptions()

        from FWCore.ParameterSet.VarParsing import VarParsing
        self.options.register('localTest'
                , False, VarParsing.multiplicity.singleton, VarParsing.varType.bool
                , 'Run a local test (input file defined in the configuration file)'
                )

options = TTWFakeCmdLine()

if options.runOnData:
    options.changeDefaults(globalTag="80X_dataRun2_2016SeptRepro_v7", process="RECO")
else:
    options.changeDefaults(globalTag="80X_mcRun2_asymptotic_2016_TrancheIV_v8", process="PAT")

import FWCore.ParameterSet.Config as cms

from cp3_llbb.Framework import Framework

from cp3_llbb.TTWAnalysis.Configuration import addTTWCategories, addTTWAnalyzer, addTTWCandidatesAnalyzer, customizeProducers

framework = Framework.Framework(options)
framework.process.framework.compressionSettings = cms.untracked.int32(207)

framework.process.fakeHltFilter = cms.EDFilter("TriggerResultsFilter",
        triggerConditions=cms.vstring("HLT_Mu17_v* OR HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*"),
        hltResults=cms.InputTag("TriggerResults", "", "HLT"),
        l1tResults=cms.InputTag(""),#cms.InputTag("gtStage2Digis"),
        l1tIgnoreMaskAndPrescale=cms.bool(False),
        throw=cms.bool(True),
        )

framework.addPreFilter(framework.process.fakeHltFilter)
##framework.process.framework.filters.hlt = cms.PSet(
##        type=cms.string("trigger"),
##        enable=cms.bool(True),
##        parameters=cms.PSet(
##            triggerSelection=cms.string("HLT_Mu17_v* OR HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*"),
##            triggerConfiguration=cms.PSet(
##                hltResults=cms.InputTag("TriggerResults", "", "HLT"),
##                l1tResults=cms.InputTag(""),#cms.InputTag("gtStage2Digis"),
##                l1tIgnoreMaskAndPrescale=cms.bool(False),
##                throw=cms.bool(True),
##                )
##            )
##        )

## ANALYZERS
addTTWCategories        (framework, addDilepton=False, addNonPrompt=True)
addTTWAnalyzer          (framework) ## make candidates
addTTWCandidatesAnalyzer(framework, addCombinations=False) ## fill most branches

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

from cp3_llbb.Framework.JetsProducer import discriminators_deepFlavour
framework.redoJEC(addBtagDiscriminators=discriminators_deepFlavour)

framework.applyMuonCorrection("rochester")
framework.applyElectronRegression()
framework.applyElectronSmearing()

##if not options.runOnData:
##    framework.smearJets(resolutionFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt', scaleFactorFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_SF_AK4PFchs.txt')
##    framework.doSystematics(['jec', 'jer'], jec={'uncertaintiesFile': 'cp3_llbb/TTWAnalysis/data/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt', 'splitBySources': True})

process = framework.create()

if options.localTest:
    if options.runOnData:
        process.source.fileNames = cms.untracked.vstring(
              "/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/5673D8AE-21EB-E611-8CEE-002590E3A0FA.root"
            , "/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/E05CF4CC-23EB-E611-B2F8-0025900EAB5E.root"
            , "/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/8228B3C8-23EB-E611-AEEC-002590E39DF4.root"
            , "/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/58B756C5-23EB-E611-A6F3-0CC47A13CD56.root"
            , "/store/data/Run2016G/SingleElectron/MINIAOD/03Feb2017-v1/50000/32C03BC6-23EB-E611-91BC-002590E3A286.root"
            , "/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/7843236A-BBEA-E611-8E2C-001E675A690A.root"
            , "/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/829FAA6D-BBEA-E611-81E6-001E67A4061D.root"
            , "/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/04B06F6B-BBEA-E611-A113-001E67A3F92F.root"
            , "/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/6274027F-BBEA-E611-9757-001E67D80528.root"
            , "/store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/881A217B-BBEA-E611-9EB1-001E67D80528.root"
            )
        process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

    else: ## MC

        ## process.source.fileNames = cms.untracked.vstring(
        ##     '/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/04D51FB4-B2B8-E511-A399-047D7B881D6A.root'
        ##     )
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

        ##process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

        ## ttH objects sync file & debug options
        ## process.source.fileNames = cms.untracked.vstring(
        ##     "/store/mc/RunIISummer16MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/0015BB42-9BAA-E611-8C7F-0CC47A7E0196.root" ## ttH-multilepton sync sample
        ##     )

        ## job1...
        process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
        process.source.fileNames = cms.untracked.vstring([
              "/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/163E57C9-7ABE-E611-A73A-0025905B857E.root",
              "/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/B48227FC-7BBE-E611-9699-0CC47A4C8F26.root",
              "/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/D6593DF9-83BE-E611-B0C9-0025905B857E.root",
              "/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/380E12A9-84BE-E611-8D12-0025905A612C.root",
              "/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/06B14BE1-88BE-E611-BF9E-0025905A608C.root",
            ])
        ## process.source.skipEvents = cms.untracked.uint32(16952)
        ## process.MessageLogger.cerr.FwkReport.reportEvery = 1

    ##process.MessageLogger = cms.Service(
    ##    "MessageLogger",
    ##    destinations = cms.untracked.vstring(
    ##        "cout"
    ##        ),
    ##    cout = cms.untracked.PSet(
    ##        threshold = cms.untracked.string("DEBUG")
    ##        ),
    ##    debugModules = cms.untracked.vstring("framework"),
    ##    categories=cms.untracked.vstring("Framework")
    ##    )
