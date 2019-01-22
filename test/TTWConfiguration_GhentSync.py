from cp3_llbb.Framework.CmdLine import CmdLine

class TTWCmdLine(CmdLine):
    def __init__(self, override=None, defaults=None):
        super(TTWCmdLine, self).__init__(override=override, defaults=defaults)
    def registerOptions(self):
        super(TTWCmdLine, self).registerOptions()

        from FWCore.ParameterSet.VarParsing import VarParsing
        self.options.register('noCuts'
                , False, VarParsing.multiplicity.singleton, VarParsing.varType.bool
                , 'Turn off all selections (True) or not (False, default)'
                )
        self.options.register('localTest'
                , False, VarParsing.multiplicity.singleton, VarParsing.varType.bool
                , 'Run a local test (input file defined in the configuration file)'
                )

options = TTWCmdLine()
options._ensureParsed()
era = options.options.era ## string

if options.runOnData: ## these are for 
    if era == "2016":
        options.changeDefaults(globalTag="94X_dataRun2_v10", process="PAT")
    elif era == "2017":
        options.changeDefaults(globalTag="94X_dataRun2_v11", process="PAT")
    elif era == "2018":
        options.changeDefaults(globalTag="102X_dataRun2_Sep2018Rereco_v1", process="PAT")
else:
    if era == "2016":
        options.changeDefaults(globalTag="94X_mcRun2_asymptotic_v3", process="PAT")
    elif era == "2017":
        options.changeDefaults(globalTag="94X_mc2017_realistic_v17", process="PAT")
    elif era == "2018":
        options.changeDefaults(globalTag="102X_upgrade2018_realistic_v12", process="PAT")

import FWCore.ParameterSet.Config as cms

from cp3_llbb.Framework import Framework

from cp3_llbb.TTWAnalysis.Configuration_Ghent import addTTWCategories, addTTWAnalyzer, addTTWCandidatesAnalyzer, customizeProducers

framework = Framework.Framework(options)
framework.process.framework.compressionSettings = cms.untracked.int32(207)

## ANALYZERS
addTTWCategories        (framework, applyFilter=(not options.noCuts), addDilepton=True, addNonPrompt=False)
addTTWAnalyzer          (framework) ## make candidates
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
#from cp3_llbb.Framework.JetsProducer import discriminators_deepFlavour
#framework.redoJEC(addBtagDiscriminators=discriminators_deepFlavour)

#framework.applyMuonCorrection("rochester")
framework.applyElectronPostReco()
#framework.applyElectronRegression()
#framework.applyElectronSmearing()

## jets, not jetsSmeared in Ghent FWK
#if not options.runOnData:
#    framework.smearJets()

##    framework.smearJets(resolutionFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt', scaleFactorFile='cp3_llbb/Framework/data/Spring16_25nsV10_MC_SF_AK4PFchs.txt')
##    framework.doSystematics(['jec', 'jer'], jec={'uncertaintiesFile': 'cp3_llbb/TTWAnalysis/data/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt', 'splitBySources': True})

process = framework.create()

## Suggested commands:
# cmsRun TTWConfiguration_GhentSync.py runOnData=0 era=2016 globalTag=94X_mcRun2_asymptotic_v3 process=PAT localTest=true

if options.localTest:
    if options.runOnData:
        process.source.fileNames = cms.untracked.vstring(
            #"/store/data/Run2016F/MuonEG/MINIAOD/03Feb2017-v1/50000/0EEAA135-F8EA-E611-9F10-0CC47A7C35A8.root"
            #"/store/data/Run2016F/MuonEG/MINIAOD/03Feb2017-v1/50000/149ED47A-0CEB-E611-871F-0CC47A78A2F6.root"
            #"/store/data/Run2016F/MuonEG/MINIAOD/03Feb2017-v1/50000/0EEAA135-F8EA-E611-9F10-0CC47A7C35A8.root"
            #"/store/data/Run2016D/MuonEG/MINIAOD/03Feb2017-v1/80000/0A40F550-05EB-E611-B3DE-0CC47A4DEF68.root"
              "/store/data/Run2016F/DoubleEG/MINIAOD/03Feb2017-v1/80000/0006AFD8-F8EA-E611-9F9D-0CC47A13D09C.root"
            , "/store/data/Run2016F/DoubleEG/MINIAOD/03Feb2017-v1/80000/00A7E4D9-F8EA-E611-A62B-002590E3A004.root"
            , "/store/data/Run2016F/DoubleEG/MINIAOD/03Feb2017-v1/80000/022C3BB9-F0EA-E611-91DC-441EA171A998.root"
            , "/store/data/Run2016F/DoubleEG/MINIAOD/03Feb2017-v1/80000/027FD1D9-04EB-E611-872F-D8D385AE8ACA.root"
            , "/store/data/Run2016F/DoubleEG/MINIAOD/03Feb2017-v1/80000/02DEA755-FEEA-E611-94A1-0CC47AD98CFA.root"
            , "/store/data/Run2016F/DoubleEG/MINIAOD/03Feb2017-v1/80000/04846305-13EB-E611-A9FB-0CC47AD98CF2.root"
            )
        process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))

    else: ## MC

        ## job1...
        if era == "2016":
            process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
            process.source.fileNames = cms.untracked.vstring([
                   "/store/mc/RunIISummer16MiniAODv3/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/00000/B020C78B-11D1-E811-B413-0242AC130002.root"
                ])
        else:
            raise RuntimeError("ERA={0}".format(era))

        process.MessageLogger.cerr.FwkReport.reportEvery = 1

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
            #categories=cms.untracked.vstring("ttW-electronID"),
            categories=cms.untracked.vstring("ttW", "GhentLeptonMVA", "TTHLeptonMVA")
            )
