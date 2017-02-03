#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "cp3_llbb/Framework/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"

namespace TTWAnalysis {
/**
 * Additional electron identification variables
 */
class DictElectronIDVars : public DictTool<pat::Electron> {
public:
  DictElectronIDVars(const edm::ParameterSet& config)
    : DictTool<pat::Electron>(config)
  {}
  virtual ~DictElectronIDVars() {}

  virtual Dict evaluate(edm::Ptr<pat::Electron> el,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    const bool valid{el.isNonnull() && el->originalObjectRef().isNonnull()};
    const double eta = valid ? el->superCluster()->eta() : 0.;

    Dict ret;
    ret.add("isEB", valid ? el->isEB() : false);
    ret.add("isEE", valid ? el->isEE() : false);

    // electron trigger emulation
    //     https://github.com/peruzzim/cmgtools-lite/blob/76X_for2016basis/TTHAnalysis/python/tools/functionsTTH.py#L10-L20
    bool passTrigEmu{false};
    if ( valid ) {
      double eInvMinusPInv = el->ecalEnergy() > 0. ? (1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) : 9.e9;
      if ( std::abs(eta) < 1.479 ) {
        passTrigEmu = (
              ( el->hadronicOverEm() < 0.10 )
           && ( std::abs(el->deltaEtaSuperClusterTrackAtVtx()) < 0.01 )
           && ( std::abs(el->deltaPhiSuperClusterTrackAtVtx()) < 0.04 )
           && ( -0.05 < eInvMinusPInv ) && ( eInvMinusPInv < 0.01 )
           && ( el->full5x5_sigmaIetaIeta() < 0.011 )
           );
      } else {
        passTrigEmu = (
              ( el->hadronicOverEm() < 0.07 )
           && ( std::abs(el->deltaEtaSuperClusterTrackAtVtx()) < 0.008 )
           && ( std::abs(el->deltaPhiSuperClusterTrackAtVtx()) < 0.07 )
           && ( -0.05 < eInvMinusPInv ) && ( eInvMinusPInv < 0.005 )
           && ( el->full5x5_sigmaIetaIeta() < 0.030 )
           );
      }
    }
    ret.add("TrigEmu", passTrigEmu);

    ret.add("NMissInner", valid && el->gsfTrack().isNonnull() ? el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) : 0 );
    ret.add("passConversionVeto", valid ? el->passConversionVeto() : false);
    ret.add("ThreeChargeAgreement", valid ? el->isGsfCtfScPixChargeConsistent() : false); // ??

    return ret;
  }
};

/**
 * Additional muon identification variables
 */
class DictMuonIDVars : public DictTool<pat::Muon> {
public:
  DictMuonIDVars(const edm::ParameterSet& config)
    : DictTool<pat::Muon>(config)
  {}
  virtual ~DictMuonIDVars() {}

  virtual Dict evaluate(edm::Ptr<pat::Muon> mu,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    const bool valid{mu.isNonnull() && mu->originalObjectRef().isNonnull()};

    Dict ret;
    ret.add("DPToPT", valid && mu->innerTrack().isNonnull() ? mu->innerTrack()->ptError()/mu->innerTrack()->pt() : 0);

    return ret;
  }
};

/**
 * Additional jet identification variables
 */
class DictJetIDVars : public DictTool<pat::Jet> {
public:
  DictJetIDVars(const edm::ParameterSet& config)
    : DictTool<pat::Jet>(config)
  {
    const auto& bTagCfg = config.getParameter<edm::ParameterSet>("bTaggers");
    for ( const auto& bTagNm : bTagCfg.getParameterNames() ) {
      m_bTaggers.emplace(bTagNm, bTagCfg.getParameter<std::string>(bTagNm));
    }
  }
  virtual ~DictJetIDVars() {}

  virtual Dict evaluate(edm::Ptr<pat::Jet> jet,
      const edm::Event* /**/, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override
  {
    Dict ret;
    const bool valid{jet.isNonnull() && jet->originalObjectRef().isNonnull()};
    ret.add("jecFactor"   , valid ? jet->jecFactor(0) : 0);
    ret.add("area"        , valid ? jet->jetArea() : -1.);
    ret.add("partonFlavor", valid ? jet->partonFlavour() : -10);
    ret.add("hadronFlavor", valid ? jet->hadronFlavour() : -10);
    ret.add("puJetID"     , valid ? jet->userFloat("pileupJetId:fullDiscriminant") : 0);
    for ( const auto& varNm : m_bTaggers ) {
      ret.add(varNm.first, valid ? jet->bDiscriminator(varNm.second) : 0);
    }
    return ret;
  }
private:
  std::map<std::string,std::string> m_bTaggers;
};
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictElectronIDVars, "ttw_electronIDVars");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictMuonIDVars    , "ttw_muonIDVars");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Jet     >::factory, TTWAnalysis::DictJetIDVars     , "ttw_jetIDVars");
