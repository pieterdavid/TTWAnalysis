#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "cp3_llbb/Framework/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"
#include "cp3_llbb/TTWAnalysis/interface/NewTypes.h"
#include "cp3_llbb/TTWAnalysis/interface/stl_helpers.h"

#include <Math/VectorUtil.h>
using namespace ROOT::Math;

namespace TTWAnalysis {

class DictElectronL1Vars : public DictTool<pat::Electron> {
public:
  DictElectronL1Vars(const edm::ParameterSet& config)
    : DictTool<pat::Electron>(config)
    , m_l1tDRCut( config.getUntrackedParameter<double>("L1DRCut", std::numeric_limits<float>::max()) )
    , m_l1tDPtCut( config.getUntrackedParameter<double>("L1DPtCut", std::numeric_limits<float>::max()) )
  {}
  virtual ~DictElectronL1Vars() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    m_l1EGammasToken = collector.consumes<l1t::EGammaBxCollection>(config.getParameter<edm::InputTag>("L1EGammaInputTag"));
  }

  virtual Dict evaluate(edm::Ptr<pat::Electron> el,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  edm::EDGetTokenT<l1t::EGammaBxCollection> m_l1EGammasToken;
  const float m_l1tDRCut, m_l1tDPtCut;
};

class DictMuonL1Vars : public DictTool<pat::Muon> {
public:
  DictMuonL1Vars(const edm::ParameterSet& config)
    : DictTool<pat::Muon>(config)
    , m_l1tDRCut( config.getUntrackedParameter<double>("L1DRCut", std::numeric_limits<float>::max()) )
    , m_l1tDPtCut( config.getUntrackedParameter<double>("L1DPtCut", std::numeric_limits<float>::max()) )
  {}
  virtual ~DictMuonL1Vars() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    m_l1MuonsToken = collector.consumes<l1t::MuonBxCollection>(config.getParameter<edm::InputTag>("L1MuonInputTag"));
  }

  virtual Dict evaluate(edm::Ptr<pat::Muon> mu,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  edm::EDGetTokenT<l1t::MuonBxCollection> m_l1MuonsToken;
  const float m_l1tDRCut, m_l1tDPtCut;
};

Dict DictElectronL1Vars::evaluate(edm::Ptr<pat::Electron> el,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  const bool valid{el.isNonnull() && el->superCluster().isNonnull()};
  // const double eta = valid ? el->superCluster()->eta() : 0.;

  std::vector<const l1t::EGamma*> l1Matches;
  const l1t::EGamma* bestL1Match{nullptr};
  if ( event ) {
    edm::Handle<l1t::EGammaBxCollection> l1EGammas;
    event->getByToken(m_l1EGammasToken, l1EGammas);

    transform_if(l1EGammas->begin(0), l1EGammas->end(0), std::back_inserter(l1Matches)
        , [el,this] ( const l1t::EGamma& l1EG ) -> bool
          { return (( VectorUtil::DeltaR(el->p4(), l1EG.p4()) <= m_l1tDRCut )
         && ( std::abs(el->p4().Pt() - l1EG.p4().Pt()) / el->p4().Pt() <= m_l1tDPtCut ) ); }
        , [] ( const l1t::EGamma& l1EG ) -> const l1t::EGamma* { return &l1EG; }
        );
  }

  if ( valid && ( ! l1Matches.empty() ) ) {
    /* std::cout << "Found " << l1Matches.size() << " matching L1 EG seeds" << std::endl;
    std::cout << "Offline : " << el->p4() << std::endl;
    for ( const l1t::EGamma* l1Match : l1Matches ) {
      std::cout << "L1T     - " << l1Match->p4() << std::endl;
    } */
    std::stable_sort(std::begin(l1Matches), std::end(l1Matches), LessOf<const l1t::EGamma*>(
          [el] ( const l1t::EGamma* l1El ) { return VectorUtil::DeltaR(el->p4(), l1El->p4()); }));
    bestL1Match = l1Matches[0];
  }

  Dict ret;
  ret.add("L1_hwQual", valid && bestL1Match ? bestL1Match->hwQual() : -1);
  ret.add("L1_isoEt" , valid && bestL1Match ? bestL1Match->isoEt() : -1);

  return ret;
}

Dict DictMuonL1Vars::evaluate(edm::Ptr<pat::Muon> mu,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  const bool valid{mu.isNonnull()}; // && mu->innerTrack().isNonnull()

  std::vector<const l1t::Muon*> l1Matches;
  const l1t::Muon* bestL1Match{nullptr};
  if ( event ) {
    edm::Handle<l1t::MuonBxCollection> l1Muons;
    event->getByToken(m_l1MuonsToken, l1Muons);

    transform_if(l1Muons->begin(0), l1Muons->end(0), std::back_inserter(l1Matches)
        , [mu,this] ( const l1t::Muon& l1Mu ) -> bool
          { return (( VectorUtil::DeltaR(mu->p4(), l1Mu.p4()) <= m_l1tDRCut )
         && ( std::abs(mu->p4().Pt() - l1Mu.p4().Pt()) / mu->p4().Pt() <= m_l1tDPtCut ) ); }
        , [] ( const l1t::Muon& l1Mu ) -> const l1t::Muon* { return &l1Mu; }
        );
  }

  if ( valid && ( ! l1Matches.empty() ) ) {
    /* std::cout << "Found " << l1Matches.size() << " matching L1 mu seeds" << std::endl;
    std::cout << "Offline : " << mu->p4() << std::endl;
    for ( const l1t::Muon* l1Match : l1Matches ) {
      std::cout << "L1T     - " << l1Match->p4() << std::endl;
    } */
    std::stable_sort(std::begin(l1Matches), std::end(l1Matches), LessOf<const l1t::Muon*>(
          [mu] ( const l1t::Muon* l1Mu ) { return VectorUtil::DeltaR(mu->p4(), l1Mu->p4()); }));
    bestL1Match = l1Matches[0];
  }

  Dict ret;
  ret.add("L1_hwQual", valid && bestL1Match ? bestL1Match->hwQual() : -1);
  ret.add("L1_hwIsoSum", valid && bestL1Match ? bestL1Match->hwIsoSum() : -1);

  return ret;
}

}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictMuonL1Vars    , "ttw_muonL1Vars");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictElectronL1Vars, "ttw_electronL1Vars");
