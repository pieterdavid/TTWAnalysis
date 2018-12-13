#include "cp3_llbb/Framework/interface/METProducer.h"
#include "cp3_llbb/TTWAnalysis/interface/CandidatesProducer.h"
#include "cp3_llbb/TTWAnalysis/interface/TTWAnalyzer.h"
#include <Math/VectorUtil.h>

#include "TTWFakeLeptonCategory.h"

using namespace TTWAnalysis;

FakeLeptonCategory::~FakeLeptonCategory() {}

void FakeLeptonCategory::configure(const edm::ParameterSet& conf)
{
  m_ttwAnaName = conf.getParameter<std::string>("ttW");
  m_isEl = conf.getParameter<bool>("elOrMu");
  m_maxMET = conf.getParameter<double>("maxMET");
  m_maxMT  = conf.getParameter<double>("maxMET");
  m_selCat = conf.getParameter<std::string>("LeptonWP");
  m_jetCat = conf.getParameter<std::string>("JetWP");
}

bool FakeLeptonCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const
{
  const auto& electrons = producers.get<CandidatesProducer<pat::Electron>>("electrons").selected();
  const auto& muons = producers.get<CandidatesProducer<pat::Muon>>("muons").selected();
  if ( electrons.empty() && muons.empty() ) { // at least one lepton
    return false;
  }
  if ( m_maxMET > 0. ) {
    const auto& metProd = producers.get<METProducer>("met");
    if ( metProd.p4.Pt() > m_maxMET ) {
      return false;
    }
  }
  return true;
}

bool FakeLeptonCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const
{
  // exactly one loose lepton
  const TTWAnalyzer& ttW = analyzers.get<TTWAnalyzer>(m_ttwAnaName);
  const auto& selLeptons = ttW.selectedLeptons(m_selCat);
  if ( selLeptons.size() != 1 ) {
    return false;
  }
  LogDebug("Framework") << "Exactly one selected lepton, idx: " << selLeptons[0] << " (nLeptons: " << ttW.getObjList<TTWAnalysis::Lepton>().size() << ")";
  const auto& candLepton = ttW.getObjList<TTWAnalysis::Lepton>()[selLeptons[0]];
  if ( candLepton.isEl() != m_isEl ) {
    return false;
  }
  LogDebug("Framework") << "Correct flavour :-)";

  if ( ! m_jetCat.empty() ) { // a jet with DeltaR > 1.
    const auto& ttwJets = ttW.getPtrList<pat::Jet>();
    const auto& selJetIdx = ttW.selectedJets(m_jetCat);
    LogDebug("Framework") << "Searching for a jet with DeltaR > 1.; there are " << ttwJets.size() << " jets saved, " << selJetIdx.size() << " of which are " << m_jetCat;
    if ( ! std::any_of(std::begin(selJetIdx), std::end(selJetIdx),
              [&candLepton,&ttwJets] ( std::uint16_t i ) {
                return ROOT::Math::VectorUtil::DeltaR(ttwJets[i]->p4(), candLepton.p4()) > 1.;
              }) )
    {
      return false;
    }
  }

  if ( m_maxMT > 0. ) { // max. MT
    const auto& metProd = producers.get<METProducer>("met");
    const auto MT = std::sqrt( 2.*candLepton.p4().Pt()*metProd.p4.Pt()*
                         (1.-std::cos(candLepton.p4().Phi()-metProd.p4.Phi())) );
    if ( MT > m_maxMT ) {
      return false;
    }
  }
  return true;
}
