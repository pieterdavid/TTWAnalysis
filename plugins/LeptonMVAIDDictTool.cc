#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "cp3_llbb/Framework/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"

#include "Helpers.h"
#include "Math/VectorUtil.h"
#include "TMVA/Reader.h"

namespace TTWAnalysis {
/**
 * (prompt) lepton multivariate ID
 *
 * See
 * - https://twiki.cern.ch/twiki/bin/viewauth/CMS/LeptonMVA
 * - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTH-Multileptons
 */
class DictTTHElectronMVA76 : public DictTool<pat::Electron>, protected DictPVHelper {
public:
  DictTTHElectronMVA76(const edm::ParameterSet& config)
    : DictTool<pat::Electron>{config}
    , DictPVHelper{config}
    , m_vars{
        {"LepGood_pt"                   , -1.}
      , {"LepGood_eta"                  , -1.}
      , {"LepGood_jetNDauChargedMVASel" , -1.}
      , {"LepGood_miniRelIsoCharged"    , -1.}
      , {"LepGood_miniRelIsoNeutral"    , -1.}
      , {"LepGood_jetPtRelv2"           , -1.}
      , {"min(LepGood_jetPtRatiov2,1.5)", -1.}
      , {"max(LepGood_jetBTagCSV,0)"    , -1.}
      , {"LepGood_sip3d"                , -1.}
      , {"log(abs(LepGood_dxy))"        , -1.}
      , {"log(abs(LepGood_dz))"         , -1.}
      , {"LepGood_mvaIdSpring15"        , -1.}
    }
    , m_maxDR{config.getParameter<double>("JetLeptonDR")}
    , m_returnAllVars{config.getUntrackedParameter<bool>("AddAllVariablesToTree", false)}
  {
    for ( auto& entry : m_vars ) {
      m_tmvaReader.AddVariable(entry.first, &(entry.second));
    }
    edm::LogInfo("TTHLeptonMVA") << "Initialising TMVA reader for TTH electron MVA (76 version)";
    m_tmvaReader.BookMVA("ttH-el-BDTG", config.getParameter<edm::FileInPath>("WeightsFile").fullPath());
  }

  virtual ~DictTTHElectronMVA76() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictPVHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    m_jetsToken = collector.consumes<std::vector<pat::Jet>>(config.getParameter<edm::InputTag>("Jets"));
  }

  virtual Dict evaluate(edm::Ptr<pat::Electron> cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  mutable std::vector<std::pair<std::string,float>> m_vars;
  float& getVar(const std::string& name) const { auto it = std::find_if(std::begin(m_vars), std::end(m_vars), [&name] ( const std::pair<std::string,float>& item ) { return name == item.first; }); return it->second; }
  mutable TMVA::Reader m_tmvaReader;
  edm::EDGetTokenT<std::vector<pat::Jet>> m_jetsToken;
  double m_maxDR;
  bool m_returnAllVars;
};

class DictTTHElectronMVA80 : public DictTool<pat::Electron>, protected DictPVHelper {
public:
  DictTTHElectronMVA80(const edm::ParameterSet& config)
    : DictTool<pat::Electron>{config}
    , DictPVHelper{config}
    , m_vars{
        {"LepGood_pt"                   , -1.}
      , {"LepGood_eta"                  , -1.}
      , {"LepGood_jetNDauChargedMVASel" , -1.}
      , {"LepGood_miniRelIsoCharged"    , -1.}
      , {"LepGood_miniRelIsoNeutral"    , -1.}
      , {"LepGood_jetPtRelv2"           , -1.}
      , {"min(LepGood_jetPtRatiov2,1.5)", -1.}
      , {"max(LepGood_jetBTagCSV,0)"    , -1.}
      , {"LepGood_sip3d"                , -1.}
      , {"log(abs(LepGood_dxy))"        , -1.}
      , {"log(abs(LepGood_dz))"         , -1.}
      , {"LepGood_mvaIdSpring16HZZ"     , -1.}
    }
    , m_maxDR{config.getParameter<double>("JetLeptonDR")}
    , m_returnAllVars{config.getUntrackedParameter<bool>("AddAllVariablesToTree", false)}
  {
    for ( auto& entry : m_vars ) {
      m_tmvaReader.AddVariable(entry.first, &(entry.second));
    }
    edm::LogInfo("TTHLeptonMVA") << "Initialising TMVA reader for TTH electron MVA (80 version)";
    m_tmvaReader.BookMVA("ttH-el-BDTG", config.getParameter<edm::FileInPath>("WeightsFile").fullPath());
  }

  virtual ~DictTTHElectronMVA80() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictPVHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    m_jetsToken = collector.consumes<std::vector<pat::Jet>>(config.getParameter<edm::InputTag>("Jets"));
  }

  virtual Dict evaluate(edm::Ptr<pat::Electron> cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  mutable std::vector<std::pair<std::string,float>> m_vars;
  float& getVar(const std::string& name) const { auto it = std::find_if(std::begin(m_vars), std::end(m_vars), [&name] ( const std::pair<std::string,float>& item ) { return name == item.first; }); return it->second; }
  mutable TMVA::Reader m_tmvaReader;
  edm::EDGetTokenT<std::vector<pat::Jet>> m_jetsToken;
  double m_maxDR;
  bool m_returnAllVars;
};

class DictTTHMuonMVA : public DictTool<pat::Muon>, protected DictPVHelper {
public:
  DictTTHMuonMVA(const edm::ParameterSet& config)
    : DictTool<pat::Muon>{config}
    , DictPVHelper{config}
    , m_vars{
        {"LepGood_pt"                   , -1.}
      , {"LepGood_eta"                  , -1.}
      , {"LepGood_jetNDauChargedMVASel" , -1.}
      , {"LepGood_miniRelIsoCharged"    , -1.}
      , {"LepGood_miniRelIsoNeutral"    , -1.}
      , {"LepGood_jetPtRelv2"           , -1.}
      , {"min(LepGood_jetPtRatiov2,1.5)", -1.}
      , {"max(LepGood_jetBTagCSV,0)"    , -1.}
      , {"LepGood_sip3d"                , -1.}
      , {"log(abs(LepGood_dxy))"        , -1.}
      , {"log(abs(LepGood_dz))"         , -1.}
      , {"LepGood_segmentCompatibility" , -1.}
    }
    , m_maxDR{config.getParameter<double>("JetLeptonDR")}
    , m_returnAllVars{config.getUntrackedParameter<bool>("AddAllVariablesToTree", false)}
    , m_uniq{config.getParameter<std::string>("UniqueName")}
  {
    for ( auto& entry : m_vars ) {
      m_tmvaReader.AddVariable(entry.first, &(entry.second));
    }
    edm::LogInfo("TTHLeptonMVA") << "Initialising TMVA reader for TTH muon MVA";
    m_tmvaReader.BookMVA("ttH-mu-BDTG", config.getParameter<edm::FileInPath>("WeightsFile").fullPath());
  }

  virtual ~DictTTHMuonMVA() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override
  {
    DictPVHelper::doConsumes(config, std::forward<edm::ConsumesCollector>(collector));
    m_jetsToken = collector.consumes<std::vector<pat::Jet>>(config.getParameter<edm::InputTag>("Jets"));
  }

  virtual Dict evaluate(edm::Ptr<pat::Muon> cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  mutable std::vector<std::pair<std::string,float>> m_vars;
  float& getVar(const std::string& name) const { auto it = std::find_if(std::begin(m_vars), std::end(m_vars), [&name] ( const std::pair<std::string,float>& item ) { return name == item.first; }); return it->second; }
  mutable TMVA::Reader m_tmvaReader;
  edm::EDGetTokenT<std::vector<pat::Jet>> m_jetsToken;
  double m_maxDR;
  bool m_returnAllVars;
  std::string m_uniq;
};
}

// implementation
namespace {
  template<typename _Arg,typename _UnaryFunction>
  struct _LessOf {
    public:
      _LessOf(_UnaryFunction fun) : m_fun(fun) {}
      bool operator() ( const _Arg& __a, const _Arg& __b ) const
      {
        return m_fun(__a) < m_fun(__b);
      }
    private:
      _UnaryFunction m_fun;
  };
  // factory method
  template<typename _Arg,typename _UnaryFunction>
  _LessOf<_Arg,_UnaryFunction> LessOf( _UnaryFunction __fun )
  { return _LessOf<_Arg,_UnaryFunction>(std::forward<_UnaryFunction>(__fun)); }

  const pat::Jet* findMatchingJet( edm::Ptr<reco::RecoCandidate> lepton, const std::vector<pat::Jet>& jets, double maxDR )
  {
    using namespace ROOT::Math;
    const auto it = std::min_element(std::begin(jets), std::end(jets),
          LessOf<const pat::Jet&>( [lepton] ( const pat::Jet& j ) {
            return VectorUtil::DeltaR(j.p4(), lepton->p4());
          } ));
    if ( ( std::end(jets) != it ) && ( VectorUtil::DeltaR(lepton->p4(), it->p4()) < maxDR ) ) {
      return &*it;
    }
    return nullptr;
  }

  // Lepton-aware JEC
  reco::Candidate::LorentzVector jetLepAwareJEC(const pat::Jet* jet, const reco::Candidate::LorentzVector& lepton, bool print=false)
  {
    if ( ( ! jet ) || ( (jet->p4()*jet->jecFactor("Uncorrected")-lepton).Rho() < 1.e-4 ) ) {
      if ( print ) { LogDebug("TTHLeptonMVA") << "No matching jet with significantly different momentum -> returning lepton"; }
      return lepton;
    } else {
      const float rawCorr{jet->jecFactor("Uncorrected")};
      const float l1Corr{jet->jecFactor("L1FastJet")};
      if ( print ) {
        LogDebug("TTHLeptonMVA") << "Raw correction: " << rawCorr << "; L1 correction: " << l1Corr;
        LogTrace("TTHLeptonMVA") << "Available JEC set: ";
        for ( const auto& jecSetName : jet->availableJECSets() ) {
          LogTrace("TTHLeptonMVA") << jecSetName << ", ";
        }
        LogTrace("TTHLeptonMVA") << std::endl;
        LogTrace("TTHLeptonMVA") << "Available JEC levels: ";
        for ( const auto& jecLevelName : jet->availableJECLevels() ) {
          LogTrace("TTHLeptonMVA") << jecLevelName << " (" << jet->jecFactor(jecLevelName) << "), ";
        }
        LogTrace("TTHLeptonMVA") << std::endl;
      }
      return (jet->p4()-lepton*(1./l1Corr))+lepton;
    }
  }

  double ptRelv2(const pat::Jet* jet, const reco::Candidate::LorentzVector& lepton)
  {
    reco::Candidate::LorentzVector m{jetLepAwareJEC(jet, lepton)};
    return ( (m-lepton).Rho() < 1.e-4 ) ? 0. : lepton.Vect().Cross((m-lepton).Vect().Unit()).R();
  }

  int numberOfChargedDaughtersMVASel(const pat::Jet* jet, const reco::RecoCandidate* cand, const reco::Vertex* pv)
  {
    using namespace ROOT::Math;
    if ( ( ! jet ) || ( ! pv ) ) { return 0; }

    return std::count_if(std::begin(jet->daughterPtrVector()), std::end(jet->daughterPtrVector()),
        [jet,pv] ( reco::CandidatePtr daug ) {
          if ( ( VectorUtil::DeltaR(daug->p4(), jet->p4()) < 0.4 ) && ( daug->charge() != 0 ) ) {
            const pat::PackedCandidate* pkCand = dynamic_cast<const pat::PackedCandidate*>(daug.get());
            if ( pkCand->fromPV() >= pat::PackedCandidate::PVTight ) {
              auto trk = pkCand->pseudoTrack();
              return ( ( trk.pt() > 1. )
                    && ( trk.hitPattern().numberOfValidHits() >= 8 )
                    && ( trk.hitPattern().numberOfValidPixelHits() >= 2 )
                    && ( trk.normalizedChi2() < 5. )
                    && ( std::fabs(trk.dxy(pv->position())) < 0.2 )
                    && ( std::fabs(trk.dz(pv->position())) < 17. )
                     ) ;
            }
          }
          return false;
        });
  }

#ifdef EDM_ML_LOGDEBUG
  // only used for debug printout below
  std::string printMomentum( const reco::Candidate::LorentzVector& mom )
  {
    std::stringstream out;
    out << "(PT=" << mom.Pt() << ",ETA=" << mom.Eta() << ",PHI=" << mom.Phi() << ",E=" << mom.E() << ")";
    return out.str();
  }
#endif
}

TTWAnalysis::Dict TTWAnalysis::DictTTHElectronMVA76::evaluate(edm::Ptr<pat::Electron> cand,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  const bool valid{cand.isNonnull() && cand->superCluster().isNonnull()};
  const double eta = valid ? cand->superCluster()->eta() : 0.;
  const reco::Vertex* pv = event ? getPV(event) : nullptr;

  edm::Handle<std::vector<pat::Jet>> jets;
  const pat::Jet* jet{nullptr};
  if ( event && valid ) {
    event->getByToken(m_jetsToken, jets);
    jet = findMatchingJet( cand, *jets, m_maxDR );
  }
  int jetNDauChargedMVASel = ( valid && jet && pv ) ? numberOfChargedDaughtersMVASel(jet, cand.get(), pv) : 0;

#ifdef EDM_ML_LOGDEBUG
  LogDebug("TTHLeptonMVA") << "\n\nElectron prompt-lepton MVA for electron with momentum " << printMomentum(valid ? cand->p4() : reco::Candidate::LorentzVector{});
  LogDebug("TTHLeptonMVA") << "Best matching jet: " << ( jet ? printMomentum(jet->p4()) : "none" );
  if ( valid && jet ) {
    reco::Candidate::LorentzVector m{jetLepAwareJEC(jet, cand->p4(), true)};
    LogDebug("TTHLeptonMVA") << "After lepton-aware JEC: " << printMomentum(m);
  }
#endif

  getVar("LepGood_pt"                   ) = valid ? cand->pt() : -1.;
  getVar("LepGood_eta"                  ) = eta;
  getVar("LepGood_jetNDauChargedMVASel" ) = jetNDauChargedMVASel;
  getVar("LepGood_miniRelIsoCharged"    ) = valid ? cand->userFloat("miniIso_AbsCharged")/cand->pt() : -1.;
  getVar("LepGood_miniRelIsoNeutral"    ) = valid ? cand->userFloat("miniIso_AbsNeutral_rhoArea")/cand->pt() : -1.;
  getVar("LepGood_jetPtRelv2"           ) = valid ? ptRelv2(jet, cand->p4()) : -1.;
  getVar("min(LepGood_jetPtRatiov2,1.5)") = valid ? std::min(cand->pt()/jetLepAwareJEC(jet, cand->p4()).Pt(), 1.5) : -1.;
  getVar("max(LepGood_jetBTagCSV,0)"    ) = std::max(jet ? jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99., 0.);
  getVar("LepGood_sip3d"                ) = valid ? std::abs(cand->userFloat("dca")) : -1;
  getVar("log(abs(LepGood_dxy))"        ) = valid ? std::log(std::abs(cand->userFloat("dxy"))) : 0.;
  getVar("log(abs(LepGood_dz))"         ) = valid ? std::log(std::abs(cand->userFloat("dz" ))) : 0.;
  getVar("LepGood_mvaIdSpring15"        ) = valid ? cand->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") : -1.;

  TTWAnalysis::Dict ret{};
  ret.add("LeptonMVA76"   , valid ? m_tmvaReader.EvaluateMVA("ttH-el-BDTG") : -1.);
  ret.add("MVAIdNonTrig25nsSpring15", getVar("LepGood_mvaIdSpring15"));
  if ( m_returnAllVars ) { // some more variables that are only calculated here
    ret.add("jetPtRelv2"    , getVar("LepGood_jetPtRelv2"));
    ret.add("jetPtRatio"    , getVar("min(LepGood_jetPtRatiov2,1.5)"));
    ret.add("jetNDaugMVASel", getVar("LepGood_jetNDauChargedMVASel"));
    ret.add("jetBTagCSV"    , getVar("max(LepGood_jetBTagCSV,0)"));
    //
    ret.add("miniRelIsoCharged", getVar("LepGood_miniRelIsoCharged"));
    ret.add("miniRelIsoNeutral", getVar("LepGood_miniRelIsoNeutral"));
  }

  return ret;
}

TTWAnalysis::Dict TTWAnalysis::DictTTHElectronMVA80::evaluate(edm::Ptr<pat::Electron> cand,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  const bool valid{cand.isNonnull() && cand->superCluster().isNonnull()};
  const double eta = valid ? cand->superCluster()->eta() : 0.;
  const reco::Vertex* pv = event ? getPV(event) : nullptr;

  edm::Handle<std::vector<pat::Jet>> jets;
  const pat::Jet* jet{nullptr};
  if ( event && valid ) {
    event->getByToken(m_jetsToken, jets);
    jet = findMatchingJet( cand, *jets, m_maxDR );
  }
  int jetNDauChargedMVASel = ( valid && jet && pv ) ? numberOfChargedDaughtersMVASel(jet, cand.get(), pv) : 0;

#ifdef EDM_ML_LOGDEBUG
  LogDebug("TTHLeptonMVA") << "\n\nElectron prompt-lepton MVA for electron with momentum " << printMomentum(valid ? cand->p4() : reco::Candidate::LorentzVector{});
  LogDebug("TTHLeptonMVA") << "Best matching jet: " << ( jet ? printMomentum(jet->p4()) : "none" );
  if ( valid && jet ) {
    reco::Candidate::LorentzVector m{jetLepAwareJEC(jet, cand->p4(), true)};
    LogDebug("TTHLeptonMVA") << "After lepton-aware JEC: " << printMomentum(m);
  }
#endif

  getVar("LepGood_pt"                   ) = valid ? cand->pt() : -1.;
  getVar("LepGood_eta"                  ) = eta;
  getVar("LepGood_jetNDauChargedMVASel" ) = jetNDauChargedMVASel;
  getVar("LepGood_miniRelIsoCharged"    ) = valid ? cand->userFloat("miniIso_AbsCharged")/cand->pt() : -1.;
  getVar("LepGood_miniRelIsoNeutral"    ) = valid ? cand->userFloat("miniIso_AbsNeutral_rhoArea")/cand->pt() : -1.;
  getVar("LepGood_jetPtRelv2"           ) = valid ? ptRelv2(jet, cand->p4()) : -1.;
  getVar("min(LepGood_jetPtRatiov2,1.5)") = valid ? std::min(cand->pt()/jetLepAwareJEC(jet, cand->p4()).Pt(), 1.5) : -1.;
  getVar("max(LepGood_jetBTagCSV,0)"    ) = std::max(jet ? jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99., 0.);
  getVar("LepGood_sip3d"                ) = valid ? std::abs(cand->userFloat("dca")) : -1;
  getVar("log(abs(LepGood_dxy))"        ) = valid ? std::log(std::abs(cand->userFloat("dxy"))) : 0.;
  getVar("log(abs(LepGood_dz))"         ) = valid ? std::log(std::abs(cand->userFloat("dz" ))) : 0.;
  getVar("LepGood_mvaIdSpring16HZZ"     ) = valid ? cand->userFloat("HZZSpring16MVAID_value") : -1.; // cand->userFloat("ElectronMVAEstimatorRun2Spring16HZZV1Values")

  TTWAnalysis::Dict ret{};
  ret.add("LeptonMVA80"     , valid ? m_tmvaReader.EvaluateMVA("ttH-el-BDTG") : -1.);
  ret.add("MVAIdSpring16HZZ", getVar("LepGood_mvaIdSpring16HZZ"));
  if ( m_returnAllVars ) { // some more variables that are only calculated here
    ret.add("jetPtRelv2"    , getVar("LepGood_jetPtRelv2"));
    ret.add("jetPtRatio"    , getVar("min(LepGood_jetPtRatiov2,1.5)"));
    ret.add("jetNDaugMVASel", getVar("LepGood_jetNDauChargedMVASel"));
    ret.add("jetBTagCSV"    , getVar("max(LepGood_jetBTagCSV,0)"));
    //
    ret.add("miniRelIsoCharged", getVar("LepGood_miniRelIsoCharged"));
    ret.add("miniRelIsoNeutral", getVar("LepGood_miniRelIsoNeutral"));
  }

  return ret;
}

TTWAnalysis::Dict TTWAnalysis::DictTTHMuonMVA::evaluate(edm::Ptr<pat::Muon> cand,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  const bool valid{cand.isNonnull()};
  const reco::Vertex* pv = event ? getPV(event) : nullptr;

  edm::Handle<std::vector<pat::Jet>> jets;
  const pat::Jet* jet{nullptr};
  if ( event && valid ) {
    event->getByToken(m_jetsToken, jets);
    jet = findMatchingJet( cand, *jets, m_maxDR );
  }
  int jetNDauChargedMVASel = ( valid && jet && pv ) ? numberOfChargedDaughtersMVASel(jet, cand.get(), pv) : 0;

#ifdef EDM_ML_LOGDEBUG
  LogDebug("TTHLeptonMVA") << "\n\nMuon prompt-lepton MVA for muon with momentum " << printMomentum(valid ? cand->p4() : reco::Candidate::LorentzVector{});
  LogDebug("TTHLeptonMVA") << "Best matching jet: " << ( jet ? printMomentum(jet->p4()) : "none" );
  if ( valid && jet ) {
    reco::Candidate::LorentzVector m{jetLepAwareJEC(jet, cand->p4(), true)};
    LogDebug("TTHLeptonMVA") << "After lepton-aware JEC: " << printMomentum(m);
  }
#endif

  getVar("LepGood_pt"                   ) = valid ? cand->pt() : -1.;
  getVar("LepGood_eta"                  ) = valid ? cand->eta() : -100.;
  getVar("LepGood_jetNDauChargedMVASel" ) = jetNDauChargedMVASel;
  getVar("LepGood_miniRelIsoCharged"    ) = valid ? cand->userFloat("miniIso_AbsCharged")/cand->pt() : -1.;
  getVar("LepGood_miniRelIsoNeutral"    ) = valid ? cand->userFloat("miniIso_AbsNeutral_rhoArea")/cand->pt() : -1.;
  getVar("LepGood_jetPtRelv2"           ) = valid ? ptRelv2(jet, cand->p4()): -1.;
  getVar("min(LepGood_jetPtRatiov2,1.5)") = valid ? std::min(cand->pt()/jetLepAwareJEC(jet, cand->p4()).Pt(), 1.5) : -1.;
  getVar("max(LepGood_jetBTagCSV,0)"    ) = std::max(jet ? jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99., 0.);
  getVar("LepGood_sip3d"                ) = valid ? std::abs(cand->dB(pat::Muon::PV3D)/cand->edB(pat::Muon::PV3D)) : -1.;
  getVar("log(abs(LepGood_dxy))"        ) = pv && valid && cand->muonBestTrack().isNonnull() ? std::log(std::abs(cand->muonBestTrack()->dxy(pv->position()))) : 0.; // TODO check definition
  getVar("log(abs(LepGood_dz))"         ) = pv && valid && cand->muonBestTrack().isNonnull() ? std::log(std::abs(cand->muonBestTrack()->dz (pv->position()))) : 0.;
  getVar("LepGood_segmentCompatibility" ) = valid ? cand->segmentCompatibility() : false;

  TTWAnalysis::Dict ret{};
  ret.add("LeptonMVA"+m_uniq, valid ? m_tmvaReader.EvaluateMVA("ttH-mu-BDTG") : -1.);
  // some more variables that are only calculated here
  if ( m_returnAllVars ) { // some more variables that are only calculated here
    ret.add("jetPtRelv2"    , getVar("LepGood_jetPtRelv2"));
    ret.add("jetPtRatio"    , getVar("min(LepGood_jetPtRatiov2,1.5)"));
    ret.add("jetNDaugMVASel", getVar("LepGood_jetNDauChargedMVASel"));
    ret.add("jetBTagCSV"    , getVar("max(LepGood_jetBTagCSV,0)"));
    //
    ret.add("miniRelIsoCharged", getVar("LepGood_miniRelIsoCharged"));
    ret.add("miniRelIsoNeutral", getVar("LepGood_miniRelIsoNeutral"));
    ret.add("segmentCompatibility", getVar("LepGood_segmentCompatibility" ));
  }

  return ret;
}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictTTHElectronMVA76, "ttw_electronMVAttH76");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Electron>::factory, TTWAnalysis::DictTTHElectronMVA80, "ttw_electronMVAttH80");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<pat::Muon    >::factory, TTWAnalysis::DictTTHMuonMVA, "ttw_muonMVAttH");
