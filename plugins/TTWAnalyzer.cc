//#include <cp3_llbb/TTWAnalysis/interface/Defines.h>
#include "cp3_llbb/TTWAnalysis/interface/TTWAnalyzer.h"
#include "cp3_llbb/TTWAnalysis/interface/stl_helpers.h"

#include <cp3_llbb/Framework/interface/METProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <Math/PtEtaPhiE4D.h>
// #include <Math/LorentzVector.h>
#include <Math/VectorUtil.h>

// To access VectorUtil::DeltaR() more easily
using namespace ROOT::Math;

using namespace TTWAnalysis;

using ElectronHybridCut = StringCutObjectSelector<pat::Electron>;
using MuonHybridCut     = StringCutObjectSelector<pat::Muon>;
using JetHybridCut      = StringCutObjectSelector<pat::Jet>;

namespace {
  template<typename _Value,typename _InputIterator,typename _UnaryFunction>
  _Value min_value( _InputIterator __first, _InputIterator __last, _UnaryFunction&& __fun )
  {
    using _ElmArg = typename std::iterator_traits<_InputIterator>::reference;
    return std::accumulate(__first, __last, std::numeric_limits<_Value>::max(),
        [&__fun] (_Value __iv, _ElmArg __elm) -> _Value {
          _Value __ov = std::forward<_UnaryFunction>(__fun)(__elm);
          return ( __iv > __ov ) ? __ov : __iv;
        });
  }

  float minDRjl( const std::vector<Lepton>& leptons, const indexlist_t& lIdxs, const pat::Jet& j )
  {
    return min_value<float>(std::begin(lIdxs), std::end(lIdxs),
        [j,leptons] ( index_t iL ) { return VectorUtil::DeltaR(j.p4(), leptons[iL].p4()); });
  }

  float minDRjl( const std::vector<Lepton>& leptons, const indexlist_t& lIdxs, const DiJet& jj )
  {
    return std::min(minDRjl(leptons, lIdxs, *(jj.first)), minDRjl(leptons, lIdxs, *(jj.second)));
  }

  LeptonCut makeLeptonCut( const edm::ParameterSet& cutPerFlavour )
  {
    return LeptonCut{
        ElectronHybridCut{cutPerFlavour.getParameter<std::string>("Electron")}
      , MuonHybridCut    {cutPerFlavour.getParameter<std::string>("Muon"    )}
    };
  }
  DiLeptonCut makeDiLeptonCut( const edm::ParameterSet& cutPerLepton )
  {
    return DiLeptonCut{
        makeLeptonCut(cutPerLepton.getParameter<edm::ParameterSet>("Leading"))
      , makeLeptonCut(cutPerLepton.getParameter<edm::ParameterSet>("SubLeading"))
    };
  }
  std::pair<std::string,JetCut> makeJetCut( const edm::ParameterSet& config )
  {
    return std::make_pair(config.getParameter<std::string>("LeptonWP")
        , JetHybridCut{config.getParameter<std::string>("Selection")});
  }
  std::pair<std::string,DiJetCut> makeDiJetCut( const edm::ParameterSet& config )
  {
    return std::make_pair(config.getParameter<std::string>("LeptonWP"), DiJetCut{
        JetHybridCut{config.getParameter<std::string>("Leading")}
      , JetHybridCut{config.getParameter<std::string>("SubLeading")}});
  }
}

float TTWAnalysis::DeltaEta(const myLorentzVector& v1, const myLorentzVector& v2) {
  return std::abs(v1.Eta() - v2.Eta());
}

TTWAnalyzer::TTWAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config)
  : Analyzer(name, tree_, config),

  // Not untracked as these parameters are mandatory
  m_electrons_producer(config.getParameter<std::string>("electronsProducer")),
  m_muons_producer(config.getParameter<std::string>("muonsProducer")),
  m_jets_producer(config.getParameter<std::string>("jetsProducer")),
  m_met_producer(config.getParameter<std::string>("metProducer")),

  m_elWP{}, m_muWP{}, m_lWP{}, m_llWP{}, m_bWP{}, m_bbWP{}, m_lljjWP{}, m_llbbWP{}, m_llbbmWP{},
  m_idxElWP{}, m_idxMuWP{}, m_idxlWP{}, m_idxllWP{}, m_idxjDRWP{}, m_idxbDRWP_PT{}, m_idxbDRWP_tag{}, m_idxjjDRWP{}, m_idxbbDRWP_PT{}, m_idxbbDRWP_tag{}, m_idxlljjDRWP{}, m_idxllbbDRWP_PT{}, m_idxllbbDRWP_tag{}, m_idxlljjmDRWP{}, m_idxllbbmDRWP_PT{}, m_idxllbbmDRWP_tag{},

  m_electronPtCut( config.getUntrackedParameter<double>("electronPtCut", 20) ),
  m_electronEtaCut( config.getUntrackedParameter<double>("electronEtaCut", 2.5) ),

  m_muonPtCut( config.getUntrackedParameter<double>("muonPtCut", 20) ),
  m_muonEtaCut( config.getUntrackedParameter<double>("muonEtaCut", 2.4) ),

  m_jetPtCut( config.getUntrackedParameter<double>("jetPtCut", 30) ),
  m_jetEtaCut( config.getUntrackedParameter<double>("jetEtaCut", 2.5) ),
  m_jetPUID( config.getUntrackedParameter<double>("jetPUID", std::numeric_limits<float>::min()) ),
  m_jetDRleptonCut( config.getUntrackedParameter<double>("jetDRleptonCut", 0.3) ),

  m_jetBTagName( config.getParameter<std::string>("bTagName") ),

  m_hltDRCut( config.getUntrackedParameter<double>("hltDRCut", std::numeric_limits<float>::max()) ),
  m_hltDPtCut( config.getUntrackedParameter<double>("hltDPtCut", std::numeric_limits<float>::max()) ),

  m_maxNumLForLL(config.getUntrackedParameter<uint32_t>("MaxNumLForLL",  5)),
  m_maxNumJForJJ(config.getUntrackedParameter<uint32_t>("MaxNumJForJJ", 10)),
  m_maxNumLLForLLJJ(config.getUntrackedParameter<uint32_t>("MaxNumLLForLLJJ", 10)),
  m_maxNumJJForLLJJ(config.getUntrackedParameter<uint32_t>("MaxNumJJForLLJJ", 10))
{
  const auto& elWPCfg = config.getParameter<edm::ParameterSet>("ElectronWP");
  for ( const auto elWPName : elWPCfg.getParameterNames() ) {
    m_elWP[elWPName] = ElectronHybridCut{config.getParameter<std::string>(elWPName)};
  }
  m_idxElWP = IdxListForWPsHolder<decltype(m_elWP)>(m_elWP, tree_, "electrons_");
  const auto& muWPCfg = config.getParameter<edm::ParameterSet>("MuonWP");
  for ( const auto muWPName : muWPCfg.getParameterNames() ) {
    m_muWP[muWPName] = MuonHybridCut{config.getParameter<std::string>(muWPName)};
  }
  m_idxMuWP = IdxListForWPsHolder<decltype(m_muWP)>(m_muWP, tree_, "muons_");

  const auto& leptonWPCfg = config.getParameter<edm::ParameterSet>("LeptonWP");
  for ( const auto leptonWPName : leptonWPCfg.getParameterNames() ) {
    m_lWP[leptonWPName] = makeLeptonCut(leptonWPCfg.getParameter<edm::ParameterSet>(leptonWPName));
  }
  m_idxlWP = IdxListForWPsHolder<decltype(m_lWP)>(m_lWP, tree, "leptons_");

  const auto& llWPCfg = config.getParameter<edm::ParameterSet>("DiLeptonWP");
  for ( const auto llWPName : llWPCfg.getParameterNames() ) {
    m_llWP[llWPName] = makeDiLeptonCut(llWPCfg.getParameter<edm::ParameterSet>(llWPName));
  }
  m_idxllWP = IdxListForWPsHolder<decltype(m_llWP)>(m_llWP, tree, "dileptons_");

  m_jetIDCut = JetHybridCut{config.getParameter<std::string>("JetID")};

  m_idxjDRWP = IdxListForWPsHolder<decltype(m_lWP)>(m_lWP, tree, "selJets_selID_DRCut_");
  m_idxjjDRWP = IdxListForWPsHolder<decltype(m_lWP)>(m_lWP, tree, "diJets_DRCut_");

  const auto& bWPCfg = config.getParameter<edm::ParameterSet>("BtagWP");
  for ( const auto & bWPName : bWPCfg.getParameterNames() ) {
    m_bWP[bWPName] = makeJetCut(bWPCfg.getParameter<edm::ParameterSet>(bWPName));
  }
  m_idxbDRWP_PT  = IdxListForWPsHolder<decltype(m_bWP)>(m_bWP, tree, "selBJets_DRCut_BWP_PtOrdered_");
  m_idxbDRWP_tag = IdxListForWPsHolder<decltype(m_bWP)>(m_bWP, tree, "selBJets_DRCut_BWP_tagOrdered_");

  const auto& bbWPCfg = config.getParameter<edm::ParameterSet>("DiBtagWP");
  for ( const auto& bbWPName : bbWPCfg.getParameterNames() ) {
    m_bbWP[bbWPName] = makeDiJetCut(bbWPCfg.getParameter<edm::ParameterSet>(bbWPName));
  }
  m_idxbbDRWP_PT  = IdxListForWPsHolder<decltype(m_bbWP)>(m_bbWP, tree, "diBJets_DRCut_BWP_PtOrdered_");
  m_idxbbDRWP_tag = IdxListForWPsHolder<decltype(m_bbWP)>(m_bbWP, tree, "diBJets_DRCut_BWP_tagOrdered_");
}


void TTWAnalyzer::doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector)
{
  m_vertices_token = collector.consumes<std::vector<reco::Vertex>>(config.getUntrackedParameter<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices")));
  const edm::ParameterSet& elVidCfg = config.getParameter<edm::ParameterSet>("electronVIDs");
  for ( const std::string& name : elVidCfg.getParameterNames() ) {
    m_el_vidTokens.emplace_back(std::make_pair(name, collector.consumes<edm::ValueMap<bool>>(elVidCfg.getParameter<edm::InputTag>(name))));
  }
}

void TTWAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers, const AnalyzersManager& analyzers, const CategoryManager& categories)
{
  using namespace TTWAnalysis;

  m_electrons.clear(); m_muons.clear(); m_leptons.clear(); m_dileptons.clear(); m_jets.clear(); m_dijets.clear(); m_dileptondijets.clear(); m_dileptondijetmets.clear();

  LogDebug("ttW") << "Begin event.";

  // get the primary vertex
  event.getByToken(m_vertices_token, m_vertices_handle);
  const reco::Vertex& pv = (*m_vertices_handle)[0];

  ///////////////////////////
  //        LEPTONS        //
  ///////////////////////////

  m_leptons.clear();

  LogDebug("ttW") << "Electrons";

  //
  std::vector<edm::Handle<edm::ValueMap<bool>>> el_vids;
  std::for_each(std::begin(m_el_vidTokens), std::end(m_el_vidTokens),
      [&event,&el_vids] ( const auto& nmAndTok ) {
        edm::Handle<edm::ValueMap<bool>> handle;
        event.getByToken(nmAndTok.second, handle);
        el_vids.emplace_back(std::move(handle));
      });
  //

  using ElectronsProducer = TTWAnalysis::CandidatesProducer<pat::Electron>;
  m_electrons = producers.get<ElectronsProducer>(m_electrons_producer).selected();
  for ( std::size_t i{0}; m_electrons.size() != i; ++i ) {
    auto iEl = m_electrons[i];
    if ( ( iEl->pt() > m_electronPtCut ) && ( std::abs(iEl->eta()) < m_electronEtaCut ) ) {
      m_leptons.emplace_back(i, iEl);
    }
    { // TODO consider putting this elsewhere
      pat::Electron* elNonConst = const_cast<pat::Electron*>(iEl.get());
      for ( std::size_t i{0}; i != el_vids.size(); ++i ) {
        try {
          elNonConst->addUserInt(m_el_vidTokens[i].first, ( (*(el_vids[i]))[iEl] ) ? 1 : 0);
        } catch ( const edm::Exception& e ) {
          std::stringstream msg;
          msg << "Problem getting id '" << m_el_vidTokens[i].first << "' for electron " << iEl.key() << " from " << iEl.id();
          throw edm::Exception(edm::errors::InvalidReference, msg.str());
        }
        LogDebug("ttW-electronID") << "Electron ID " << m_el_vidTokens[i].first << " for electron " << iEl.key() << " from " << iEl.id() << " : " << (*(el_vids[i]))[iEl] << " (" << el_vids[i]->size() << " in map)";
      }
    }
  }
  for ( decltype(m_idxElWP)::iterator ielWP{m_idxElWP.begin()}; m_idxElWP.end() != ielWP; ++ielWP ) {
    for ( std::size_t i{0}; m_electrons.size() != i; ++i ) {
      if ( ielWP.cut()(*(m_electrons[i])) ) {
        ielWP.idxList()->push_back(i);
      }
    }
  }

  LogDebug("ttW") << "Muons";

  using MuonsProducer = TTWAnalysis::CandidatesProducer<pat::Muon>;
  m_muons = producers.get<MuonsProducer>(m_muons_producer).selected();
  for ( std::size_t i{0}; m_muons.size() != i; ++i ) {
    const auto iMu = m_muons[i];
    if ( ( iMu->pt() > m_muonPtCut ) && ( std::abs(iMu->eta()) < m_muonEtaCut ) ) {
      m_leptons.emplace_back(i, iMu);
    }
    { // TODO consider putting this elsewhere
      pat::Muon* muNonConst = const_cast<pat::Muon*>(iMu.get());
      muNonConst->addUserInt("tightMuonID", iMu->isTightMuon(pv));
    }
  }
  for ( decltype(m_idxMuWP)::iterator imuWP{m_idxMuWP.begin()}; m_idxMuWP.end() != imuWP; ++imuWP ) {
    for ( std::size_t i{0}; m_muons.size() != i; ++i ) {
      if ( imuWP.cut()(*(m_muons[i])) ) {
        imuWP.idxList()->push_back(i);
      }
    }
  }

  std::sort(m_leptons.begin(), m_leptons.end(), MoreOf<const Lepton&>([] ( const Lepton& l ) { return l.recoCand()->pt(); }));

  for ( decltype(m_idxlWP)::iterator ilWP{m_idxlWP.begin()}; m_idxlWP.end() != ilWP; ++ilWP ) {
    for ( std::size_t iL{0}; m_leptons.size() != iL; ++iL ) {
      if ( ilWP.cut()(m_leptons[iL]) ) {
        ilWP.idxList()->push_back(iL);
      }
    }
  }

  ///////////////////////////
  //       TRIGGER         //
  ///////////////////////////
  LogDebug("ttW") << "Trigger";
  if (producers.exists("hlt")) {
    const HLTProducer& hlt = producers.get<HLTProducer>("hlt");
    if (hlt.paths.empty()) {
      LogDebug("ttW") << "No HLT path triggered for this event. Skipping HLT matching.";
    } else {
      LogDebug("ttW") << "HLT path triggered for this event:";
      for (const std::string& path: hlt.paths) {
        LogDebug("ttW") << "\t" << path;
      }
      for ( auto& lepton : m_leptons ) {
        lepton.hlt_idx = matchOfflineLepton(lepton, hlt);
      }
    }
  }

  ///////////////////////////
  //       DILEPTONS       //
  ///////////////////////////

  if ( m_maxNumLForLL ) {
    LogDebug("ttW") << "Dileptons";

    const std::size_t nLForLL = std::min(m_maxNumLForLL, m_leptons.size());
    if ( nLForLL != m_leptons.size() ) {
      edm::LogWarning("ttW") << "Limiting the number of leptons for dileptons to " << nLForLL;
    }
    for(uint16_t i1 = 0; i1 < nLForLL; i1++) {
      for(uint16_t i2 = i1 + 1; i2 < nLForLL; i2++) {
        const Lepton& l1 = m_leptons[i1];
        const Lepton& l2 = m_leptons[i2];

        m_dileptons.emplace_back(i1,l1,i2,l2);
      }
    }

    for ( decltype(m_idxllWP)::iterator illWP{m_idxllWP.begin()}; m_idxllWP.end() != illWP; ++illWP ) {
      for ( std::size_t iLL{0}; m_dileptons.size() != iLL; ++iLL ) {
        if ( illWP.cut()(m_dileptons[iLL]) ) {
          illWP.idxList()->push_back(iLL);
        }
      }
    }
  }

  ///////////////////////////
  //       JETS            //
  ///////////////////////////

  LogDebug("ttW") << "Jets";

  using JetsProducer = TTWAnalysis::CandidatesProducer<pat::Jet>;
  const auto& jets = producers.get<JetsProducer>(m_jets_producer).selected();
  // std::cout << "Retrieved " << jets.size() << " jets" << std::endl;
  for ( std::size_t i{0}; jets.size() != i; ++i ) {
    if ( ( jets[i]->pt() > m_jetPtCut ) && ( std::abs(jets[i]->eta()) < m_jetEtaCut ) ) {
      m_jets.push_back(jets[i]);
    }
  }
  indexlist_t selJets_selID;
  for ( std::size_t iJ{0}; m_jets.size() != iJ; ++iJ ) {
    if ( m_jetIDCut(*(m_jets[iJ])) ) {
      selJets_selID.push_back(iJ);
    }
  }
  // std::cout << "ttW: Selected " << selJets_selID.size() << " jets (after jet ID)" << std::endl;

  for ( decltype(m_idxjDRWP)::iterator iWP{m_idxjDRWP.begin()}; m_idxjDRWP.end() != iWP; ++iWP ) {
    for ( auto iJ : selJets_selID ) {
      if ( minDRjl(m_leptons, *(m_idxlWP.at(iWP.index()).idxList()), m_jets[iJ]) > m_jetDRleptonCut ) {
        iWP.idxList()->push_back(iJ);
      }
    }
  }
  for ( decltype(m_idxbDRWP_PT)::iterator iWP{m_idxbDRWP_PT.begin()}; m_idxbDRWP_PT.end() != iWP; ++iWP ) {
    for ( auto iJ : *(m_idxjDRWP.find(iWP.cut().first).idxList()) ) {
      if ( iWP.cut().second(m_jets[iJ]) ) {
        iWP.idxList()->push_back(iJ);
      }
    }
  }

  for ( std::size_t i{0}; m_idxbDRWP_PT.size() != i; ++i ) {
    indexlist_t& tagIdx = *(m_idxbDRWP_tag.at(i).idxList());
    tagIdx = *(m_idxbDRWP_PT.at(i).idxList());
    std::stable_sort(std::begin(tagIdx), std::end(tagIdx),
        MoreOf<index_t>([this] ( index_t iJ ) { return m_jets[iJ]->bDiscriminator(m_jetBTagName); }) );
  }

  ///////////////////////////
  //       DIJETS          //
  ///////////////////////////

  if ( m_maxNumJForJJ ) {
    LogDebug("ttW") << "Dijets";

    // Next, construct DiJets out of selected jets with selected ID (not accounting for minDRjl here)
    const std::size_t nJForJJ = std::min(m_maxNumJForJJ, selJets_selID.size());
    if ( nJForJJ != selJets_selID.size() ) {
      edm::LogWarning("ttW") << "Limiting the number of jets for dijets to " << nJForJJ;
    }
    for ( index_t i1{0}; nJForJJ != i1; ++i1 ) {
      index_t iJ1 = selJets_selID[i1];
      for ( index_t i2{static_cast<index_t>(i1+1)}; nJForJJ != i2; ++i2 ) {
        index_t iJ2 = selJets_selID[i2];
        m_dijets.emplace_back(iJ1, m_jets[iJ1], iJ2, m_jets[iJ2]);
      }
    }
    for ( decltype(m_idxjjDRWP)::iterator iWP{m_idxjjDRWP.begin()}; m_idxjjDRWP.end() != iWP; ++iWP ) {
      for ( index_t iJJ{0}; m_dijets.size() != iJJ; ++iJJ ) {
        if ( minDRjl(m_leptons, *(m_idxlWP.at(iWP.index()).idxList()), m_dijets[iJJ]) > m_jetDRleptonCut ) {
          iWP.idxList()->push_back(iJJ);
        }
      }
    }
    for ( decltype(m_idxbbDRWP_PT)::iterator iWP{m_idxbbDRWP_PT.begin()}; m_idxbbDRWP_PT.end() != iWP; ++iWP ) {
      for ( auto iJJ : *(m_idxjjDRWP.find(iWP.cut().first).idxList()) ) {
        if ( iWP.cut().second(m_dijets[iJJ]) ) {
          iWP.idxList()->push_back(iJJ); // NOTE di-b working points may include an ETA cut
        }
      }
    }

    for ( std::size_t i{0}; m_idxbbDRWP_PT.size() != i; ++i ) {
      indexlist_t& tagIdx = *(m_idxbbDRWP_tag.at(i).idxList());
      tagIdx = *(m_idxbbDRWP_PT.at(i).idxList());
      std::stable_sort(std::begin(tagIdx), std::end(tagIdx),
          MoreOf<index_t>([this] ( index_t iJJ ) {
            return m_dijets[iJJ].first->bDiscriminator(m_jetBTagName)
                 + m_dijets[iJJ].second->bDiscriminator(m_jetBTagName);
          }));
    }
  }

  ///////////////////////////
  //    EVENT VARIABLES    //
  ///////////////////////////

  if ( m_maxNumLLForLLJJ && m_maxNumJJForLLJJ ) {
    LogDebug("ttW") << "Dileptons-dijets";

    const std::size_t nLLForLLJJ = std::min(m_maxNumLLForLLJJ, m_dileptons.size());
    const std::size_t nJJForLLJJ = std::min(m_maxNumJJForLLJJ, m_dijets.size());
    if ( nLLForLLJJ != m_dileptons.size() ) {
      edm::LogWarning("ttW") << "Limited the number of dileptons for lljj combinations to " << nLLForLLJJ;
    }
    if ( nJJForLLJJ != m_dijets.size() ) {
      edm::LogWarning("ttW") << "Limited the number of dijets for lljj combinations to " << nJJForLLJJ;
    }
    // leptons-(b-)jets
    for ( index_t iLL{0}; nLLForLLJJ != iLL; ++iLL ) {
      for ( index_t iJJ{0}; nJJForLLJJ != iJJ; ++iJJ ) {
        m_dileptondijets.emplace_back(iLL, m_dileptons[iLL], iJJ, m_dijets[iJJ]);
        // TODO add delta R, delta Eta and DeltaPhi combinatorics (between jet and lepton)
      }
    }
    for ( decltype(m_idxlljjDRWP)::iterator iWP{m_idxlljjDRWP.begin()}; m_idxlljjDRWP.end() != iWP; ++iWP ) {
      for ( index_t iLLJJ{0}; m_dileptondijets.size() != iLLJJ; ++iLLJJ ) {
        if ( ( minDRjl(m_leptons, *(m_idxlWP.find(iWP.cut().first).idxList()), *(m_dileptondijets[iLLJJ].jj)) > m_jetDRleptonCut ) && ( iWP.cut().second(*(m_dileptondijets[iLLJJ].ll)) ) ) {
          iWP.idxList()->push_back(iLLJJ);
        }
      }
    }
    for ( decltype(m_idxllbbDRWP_PT)::iterator iWP{m_idxllbbDRWP_PT.begin()}; m_idxllbbDRWP_PT.end() != iWP; ++iWP ) {
      for ( auto iLLJJ : *(m_idxlljjDRWP.find(iWP.cut().first).idxList()) ) {
        if ( iWP.cut().second(*(m_dileptondijets[iLLJJ].jj)) ) {
          iWP.idxList()->push_back(iLLJJ); // NOTE di-b working points may include an ETA cut
        }
      }
    }

    // Order selected di-lepton-di-b-jets according to decreasing CSVv2 discriminant
    for ( std::size_t i{0}; m_idxllbbDRWP_PT.size() != i; ++i ) {
      indexlist_t& tagIdx = *(m_idxllbbDRWP_tag.at(i).idxList());
      tagIdx = *(m_idxllbbDRWP_PT.at(i).idxList());
      std::stable_sort(std::begin(tagIdx), std::end(tagIdx),
          MoreOf<index_t>([this] ( index_t iLLJJ ) {
            return m_dileptondijets[iLLJJ].jj->first->bDiscriminator(m_jetBTagName)
                 + m_dileptondijets[iLLJJ].jj->second->bDiscriminator(m_jetBTagName);
          }));
    }

    // leptons-(b-)jets-MET

    LogDebug("ttW") << "Dileptons-Dijets-MET";

    // Using regular MET
    const myLorentzVector met{producers.get<METProducer>(m_met_producer).p4};
    for ( index_t illjj{0}; m_dileptondijets.size() != illjj; ++illjj ) {
      m_dileptondijetmets.emplace_back(illjj, m_dileptondijets[illjj], met);
      // TODO delta(R,eta,phi) l-met, j-met, min&max
    }
    for ( decltype(m_idxlljjmDRWP)::iterator iWP{m_idxlljjmDRWP.begin()}; m_idxlljjmDRWP.end() != iWP; ++iWP ) {
      for ( index_t iLLJJM{0}; m_dileptondijetmets.size() != iLLJJM; ++iLLJJM ) {
        if ( ( minDRjl(m_leptons, *(m_idxlWP.find(iWP.cut().first).idxList()), *(m_dileptondijetmets[iLLJJM].lljj->jj)) > m_jetDRleptonCut ) && ( iWP.cut().second(*(m_dileptondijetmets[iLLJJM].lljj->ll)) ) ) {
          iWP.idxList()->push_back(iLLJJM);
        }
      }
    }
    for ( decltype(m_idxllbbmDRWP_PT)::iterator iWP{m_idxllbbmDRWP_PT.begin()}; m_idxllbbmDRWP_PT.end() != iWP; ++iWP ) {
      for ( auto iLLJJM : *(m_idxlljjmDRWP.find(iWP.cut().first).idxList()) ) {
        if ( iWP.cut().second(*(m_dileptondijetmets[iLLJJM].lljj->jj)) ) {
          iWP.idxList()->push_back(iLLJJM); // NOTE di-b working points may include an ETA cut
        }
      }
    }
    // Store objects according to CSVv2
    // First regular MET
    for ( std::size_t illbbmWP{0}; m_idxllbbmDRWP_PT.size() != illbbmWP; ++illbbmWP ) {
      indexlist_t& tagIdx = *(m_idxllbbmDRWP_tag.at(illbbmWP).idxList());
      tagIdx = *(m_idxllbbmDRWP_PT.at(illbbmWP).idxList());
      std::stable_sort(std::begin(tagIdx), std::end(tagIdx),
          MoreOf<index_t>([this] ( index_t iLLJJM ) {
            return m_dileptondijetmets[iLLJJM].lljj->jj->first->bDiscriminator(m_jetBTagName)
                 + m_dileptondijetmets[iLLJJM].lljj->jj->second->bDiscriminator(m_jetBTagName);
          }));
    }
  }

  LogDebug("ttW") << "End event.";
}

/**
 * Try to match `lepton` with an online object, using a deltaR and a deltaPt cut
 * Returns the index inside the HLTProducer collection, or -1 if no match is found.
 */
TTWAnalysis::sindex_t TTWAnalyzer::matchOfflineLepton( const TTWAnalysis::Lepton& lepton, const HLTProducer& hlt ) const
{
  LogDebug("ttW") << "Trying to match offline lepton: ";
  LogDebug("ttW") << "\tMuon? " << lepton.isMu() << " ; Pt: " << lepton.p4().Pt() << " ; Eta: " << lepton.p4().Eta() << " ; Phi: " << lepton.p4().Phi() << " ; E: " << lepton.p4().E() << std::endl;

  float min_dr = std::numeric_limits<float>::max();

  TTWAnalysis::sindex_t index{-1};
  for (size_t iHlt{0}; hlt.object_p4.size() != iHlt; ++iHlt) {
    float dr = VectorUtil::DeltaR(lepton.p4(), hlt.object_p4[iHlt]);
    float dpt_over_pt = std::abs(lepton.p4().Pt() - hlt.object_p4[iHlt].Pt()) / lepton.p4().Pt();

    if ( ( dr <= m_hltDRCut ) && ( dpt_over_pt <= m_hltDPtCut )
      && ( dr < min_dr ) )
    {
      min_dr = dr;
      index = iHlt;
    }
  }

  if (index != -1) {
    LogDebug("ttW") << "\033[32mMatched with online object:\033[00m";
    LogDebug("ttW") << "\tPDG Id: " << hlt.object_pdg_id[index] << " ; Pt: " << hlt.object_p4[index].Pt() << " ; Eta: " << hlt.object_p4[index].Eta() << " ; Phi: " << hlt.object_p4[index].Phi() << " ; E: " << hlt.object_p4[index].E();
    LogDebug("ttW") << "\tΔR: " << min_dr << " ; ΔPt / Pt: " << std::abs(lepton.p4().Pt() - hlt.object_p4[index].Pt()) / lepton.p4().Pt();
  } else {
    LogDebug("ttW") << "\033[31mNo match found\033[00m";
  }

  return index;
}

#include <FWCore/PluginManager/interface/PluginFactory.h>
DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, TTWAnalyzer, "ttw_analyzer");
