#include "cp3_llbb/Framework/interface/MuonsProducer.h"
#include "cp3_llbb/Framework/interface/ElectronsProducer.h"
#include "cp3_llbb/Framework/interface/HLTProducer.h"
#include "cp3_llbb/TTWAnalysis/interface/TTWAnalyzer.h"

#include "cp3_llbb/TTWAnalysis/interface/TTWDileptonCategory.h"

#include <cassert>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>


using namespace TTWAnalysis;

namespace {
  std::vector<std::string> splitBySubString( const std::string& input, const std::string& delim )
  {
    std::vector<std::string> parts;
    boost::iter_split(parts, input, boost::algorithm::first_finder(delim));
    return parts;
  }
  std::pair<std::string,std::string> parseNamedCut( const std::string& val )
  {
    auto parts = splitBySubString(val, ":__:");
    assert(parts.size() == 2);
    return std::make_pair(parts[0], parts[1]);
  }
}

GDileptonCategory::GDileptonCategory()
 : Category()
 , m_llIsInCateg("")
{}

GDileptonCategory::~GDileptonCategory() {}

void GDileptonCategory::configure(const edm::ParameterSet& conf)
{
  // pre-ana in category
  m_nReqEl = conf.getParameter<unsigned int>("NElectrons");
  m_nReqMu = conf.getParameter<unsigned int>("NMuons");
  m_reqCh  = conf.getParameter<int>("Charge");
  // post-ana in category
  m_llIsInCateg = DiLeptonHybridCut(conf.getParameter<std::string>("Category"));
  // post-ana cuts
  m_postCuts.clear();
  m_postCuts.emplace_back(std::pair<std::string,DiLeptonCut>("Category", m_llIsInCateg));
  m_hltMatcher.addRegex(conf.getParameter<std::vector<std::string>>("HLT"));
  m_postCuts.emplace_back(std::pair<std::string,DiLeptonCut>("DiLeptonTriggerMatch", std::reference_wrapper<DiLeptonHLTMatch>(m_hltMatcher)));
  for ( const auto& otherCutStr : conf.getParameter<std::vector<std::string>>("Cuts") ) {
    auto parsedCut = parseNamedCut(otherCutStr);
    m_postCuts.emplace_back(std::pair<std::string,DiLeptonCut>(parsedCut.first, DiLeptonHybridCut(parsedCut.second)));
  }

  // construct IDIso combinations (and postfix strings)
  m_idIsoCombAndPostfix.clear();
  for ( const LepID::LepID& id1 : LepID::it ) {
    for ( const LepID::LepID& id2 : LepID::it ) {
      for ( const LepIso::LepIso& iso1 : LepIso::it ) {
        for ( const LepIso::LepIso& iso2 : LepIso::it ) {
          m_idIsoCombAndPostfix.push_back(std::make_pair(LepLepIDIso(id1, iso1, id2, iso2), "_"+LepLepIDIsoStr(id1, iso1, id2, iso2)));
        }
      }
    }
  }
}

bool GDileptonCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const
{
  unsigned int nPass;
  if ( m_nReqEl > 0 ) {
    const auto& electrons = producers.get<ElectronsProducer>("electrons");
    if ( m_reqCh != 0 ) {
      nPass = std::count_if(std::begin(electrons.charge), std::end(electrons.charge), [this] ( int8_t ch ) { return m_reqCh == ch; });
    } else {
      nPass = electrons.p4.size();
    }
    if ( nPass < m_nReqEl ) {
      return false;
    }
  }
  if ( m_nReqMu > 0 ) {
    const auto& muons = producers.get<MuonsProducer>("muons");
    if ( m_reqCh != 0 ) {
      nPass = std::count_if(std::begin(muons.charge), std::end(muons.charge), [this] ( int8_t ch ) { return m_reqCh == ch; });
    } else {
      nPass = muons.p4.size();
    }
    if ( nPass < m_nReqMu ) {
      return false;
    }
  }
  return true;
}

bool GDileptonCategory::diLeptonIsInCategory( const TTWAnalyzer* ttW, uint16_t dilepIdx ) const
{
  if ( m_llIsInCateg(ttW->diLeptons[dilepIdx]) ) {
    if ( m_reqCh == 0 ) {
      return true;
    } else {
      if ( ( m_reqCh == ttW->leptons[ttW->diLeptons[dilepIdx].lidxs.first].charge )
        && ( m_reqCh == ttW->leptons[ttW->diLeptons[dilepIdx].lidxs.second].charge ) ) {
        return true;
      }
    }
  }
  return false;
}

bool GDileptonCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const
{
  const TTWAnalyzer& ttW = analyzers.get<TTWAnalyzer>("ttW");

  for ( const auto& combAndPostfix : m_idIsoCombAndPostfix ) {
    const auto& catDiLeptons = ttW.diLeptons_IDIso[combAndPostfix.first];
    if ( ( ! catDiLeptons.empty() )
      && ( std::any_of(std::begin(catDiLeptons), std::end(catDiLeptons),
                  std::bind(&GDileptonCategory::diLeptonIsInCategory, this, &ttW, std::placeholders::_1)
                  ) ) ) {
      return true;
    }
  }
  return false;
}

void GDileptonCategory::register_cuts(CutManager& manager)
{
  for ( const auto& combAndPostfix : m_idIsoCombAndPostfix ) {
    for ( const auto& cut : m_postCuts ) {
      manager.new_cut(cut.first+combAndPostfix.second, cut.first+combAndPostfix.second);
    }
  }
}

void GDileptonCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {

  const TTWAnalyzer& ttW = analyzers.get<TTWAnalyzer>("ttW");
  m_hltMatcher.setHLTProd(&producers.get<HLTProducer>("hlt"));

  for ( const auto& combAndPostfix : m_idIsoCombAndPostfix ) {
    if ( ! ttW.diLeptons_IDIso[combAndPostfix.first].empty() ) {
      std::vector<uint16_t> sel_dileptons;
      sel_dileptons.reserve(ttW.diLeptons.size());
      std::copy_if(std::begin(ttW.diLeptons_IDIso[combAndPostfix.first]), std::end(ttW.diLeptons_IDIso[combAndPostfix.first]),
                   std::back_inserter(sel_dileptons), std::bind(&GDileptonCategory::diLeptonIsInCategory, this, &ttW, std::placeholders::_1));
      for ( const auto& cut : m_postCuts ) {
        std::vector<uint16_t> passSel_dileptons;
        passSel_dileptons.reserve(sel_dileptons.size());
        std::copy_if(std::begin(sel_dileptons), std::end(sel_dileptons), std::back_inserter(passSel_dileptons),
                    [&ttW,&cut] ( uint16_t idx ) -> bool { return cut.second(ttW.diLeptons[idx]); }
                );
        if ( ! passSel_dileptons.empty() ) {
          manager.pass_cut(cut.first+combAndPostfix.second);
        }
      }
    }
  }
}

bool GDileptonCategory::DiLeptonHLTMatch::operator() ( const TTWAnalysis::DiLepton& diLep ) const
{
  if ( ( ! m_hltProd ) || ( diLep.hlt_idxs.first == -1 ) || ( diLep.hlt_idxs.second == -1 ) ) { return false; }

  std::vector<std::string> mObj1, mObj2, mComm;
  for ( const auto& paths : m_regex ) {
    const auto& objp1 = m_hltProd->object_paths[diLep.hlt_idxs.first];
    std::copy_if(std::begin(objp1), std::end(objp1), std::back_inserter(mObj1),
        [&paths] ( const std::string& obj ) { return boost::regex_match(obj, paths); } );
    const auto& objp2 = m_hltProd->object_paths[diLep.hlt_idxs.second];
    std::copy_if(std::begin(objp2), std::end(objp2), std::back_inserter(mObj2),
        [&paths] ( const std::string& obj ) { return boost::regex_match(obj, paths); } );
  }
  std::set_intersection(std::begin(mObj1), std::end(mObj1), std::begin(mObj2), std::end(mObj2), std::back_inserter(mComm));
  return ! mComm.empty();
}
