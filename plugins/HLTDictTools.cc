#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"

#include "cp3_llbb/Framework/interface/Types.h"
#include "cp3_llbb/TTWAnalysis/interface/DictTool.h"
#include "cp3_llbb/TTWAnalysis/interface/NewTypes.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include <cp3_llbb/Framework/interface/HLTService.h>

#include <Math/VectorUtil.h>
using namespace ROOT::Math;

#include <boost/regex.hpp>

namespace TTWAnalysis {

template<class Candidate>
class DictHLTMatchv2 : public DictTool<Candidate> {
public:
  DictHLTMatchv2(const edm::ParameterSet& config)
    : DictTool<Candidate>(config)
    , m_hltDRCut( config.getUntrackedParameter<double>("hltDRCut", std::numeric_limits<float>::max()) )
    , m_hltDPtCut( config.getUntrackedParameter<double>("hltDPtCut", std::numeric_limits<float>::max()) )
    //, m_hasPrinted(0)
  {
    if (config.exists("triggers")) {
      m_hlt_service.reset(new HLTService(config.getUntrackedParameter<edm::FileInPath>("triggers").fullPath()));
      std::cout << "HLT configuration: " << std::endl;
      m_hlt_service->print();
      std::cout << std::endl;
    }

    const auto& selsConfig = config.getParameter<edm::ParameterSet>("Selections_PathRegex");
    for ( const auto& iName : selsConfig.getParameterNames() ) {
      const auto& sels = selsConfig.getParameter<std::vector<std::string>>(iName);
      std::vector<boost::regex> patterns; patterns.reserve(sels.size());
      std::transform(std::begin(sels), std::end(sels), std::back_inserter(patterns),
          [] ( const std::string& ireg ) { return boost::regex(ireg, boost::regex_constants::icase); });
      m_selections_pathregex.emplace(iName, patterns);
    }
    const auto& filtersConfig = config.getParameter<edm::ParameterSet>("Selections_Filter");
    for ( const auto& iName : filtersConfig.getParameterNames() ) {
      const auto& filterName = filtersConfig.getParameter<std::string>(iName);
      m_selections_filter.emplace(iName, filtersConfig.getParameter<std::string>(iName));
      m_filters.push_back(filterName);
    }
  }
  virtual ~DictHLTMatchv2() {}

  virtual void doConsumes(const edm::ParameterSet& config, edm::ConsumesCollector&& collector) override {
    m_hlt_token = collector.consumes<edm::TriggerResults>(config.getUntrackedParameter<edm::InputTag>("hlt", edm::InputTag("TriggerResults", "", "HLT")));
    // m_prescales_token = collector.consumes<pat::PackedTriggerPrescales>(config.getUntrackedParameter<edm::InputTag>("prescales", edm::InputTag("patTrigger")));
    m_trigger_objects_token = collector.consumes<pat::TriggerObjectStandAloneCollection>(config.getUntrackedParameter<edm::InputTag>("objects", edm::InputTag("selectedPatTrigger")));
  }

  struct SelectedObject {
    pat::TriggerObjectStandAlone obj;
    std::vector<std::string> paths;
    std::vector<std::string> filters;
    SelectedObject(const pat::TriggerObjectStandAlone& theObj,
                   std::vector<std::string>&& thePaths,
                   std::vector<std::string>&& theFilters)
      : obj(theObj)
      , paths(thePaths)
      , filters(theFilters)
    {}

    bool hasMatchingPath( const boost::regex& sel ) const
    {
      return std::any_of(std::begin(paths), std::end(paths),
        [&sel] ( const std::string& trig ) {
          return boost::regex_match(trig, sel);
        });
    }
    bool kinMatch( const Lepton& lepton, float maxDR, float maxDPT ) const
    {
      return (( VectorUtil::DeltaR(lepton.p4(), obj.p4()) <= maxDR )
           && ( std::abs(lepton.p4().Pt() - obj.p4().Pt()) / lepton.p4().Pt() <= maxDPT ) );
    }
    bool hasFilter( const std::string& filterLabel ) const
    {
      return std::end(filters) != std::find(std::begin(filters), std::end(filters), filterLabel);
    }
  };

  virtual Dict evaluate(const Candidate& cand,
      const edm::Event* event, const edm::EventSetup* /**/,
      const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const override;

private:
  std::map<std::string,std::vector<boost::regex>> m_selections_pathregex;
  std::map<std::string,std::string> m_selections_filter;
  // Tokens
  edm::EDGetTokenT<edm::TriggerResults> m_hlt_token;
  // edm::EDGetTokenT<pat::PackedTriggerPrescales> m_prescales_token;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> m_trigger_objects_token;
  // Service
  std::shared_ptr<HLTService> m_hlt_service;

  const float m_hltDRCut, m_hltDPtCut;

  mutable edm::EventID m_lastEvent;
  mutable std::vector<std::string> m_paths;
  // mutable std::vector<uint16_t> m_prescales;
  std::vector<std::string> m_filters;
  mutable std::vector<SelectedObject> m_triggerObjects;

  //mutable std::size_t m_hasPrinted;

  void readTriggerNamesAndPrescales( const edm::Event* event ) const
  {
    m_paths.clear(); // m_prescales.clear();

    edm::Handle<edm::TriggerResults> hlt;
    event->getByToken(m_hlt_token, hlt);
    // edm::Handle<pat::PackedTriggerPrescales> prescales;
    //event->getByToken(m_prescales_token, prescales);

    const edm::TriggerNames& triggerNames = event->triggerNames(*hlt);

    bool filter = m_hlt_service.get() != nullptr;
    const HLTService::PathVector* valid_paths = nullptr;
    if (filter) {
      valid_paths = &m_hlt_service->getPaths(event->id().run());
    }

    // get paths and prescales
    for ( size_t i = 0 ; hlt->size() != i; ++i ) {
      if (hlt->accept(i)) {
        std::string triggerName = triggerNames.triggerName(i);
        if (triggerName == "HLTriggerFinalPath")
          continue; // This one is pretty useless...
        if (triggerName[0] == 'A')
          continue; // Remove AlCa HLT paths

        if ( (!filter) || std::any_of(valid_paths->begin(), valid_paths->end(),
          [&triggerName] ( const HLTService::PathName& regex ) { return boost::regex_match(triggerName, regex); } ) )
        {
          m_paths.push_back(triggerName);
          // m_prescales.push_back(prescales.isValid() ? prescales->getPrescaleForIndex(i) : 1.);
        }
      }
    }
    std::sort(std::begin(m_paths), std::end(m_paths));
    //if ( ( ! m_paths.empty() ) && ( m_hasPrinted < 10 ) ) {
    //  ++m_hasPrinted;
    //  std::cout << "Trigger paths matching any of the selected: ";
    //  for ( const auto& iPath : m_paths ) {
    //    std::cout << iPath << ", ";
    //  }
    //  std::cout << std::endl;
    //}
  }

  void collectTriggerObjects( const edm::Event* event ) const
  {
    m_triggerObjects.clear();

    edm::Handle<edm::TriggerResults> hlt;
    event->getByToken(m_hlt_token, hlt);
    const bool filter = m_hlt_service.get() != nullptr;
    const edm::TriggerNames& triggerNames = event->triggerNames(*hlt);
    edm::Handle<pat::TriggerObjectStandAloneCollection> objects;
    event->getByToken(m_trigger_objects_token, objects);

    if ( objects.isValid() ) {
      for ( pat::TriggerObjectStandAlone obj : *objects ) {
        obj.unpackPathNames(triggerNames);
        //obj.unpackNamesAndLabels(event, hlt); // for newer cmssw

        std::vector<std::string> object_paths = obj.pathNames(false);
        const std::vector<std::string>& object_filters = obj.filterLabels();

        // Check if this object has triggered at least one of the paths or filters we are interesting in
        if ( (!filter) || std::any_of(std::begin(object_paths), std::end(object_paths),
          [this] ( const std::string& path ) {
            return std::end(m_paths) != std::find(std::begin(m_paths), std::end(m_paths), path);
          } )          || std::any_of(std::begin(object_filters), std::end(object_filters),
          [this] ( const std::string& filter ) {
            return std::end(m_filters) != std::find(std::begin(m_filters), std::end(m_filters), filter);
          } ) )
        {
          std::sort(std::begin(object_paths), std::end(object_paths));

          std::vector<std::string> sel_paths;
          std::set_intersection(m_paths.begin(), m_paths.end(), object_paths.begin(), object_paths.end(), std::back_inserter(sel_paths));

          std::vector<std::string> sel_filters;
          std::set_intersection(std::begin(m_filters), std::end(m_filters), std::begin(object_filters), std::end(object_filters), std::back_inserter(sel_filters));

          m_triggerObjects.emplace_back(std::move(obj), std::move(sel_paths), std::move(sel_filters));
        }
      }
    }
    //std::cout << "Collected trigger objects: " << m_triggerObjects.size() << std::endl;
    //for ( const auto& trigObj : m_triggerObjects ) {
    //  std::cout << "  - (PT=" << trigObj.obj.pt() << ",ETA=" << trigObj.obj.eta() << ",PHI=" << trigObj.obj.phi() << ") PATHS ";
    //  for ( const auto& trigPath : trigObj.paths   ) { std::cout << trigPath << ", "; }
    //  std::cout << " FILTERS ";
    //  for ( const auto& trigFilt : trigObj.filters ) { std::cout << trigFilt << ", "; }
    //  std::cout << std::endl;
    //}
  }
};

bool matchHLTCandidate( const Lepton& lepton, const pat::TriggerObjectStandAlone& trigObj, float maxDR, float maxDPT )
{
  return (( VectorUtil::DeltaR(lepton.p4(), trigObj.p4()) <= maxDR )
       && ( std::abs(lepton.p4().Pt() - trigObj.p4().Pt()) / lepton.p4().Pt() <= maxDPT ) );

}

template<class Candidate>
bool hasHLTCandidateMatchPath( const Candidate& cand, const std::vector<typename DictHLTMatchv2<Candidate>::SelectedObject>& trigObjects, const std::vector<boost::regex>& selections, float maxDR, float maxDPT );

template<>
bool hasHLTCandidateMatchPath<Lepton>( const Lepton& lepton, const std::vector<typename DictHLTMatchv2<Lepton>::SelectedObject>& trigObjects, const std::vector<boost::regex>& selections, float maxDR, float maxDPT )
{
  return std::any_of(std::begin(trigObjects), std::end(trigObjects),
      [&,maxDR,maxDPT] ( const DictHLTMatchv2<Lepton>::SelectedObject& obj ) {
        return std::any_of(std::begin(selections), std::end(selections),
          [&,maxDR,maxDPT] ( const boost::regex& sel ) {
            return obj.hasMatchingPath(sel) && obj.kinMatch(lepton, maxDR, maxDPT);
          });
        });
}

template<>
bool hasHLTCandidateMatchPath<DiLepton>( const DiLepton& dilepton, const std::vector<typename DictHLTMatchv2<DiLepton>::SelectedObject>& trigObjects, const std::vector<boost::regex>& selections, float maxDR, float maxDPT )
{
  std::vector<std::string> mObj1, mObj2, mComm;
  for ( const auto& sel : selections ) {
    for ( const auto& obj : trigObjects ) {
      if ( dilepton.first  && obj.hasMatchingPath(sel) && obj.kinMatch(*(dilepton.first ), maxDR, maxDPT ) ) {
        std::copy_if(std::begin(obj.paths), std::end(obj.paths), std::back_inserter(mObj1),
            [&sel] ( const std::string& oSel ) { return boost::regex_match(oSel, sel); } );
      }
      if ( dilepton.second && obj.hasMatchingPath(sel) && obj.kinMatch(*(dilepton.second), maxDR, maxDPT ) ) {
        std::copy_if(std::begin(obj.paths), std::end(obj.paths), std::back_inserter(mObj2),
            [&sel] ( const std::string& oSel ) { return boost::regex_match(oSel, sel); } );
      }
    }
  }
  std::set_intersection(std::begin(mObj1), std::end(mObj1), std::begin(mObj2), std::end(mObj2), std::back_inserter(mComm));
  //std::cout << "HLT-v2 dilepton matching for ";
  //for ( const auto& res : selections ) {
  //  std::cout << res.str() << ", ";
  //}
  //if ( ( ! mObj1.empty() ) || ( ! mObj2.empty() ) || ( ! mComm.empty() ) ) {
  //  std::cout << "First lepton matches ";
  //  for ( const auto& sel : mObj1 ) { std::cout << sel << ", "; }
  //  std::cout << " second lepton matches ";
  //  for ( const auto& sel : mObj2 ) { std::cout << sel << ", "; }
  //  std::cout << " --> common matches are ";
  //  for ( const auto& sel : mComm ) { std::cout << sel << ", "; }
  //  std::cout << std::endl;
  //} else {
  //  std::cout << "  NO match" << std::endl;
  //}
  return ! mComm.empty();
}

template<class Candidate>
bool hasHLTCandidateMatchFilter( const Candidate& cand, const std::vector<typename DictHLTMatchv2<Candidate>::SelectedObject>& trigObjects, const std::string& filter, float maxDR, float maxDPT );

template<>
bool hasHLTCandidateMatchFilter<Lepton>( const Lepton& lepton, const std::vector<typename DictHLTMatchv2<Lepton>::SelectedObject>& trigObjects, const std::string& filter, float maxDR, float maxDPT )
{
  return std::any_of(std::begin(trigObjects), std::end(trigObjects),
      [&,maxDR,maxDPT] ( const DictHLTMatchv2<Lepton>::SelectedObject& obj ) {
        return obj.hasFilter(filter) && obj.kinMatch(lepton, maxDR, maxDPT);
      });
}
template<>
bool hasHLTCandidateMatchFilter<DiLepton>( const DiLepton& dilepton, const std::vector<typename DictHLTMatchv2<DiLepton>::SelectedObject>& trigObjects, const std::string& filter, float maxDR, float maxDPT )
{ // FIXME untested, but need implementation to compile
  return  ( std::any_of(std::begin(trigObjects), std::end(trigObjects),
      [&,maxDR,maxDPT] ( const DictHLTMatchv2<DiLepton>::SelectedObject& obj ) {
        return obj.hasFilter(filter) && obj.kinMatch(*(dilepton.first ), maxDR, maxDPT);
      }) && std::any_of(std::begin(trigObjects), std::end(trigObjects),
      [&,maxDR,maxDPT] ( const DictHLTMatchv2<DiLepton>::SelectedObject& obj ) {
        return obj.hasFilter(filter) && obj.kinMatch(*(dilepton.second), maxDR, maxDPT);
      }) );
}

template<class Candidate>
Dict DictHLTMatchv2<Candidate>::evaluate(const Candidate& cand,
    const edm::Event* event, const edm::EventSetup* /**/,
    const ProducersManager* /**/, const AnalyzersManager* /**/, const CategoryManager* /**/) const
{
  if ( event && ( m_lastEvent != event->id() ) ) {
    readTriggerNamesAndPrescales(event);
    collectTriggerObjects(event);
    m_lastEvent = event->id();
  }

  Dict ret{};
  for ( const auto& selection : m_selections_pathregex ) {
    ret.add(selection.first, hasHLTCandidateMatchPath<Candidate>(cand, m_triggerObjects, selection.second, m_hltDRCut, m_hltDPtCut));
  }
  for ( const auto& filter : m_selections_filter ) {
    ret.add(filter.first, hasHLTCandidateMatchFilter<Candidate>(cand, m_triggerObjects, filter.second, m_hltDRCut, m_hltDPtCut));
  }

  return ret;
}

}

DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::Lepton  >::factory, TTWAnalysis::DictHLTMatchv2<TTWAnalysis::Lepton>  , "ttw_leptonHLTMatchv2");
DEFINE_EDM_PLUGIN(TTWAnalysis::DictTool<TTWAnalysis::DiLepton>::factory, TTWAnalysis::DictHLTMatchv2<TTWAnalysis::DiLepton>, "ttw_dileptonHLTMatchv2");
