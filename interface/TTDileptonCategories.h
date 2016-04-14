#pragma once

#include <boost/regex.hpp>
#include <vector>
#include <string>

#include <cp3_llbb/Framework/interface/Category.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

namespace TTWAnalysis{

class DileptonCategory: public Category {
  public:
    virtual void configure(const edm::ParameterSet& conf) override {
      m_MllCutSF = conf.getUntrackedParameter<double>("MllCutSF", 20);
      m_MllCutDF = conf.getUntrackedParameter<double>("MllCutDF", 20);
      m_MllZVetoCutLow = conf.getUntrackedParameter<double>("MllZVetoCutLow", 86);
      m_MllZVetoCutHigh = conf.getUntrackedParameter<double>("MllZVetoCutHigh", 116);
      m_HLTDoubleMuon = conf.getUntrackedParameter<std::vector<std::string>>("HLTDoubleMuon");
      m_HLTDoubleEG = conf.getUntrackedParameter<std::vector<std::string>>("HLTDoubleEG");
      m_HLTMuonEG = conf.getUntrackedParameter<std::vector<std::string>>("HLTMuonEG");

      for(const auto& hlt: m_HLTDoubleMuon)
        m_HLTDoubleMuonRegex.push_back( boost::regex(hlt, boost::regex_constants::icase) );
      for(const auto& hlt: m_HLTDoubleEG)
        m_HLTDoubleEGRegex.push_back( boost::regex(hlt, boost::regex_constants::icase) );
      for(const auto& hlt: m_HLTMuonEG)
        m_HLTMuonEGRegex.push_back( boost::regex(hlt, boost::regex_constants::icase) );
    }

    DileptonCategory():
      baseStrCategory("Category"),
      baseStrExtraDiLeptonVeto("ExtraDiLeptonVeto"),
      baseStrDiLeptonTriggerMatch("DiLeptonTriggerMatch"),
      baseStrMllCut("Mll"),
      baseStrMllZVetoCut("MllZVeto"),
      baseStrDiLeptonIsOS("DiLeptonIsOS")
      {}

  protected:
    float m_MllCutSF, m_MllCutDF, m_MllZVetoCutLow, m_MllZVetoCutHigh;

    std::vector<std::string> m_HLTDoubleMuon;
    std::vector<std::string> m_HLTDoubleEG;
    std::vector<std::string> m_HLTMuonEG;

    std::string baseStrCategory;
    std::string baseStrExtraDiLeptonVeto;
    std::string baseStrDiLeptonTriggerMatch;
    std::string baseStrMllCut;
    std::string baseStrMllZVetoCut;
    std::string baseStrDiLeptonIsOS;

    std::vector<boost::regex> m_HLTDoubleMuonRegex;
    std::vector<boost::regex> m_HLTDoubleEGRegex;
    std::vector<boost::regex> m_HLTMuonEGRegex;

    enum class HLT { DoubleMuon, DoubleEG, MuonEG };

    // Check that the hlt objects at indices hltIdx1, hltIdx2 have fired at least one and the same 
    // of the trigger paths in the group specified by pathGroup.
    bool checkHLT(const HLTProducer& hlt, uint16_t hltIdx1, uint16_t hltIdx2, HLT pathGroup) const {
      const std::vector<boost::regex>* chosenPathGroup(nullptr);

      switch(pathGroup){
        
        case HLT::DoubleMuon:
          chosenPathGroup = &m_HLTDoubleMuonRegex;
          break;

        case HLT::DoubleEG:
          chosenPathGroup = &m_HLTDoubleEGRegex;
          break;

        case HLT::MuonEG:
          chosenPathGroup = &m_HLTMuonEGRegex;
          break;

        default:
          break;
      
      }

      if( !chosenPathGroup || hltIdx1 >= hlt.object_paths.size() || hltIdx2 >= hlt.object_paths.size() )
        return false;

      std::vector<std::string> matchedTriggersObj1, matchedTriggersObj2, commonMatchedTriggers;

      for(const auto& paths: *chosenPathGroup){
        
        for(const auto& objHlt: hlt.object_paths[hltIdx1]){
          if( boost::regex_match(objHlt, paths) )
            matchedTriggersObj1.push_back(objHlt);
        }
        
        for(const auto& objHlt: hlt.object_paths[hltIdx2]){
          if( boost::regex_match(objHlt, paths) )
            matchedTriggersObj2.push_back(objHlt);
        }
      
      }

      std::set_intersection( matchedTriggersObj1.begin(), matchedTriggersObj1.end(), matchedTriggersObj2.begin(), matchedTriggersObj2.end(), std::back_inserter(commonMatchedTriggers));

      return commonMatchedTriggers.size() > 0;
    }
      
};

class ElElCategory: public DileptonCategory {
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class ElMuCategory: public DileptonCategory {
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

 class MuElCategory: public DileptonCategory {
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

class MuMuCategory: public DileptonCategory {
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
};

}

