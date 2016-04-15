#pragma once

#include <utility>
#include <vector>
#include <string>
#include <functional>
#include <boost/regex.hpp>

#include "cp3_llbb/Framework/interface/Category.h"
#include "cp3_llbb/TTWAnalysis/interface/Types.h"

// fwd
class HLTProducer;
class TTWAnalyzer;

namespace TTWAnalysis {

class GDileptonCategory : public Category {
  public:
    GDileptonCategory();
    virtual ~GDileptonCategory();
    virtual void configure(const edm::ParameterSet& conf) override;
  private:
    typedef std::function<bool(const TTWAnalysis::DiLepton&)> DiLeptonCut;
    typedef StringCutObjectSelector<TTWAnalysis::DiLepton> DiLeptonHybridCut;
    // for in_category_pre_analyzer
    unsigned int m_nReqEl, m_nReqMu;
    int m_reqCh;
    // for in_category_post_analyzer
    DiLeptonHybridCut m_llIsInCateg;
    // registered cuts
    std::vector<std::pair<std::string,std::function<bool(const TTWAnalysis::DiLepton&)>>> m_postCuts;
    // cache ID-Iso combinations
    std::vector<std::pair<uint16_t,std::string>> m_idIsoCombAndPostfix;
  public:
    virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
    virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
    virtual void register_cuts(CutManager& manager) override;
    virtual void evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
  private:
    /*
     * Helper functor for HLT candidate matching
     */
    class DiLeptonHLTMatch {
      public:
        DiLeptonHLTMatch() : m_hltProd(nullptr) {}

        const HLTProducer* hltProd() const { return m_hltProd; }
        void setHLTProd( const HLTProducer* hltProd ) { m_hltProd = hltProd; }
        void addRegex( const std::vector<std::string>& regex )
        {
          std::transform(std::begin(regex), std::end(regex), std::back_inserter(m_regex),
              [] ( const std::string& ireg ) { return boost::regex(ireg, boost::regex_constants::icase); } );
        }

        bool operator() ( const TTWAnalysis::DiLepton& diLep ) const;

      private:
        const HLTProducer* m_hltProd;
        std::vector<boost::regex> m_regex;
    };
    mutable DiLeptonHLTMatch m_hltMatcher;

    bool diLeptonIsInCategory( const TTWAnalyzer* ttW, uint16_t dilepIdx ) const;
};

}
