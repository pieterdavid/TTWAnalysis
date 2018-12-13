#include "cp3_llbb/Framework/interface/Analyzer.h"
#include "cp3_llbb/TTWAnalysis/interface/TTWDileptonCategory.h"
#include "TTWFakeLeptonCategory.h"

class TTWCategoriesFactory : public Framework::Analyzer {
public:
  TTWCategoriesFactory(const std::string& name, const ROOT::TreeGroup& tree, const edm::ParameterSet& config)
    : Framework::Analyzer(name, tree, config) {}

  virtual ~TTWCategoriesFactory() {}

  void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const AnalyzersManager&, const CategoryManager&) override final {}

  virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet&) override;
};

void TTWCategoriesFactory::registerCategories(CategoryManager& manager, const edm::ParameterSet& config)
{
  for ( const std::string& catName : config.getParameterNames() ) {
    if ( ( catName.size() > 4 ) && ( "fake" == catName.substr(0,4) ) ) {
      manager.new_category<TTWAnalysis::FakeLeptonCategory>(catName, catName+" category", config.getParameter<edm::ParameterSet>(catName));
    } else {
      manager.new_category<TTWAnalysis::GDileptonCategory>(catName, catName+" category", config.getParameter<edm::ParameterSet>(catName));
    }
  }
}

#include <FWCore/PluginManager/interface/PluginFactory.h>
DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, TTWCategoriesFactory, "ttw_categories");
