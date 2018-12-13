#pragma once

#include "cp3_llbb/Framework/interface/Category.h"

namespace TTWAnalysis {
class FakeLeptonCategory : public Category {
public:
  FakeLeptonCategory() : Category() {}
  virtual ~FakeLeptonCategory();
  virtual void configure(const edm::ParameterSet& conf) override;
  virtual bool event_in_category_pre_analyzers(const ProducersManager& producers) const override;
  virtual bool event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const override;
private:
  std::string m_ttwAnaName, m_selCat, m_jetCat;
  bool m_isEl;
  double m_maxMET, m_maxMT;
};
}
