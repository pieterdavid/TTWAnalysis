#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>

#include <cp3_llbb/TTWAnalysis/interface/TTDileptonCategories.h>
#include <cp3_llbb/TTWAnalysis/interface/TTAnalyzer.h>

#include <cp3_llbb/TTWAnalysis/interface/Types.h>
#include <cp3_llbb/TTWAnalysis/interface/Indices.h>

using namespace TTWAnalysis;

// ***** ***** *****
// Dilepton El-El category
// ***** ***** *****
bool ElElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
  return electrons.p4.size() >= 2;
}

bool ElElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");

  // It at least one DiLepton of highest Pt and of type ElEl among all ID pairs is found, keep event in this category

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
          if(tt.diLeptons_IDIso[comb].size() >= 1) {
            if( tt.diLeptons[ tt.diLeptons_IDIso[comb][0] ].isElEl )
              return true;
          }
        
        }
      }
    }
  }

  return false;
}

void ElElCategory::register_cuts(CutManager& manager) {
  
  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          
          manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
          manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
          manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
          manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
          manager.new_cut(baseStrMllZVetoCut + postFix, baseStrMllZVetoCut + postFix);
          manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);

        }
      }
    }
  }

}

void ElElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
    
          if(tt.diLeptons_IDIso[comb].size() >= 1) {
            const DiLepton& m_diLepton = tt.diLeptons[ tt.diLeptons_IDIso[comb][0] ];
            
            if(m_diLepton.isElEl) {
              manager.pass_cut(baseStrCategory + postFix);

              if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
                // We have fired a trigger. Now, check that it is actually a DoubleEG trigger
                if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::DoubleEG) )
                  manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
              }
              
              if(m_diLepton.p4.M() > m_MllCutSF)
                manager.pass_cut(baseStrMllCut + postFix);
              
              if(m_diLepton.p4.M() < m_MllZVetoCutLow || m_diLepton.p4.M() > m_MllZVetoCutHigh)
                manager.pass_cut(baseStrMllZVetoCut + postFix);
              
              if(m_diLepton.isOS)
                manager.pass_cut(baseStrDiLeptonIsOS + postFix);
            }
          }
          
          // For electrons, in principe only veto using VetoID.
          // But since the user can access any cut he wants, he can take the IDVV_IsoWhatever cut.
          if(tt.diLeptons_IDIso[comb].size() >= 2) { 
            manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
          }

        }
      }
    }
  }

}

// ***** ***** *****
// Dilepton El-Mu category
// ***** ***** *****
bool ElMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
  const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
  return electrons.p4.size() >= 1 && muons.p4.size() >= 1;
}

bool ElMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");

  // It at least one DiLepton of highest Pt and of type ElMu among all ID pairs is found, keep event in this category

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
          if(tt.diLeptons_IDIso[comb].size() >= 1) {
            if( tt.diLeptons[ tt.diLeptons_IDIso[comb][0] ].isElMu )
              return true;
          }
        
        }
      }
    }
  }

  return false;
}

void ElMuCategory::register_cuts(CutManager& manager) {
  
  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          
          manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
          manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
          manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
          manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
          manager.new_cut(baseStrMllZVetoCut + postFix, baseStrMllZVetoCut + postFix);
          manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);

        }
      }
    }
  }

}

void ElMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
    
          if(tt.diLeptons_IDIso[comb].size() >= 1) {
            const DiLepton& m_diLepton = tt.diLeptons[ tt.diLeptons_IDIso[comb][0] ];
            
            if(m_diLepton.isElMu) {
              manager.pass_cut(baseStrCategory + postFix);

              if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
                // We have fired a trigger. Now, check that it is actually a MuonEG trigger
                if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::MuonEG) )
                  manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
              }
              
              if(m_diLepton.p4.M() > m_MllCutDF)
                manager.pass_cut(baseStrMllCut + postFix);
              
              if(m_diLepton.p4.M() < m_MllZVetoCutLow || m_diLepton.p4.M() > m_MllZVetoCutHigh)
                manager.pass_cut(baseStrMllZVetoCut + postFix);
              
              if(m_diLepton.isOS)
                manager.pass_cut(baseStrDiLeptonIsOS + postFix);
            }
          }
          
          // For electrons, in principe only veto using VetoID.
          // But since the user can access any cut he wants, he can take the IDVV_IsoWhatever cut.
          if(tt.diLeptons_IDIso[comb].size() >= 2) { 
            manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
          }

        }
      }
    }
  }

}

// ***** ***** *****
// Dilepton Mu-El category
// ***** ***** *****
bool MuElCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const ElectronsProducer& electrons = producers.get<ElectronsProducer>("electrons");
  const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
  return electrons.p4.size() >= 1 && muons.p4.size() >= 1;
}

bool MuElCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");

  // It at least one DiLepton of highest Pt and of type MuEl among all ID pairs is found, keep event in this category

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
          if(tt.diLeptons_IDIso[comb].size() >= 1) {
            if( tt.diLeptons[ tt.diLeptons_IDIso[comb][0] ].isMuEl )
              return true;
          }
        
        }
      }
    }
  }

  return false;
}

void MuElCategory::register_cuts(CutManager& manager) {
  
  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          
          manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
          manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
          manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
          manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
          manager.new_cut(baseStrMllZVetoCut + postFix, baseStrMllZVetoCut + postFix);
          manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);

        }
      }
    }
  }

}

void MuElCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
    
          if(tt.diLeptons_IDIso[comb].size() >= 1) {
            const DiLepton& m_diLepton = tt.diLeptons[ tt.diLeptons_IDIso[comb][0] ];
            
            if(m_diLepton.isMuEl) {
              manager.pass_cut(baseStrCategory + postFix);

              if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
                // We have fired a trigger. Now, check that it is actually a MuonEG trigger
                if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::MuonEG) )
                  manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
              }
              
              if(m_diLepton.p4.M() > m_MllCutDF)
                manager.pass_cut(baseStrMllCut + postFix);
              
              if(m_diLepton.p4.M() < m_MllZVetoCutLow || m_diLepton.p4.M() > m_MllZVetoCutHigh)
                manager.pass_cut(baseStrMllZVetoCut + postFix);
              
              if(m_diLepton.isOS)
                manager.pass_cut(baseStrDiLeptonIsOS + postFix);
            }
          }
          
          // For electrons, in principe only veto using VetoID.
          // But since the user can access any cut he wants, he can take the IDVV_IsoWhatever cut.
          if(tt.diLeptons_IDIso[comb].size() >= 2) { 
            manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
          }

        }
      }
    }
  }

}

// ***** ***** *****
// Dilepton Mu-Mu category
// ***** ***** *****
bool MuMuCategory::event_in_category_pre_analyzers(const ProducersManager& producers) const {
  const MuonsProducer& muons = producers.get<MuonsProducer>("muons");
  return muons.p4.size() >= 2;
}

bool MuMuCategory::event_in_category_post_analyzers(const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");

  // It at least one DiLepton of highest Pt and of type MuMu among all ID pairs is found, keep event in this category

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
          if(tt.diLeptons_IDIso[comb].size() >= 1) {
            if( tt.diLeptons[ tt.diLeptons_IDIso[comb][0] ].isMuMu )
              return true;
          }
        
        }
      }
    }
  }

  return false;
}

void MuMuCategory::register_cuts(CutManager& manager) {
  
  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          
          manager.new_cut(baseStrCategory + postFix, baseStrCategory + postFix);
          manager.new_cut(baseStrExtraDiLeptonVeto + postFix, baseStrExtraDiLeptonVeto + postFix);
          manager.new_cut(baseStrDiLeptonTriggerMatch + postFix, baseStrDiLeptonTriggerMatch + postFix);
          manager.new_cut(baseStrMllCut + postFix, baseStrMllCut + postFix);
          manager.new_cut(baseStrMllZVetoCut + postFix, baseStrMllZVetoCut + postFix);
          manager.new_cut(baseStrDiLeptonIsOS + postFix, baseStrDiLeptonIsOS + postFix);

        }
      }
    }
  }

}

void MuMuCategory::evaluate_cuts_post_analyzers(CutManager& manager, const ProducersManager& producers, const AnalyzersManager& analyzers) const {
  
  const TTAnalyzer& tt = analyzers.get<TTAnalyzer>("tt");
  const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

  for(const LepID::LepID& id1: LepID::it) {
    for(const LepID::LepID& id2: LepID::it) {
      for(const LepIso::LepIso& iso1: LepIso::it) {
        for(const LepIso::LepIso& iso2: LepIso::it) {
          
          std::string postFix("_");
          postFix += LepLepIDIsoStr(id1, iso1, id2, iso2);
          uint16_t comb = LepLepIDIso(id1, iso1, id2, iso2);
    
          if(tt.diLeptons_IDIso[comb].size() >= 1) {
            const DiLepton& m_diLepton = tt.diLeptons[ tt.diLeptons_IDIso[comb][0] ];
            
            if(m_diLepton.isMuMu) {
              manager.pass_cut(baseStrCategory + postFix);

              if(m_diLepton.hlt_idxs.first >= 0 && m_diLepton.hlt_idxs.second >= 0){
                // We have fired a trigger. Now, check that it is actually a DoubleMuon trigger
                if( checkHLT(hlt, m_diLepton.hlt_idxs.first, m_diLepton.hlt_idxs.second, HLT::DoubleMuon) )
                  manager.pass_cut(baseStrDiLeptonTriggerMatch + postFix);
              }
              
              if(m_diLepton.p4.M() > m_MllCutSF)
                manager.pass_cut(baseStrMllCut + postFix);
              
              if(m_diLepton.p4.M() < m_MllZVetoCutLow || m_diLepton.p4.M() > m_MllZVetoCutHigh)
                manager.pass_cut(baseStrMllZVetoCut + postFix);
              
              if(m_diLepton.isOS)
                manager.pass_cut(baseStrDiLeptonIsOS + postFix);
            }
          }
          
          if(tt.diLeptons_IDIso[comb].size() >= 2) { 
            manager.pass_cut(baseStrExtraDiLeptonVeto + postFix);
          }

        }
      }
    }
  }

}

