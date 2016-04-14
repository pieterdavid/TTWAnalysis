#pragma once

#include <cp3_llbb/TTAnalysis/interface/Types.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>

namespace TTAnalysis {
  
  float DeltaEta(const myLorentzVector &v1, const myLorentzVector &v2);
  
  // Used by std::sort to sort jets according to decreasing b-tagging discriminant value
  class jetBTagDiscriminantSorter {
    
    public:
 
      // Either use indices to jets producer
      jetBTagDiscriminantSorter(const JetsProducer& jets, const std::string& taggerName): 
        m_jetsProducer(&jets),
        m_taggerName(taggerName),
        m_jetsArray(nullptr)
        {}

      // Or use indices to Jets array
      jetBTagDiscriminantSorter(const JetsProducer& jets, const std::string& taggerName, const std::vector<Jet>& tt_jets): 
        m_jetsProducer(&jets),
        m_taggerName(taggerName),
        m_jetsArray(&tt_jets)
        {}
      
      bool operator()(uint16_t idxJet1, uint16_t idxJet2){
        if(!m_jetsArray)
          return m_jetsProducer->getBTagDiscriminant(idxJet1, m_taggerName) > m_jetsProducer->getBTagDiscriminant(idxJet2, m_taggerName);
        else
          return m_jetsProducer->getBTagDiscriminant(m_jetsArray->at(idxJet1).idx, m_taggerName) > m_jetsProducer->getBTagDiscriminant(m_jetsArray->at(idxJet2).idx, m_taggerName);
      }
  
    private:
  
      const JetsProducer* const m_jetsProducer;
      const std::string m_taggerName;
      const std::vector<Jet>* const m_jetsArray;

  };

  // Used by std::sort to sort DiJets according to decreasing b-tagging discriminant value
  class diJetBTagDiscriminantSorter {
    
    public:
 
      // Either work directly on DiJet objects

      diJetBTagDiscriminantSorter(const JetsProducer& jets, const std::string& taggerName): 
        m_jetsProducer(jets), 
        m_taggerName(taggerName) 
        {}
      
      bool operator()(const DiJet& diJet1, const DiJet& diJet2){
        return 
          ( m_jetsProducer.getBTagDiscriminant(diJet1.idxs.first, m_taggerName) + m_jetsProducer.getBTagDiscriminant(diJet1.idxs.second, m_taggerName) ) > 
          ( m_jetsProducer.getBTagDiscriminant(diJet2.idxs.first, m_taggerName) + m_jetsProducer.getBTagDiscriminant(diJet2.idxs.second, m_taggerName) );
      }

      // Or work on vectors containing indices pointing to jets themselves
      // 1) Using DiJets
      diJetBTagDiscriminantSorter(const JetsProducer& jets, const std::string& taggerName, const std::vector<DiJet>& diJets):  
        m_jetsProducer(jets), 
        m_taggerName(taggerName), 
        m_diJets(&diJets), 
        m_diLepDiJets(nullptr), 
        m_diLepDiJetsMet(nullptr) 
        {}
      // 2) Using DiLepDiJets
      diJetBTagDiscriminantSorter(const JetsProducer& jets, const std::string& taggerName, const std::vector<DiLepDiJet>& diLepDiJets):  
        m_jetsProducer(jets), 
        m_taggerName(taggerName), 
        m_diJets(nullptr), 
        m_diLepDiJets(&diLepDiJets), 
        m_diLepDiJetsMet(nullptr) 
        {}
      // 3) Using DiLepDiJetsMet
      diJetBTagDiscriminantSorter(const JetsProducer& jets, const std::string& taggerName, const std::vector<DiLepDiJetMet>& diLepDiJetsMet):  
        m_jetsProducer(jets), 
        m_taggerName(taggerName), 
        m_diJets(nullptr), 
        m_diLepDiJets(nullptr), 
        m_diLepDiJetsMet(&diLepDiJetsMet) 
        {}
      
      bool operator()(const uint16_t idx1, const uint16_t idx2){
        
        if(m_diJets){
          return 
            ( m_jetsProducer.getBTagDiscriminant((*m_diJets)[idx1].idxs.first, m_taggerName) + m_jetsProducer.getBTagDiscriminant((*m_diJets)[idx1].idxs.second, m_taggerName) ) > 
            ( m_jetsProducer.getBTagDiscriminant((*m_diJets)[idx2].idxs.first, m_taggerName) + m_jetsProducer.getBTagDiscriminant((*m_diJets)[idx2].idxs.second, m_taggerName) );
        
        }else if(m_diLepDiJets){
          return 
            ( m_jetsProducer.getBTagDiscriminant((*m_diLepDiJets)[idx1].diJet->idxs.first, m_taggerName) + m_jetsProducer.getBTagDiscriminant((*m_diLepDiJets)[idx1].diJet->idxs.second, m_taggerName) ) > 
            ( m_jetsProducer.getBTagDiscriminant((*m_diLepDiJets)[idx2].diJet->idxs.first, m_taggerName) + m_jetsProducer.getBTagDiscriminant((*m_diLepDiJets)[idx2].diJet->idxs.second, m_taggerName) );
        
        }else if(m_diLepDiJetsMet){
          return 
            ( m_jetsProducer.getBTagDiscriminant((*m_diLepDiJetsMet)[idx1].diJet->idxs.first, m_taggerName) + m_jetsProducer.getBTagDiscriminant((*m_diLepDiJetsMet)[idx1].diJet->idxs.second, m_taggerName) ) > 
            ( m_jetsProducer.getBTagDiscriminant((*m_diLepDiJetsMet)[idx2].diJet->idxs.first, m_taggerName) + m_jetsProducer.getBTagDiscriminant((*m_diLepDiJetsMet)[idx2].diJet->idxs.second, m_taggerName) );
        
        }else{
          return false;
        }
      }
    
    private:
  
      const JetsProducer& m_jetsProducer;
      const std::string m_taggerName;
      const std::vector<DiJet>* m_diJets; 
      const std::vector<DiLepDiJet>* m_diLepDiJets; 
      const std::vector<DiLepDiJetMet>* m_diLepDiJetsMet; 
  };

}

