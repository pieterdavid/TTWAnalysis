#include <cp3_llbb/TTWAnalysis/interface/Defines.h>
#include <cp3_llbb/TTWAnalysis/interface/Types.h>
#include <cp3_llbb/TTWAnalysis/interface/Tools.h>
#include <cp3_llbb/TTWAnalysis/interface/GenStatusFlags.h>
#include <cp3_llbb/TTWAnalysis/interface/TTAnalyzer.h>
#include <cp3_llbb/TTWAnalysis/interface/TTDileptonCategories.h>

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/ElectronsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/METProducer.h>
#include <cp3_llbb/Framework/interface/HLTProducer.h>
#include <cp3_llbb/Framework/interface/GenParticlesProducer.h>

#include <Math/PtEtaPhiE4D.h>
#include <Math/LorentzVector.h>
#include <Math/VectorUtil.h>

// To access VectorUtil::DeltaR() more easily
using namespace ROOT::Math;

using namespace TTWAnalysis;

float TTWAnalysis::DeltaEta(const myLorentzVector& v1, const myLorentzVector& v2) {
  return std::abs(v1.Eta() - v2.Eta());
}

void TTAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup, const ProducersManager& producers, const AnalyzersManager& analyzers, const CategoryManager& categories) {
  
  #ifdef _TT_DEBUG_
    std::cout << "Begin event." << std::endl;
  #endif

  // Initizalize vectors depending on IDs/WPs to the right lengths
  // Only a resize() is needed (and no assign()), since TreeWrapper clears the vectors after each event.

  electrons_IDIso.resize( LepID::Count * LepIso::Count );
  muons_IDIso.resize( LepID::Count * LepIso::Count );
  leptons_IDIso.resize( LepID::Count * LepIso::Count );

  diLeptons_IDIso.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count );
  
  selJets_selID_DRCut.resize( LepID::Count * LepIso::Count );
  selBJets_DRCut_BWP_PtOrdered.resize( LepID::Count * LepIso::Count * BWP::Count );
  selBJets_DRCut_BWP_CSVv2Ordered.resize( LepID::Count * LepIso::Count * BWP::Count );

  diJets_DRCut.resize( LepID::Count * LepIso::Count );
  diBJets_DRCut_BWP_PtOrdered.resize( LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  diBJets_DRCut_BWP_CSVv2Ordered.resize( LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  
  diLepDiJets_DRCut.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count );
  diLepDiBJets_DRCut_BWP_PtOrdered.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  diLepDiBJets_DRCut_BWP_CSVv2Ordered.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  
  diLepDiJetsMet_DRCut.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count );
  diLepDiBJetsMet_DRCut_BWP_PtOrdered.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  diLepDiBJetsMet_DRCut_BWP_CSVv2Ordered.resize( LepID::Count * LepIso::Count * LepID::Count * LepIso::Count * BWP::Count * BWP::Count );
  
  gen_matched_b.resize( LepID::Count * LepIso::Count , -1);
  gen_matched_b_beforeFSR.resize( LepID::Count * LepIso::Count , -1);
  gen_matched_bbar.resize( LepID::Count * LepIso::Count , -1);
  gen_matched_bbar_beforeFSR.resize( LepID::Count * LepIso::Count , -1);
  gen_b_deltaR.resize( LepID::Count * LepIso::Count );
  gen_b_beforeFSR_deltaR.resize( LepID::Count * LepIso::Count );
  gen_bbar_deltaR.resize( LepID::Count * LepIso::Count );
  gen_bbar_beforeFSR_deltaR.resize( LepID::Count * LepIso::Count );

  ///////////////////////////
  //       ELECTRONS       //
  ///////////////////////////

  #ifdef _TT_DEBUG_
    std::cout << "Electrons" << std::endl;
  #endif

  const ElectronsProducer& electrons = producers.get<ElectronsProducer>(m_electrons_producer);

  for(uint16_t ielectron = 0; ielectron < electrons.p4.size(); ielectron++){
    if( electrons.p4[ielectron].Pt() > m_electronPtCut && std::abs(electrons.p4[ielectron].Eta()) < m_electronEtaCut ){
      
      Lepton m_lepton(
          electrons.p4[ielectron], 
          ielectron, 
          electrons.charge[ielectron], 
          true, false,
          electrons.ids[ielectron][m_electronVetoIDName],
          electrons.ids[ielectron][m_electronLooseIDName],
          electrons.ids[ielectron][m_electronMediumIDName],
          electrons.ids[ielectron][m_electronTightIDName],
          electrons.relativeIsoR03_withEA[ielectron]
      );
      
      for(const LepID::LepID& id: LepID::it){
        for(const LepIso::LepIso& iso: LepIso::it){ // loop over Iso not really needed since not considered for electrons (for the moment)
          uint16_t idx = LepIDIso(id, iso);
          if( m_lepton.ID[id] && m_lepton.iso[iso] )
            electrons_IDIso[idx].push_back(ielectron);
        }
      }
      
      leptons.push_back(m_lepton);
    }
  }

  ///////////////////////////
  //       MUONS           //
  ///////////////////////////
  
  #ifdef _TT_DEBUG_
    std::cout << "Muons" << std::endl;
  #endif

  const MuonsProducer& muons = producers.get<MuonsProducer>(m_muons_producer);

  for(uint16_t imuon = 0; imuon < muons.p4.size(); imuon++){
    if(muons.p4[imuon].Pt() > m_muonPtCut && std::abs(muons.p4[imuon].Eta()) < m_muonEtaCut ){
      
      Lepton m_lepton(
          muons.p4[imuon], 
          imuon,
          muons.charge[imuon], 
          false, true,
          muons.isLoose[imuon], // isVeto => for muons, re-use isLoose
          muons.isLoose[imuon],
          muons.isMedium[imuon],
          muons.isTight[imuon],
          muons.relativeIsoR04_deltaBeta[imuon],
          muons.relativeIsoR04_deltaBeta[imuon] < m_muonLooseIsoCut,
          muons.relativeIsoR04_deltaBeta[imuon] < m_muonTightIsoCut
      );

      for(const LepID::LepID& id: LepID::it){
        for(const LepIso::LepIso& iso: LepIso::it){
          uint16_t idx = LepIDIso(id, iso);
          if( m_lepton.ID[id] && m_lepton.iso[iso] )
            muons_IDIso[idx].push_back(imuon);
        }
      }

      leptons.push_back(m_lepton);
    }
  }

  // Sort the leptons vector according to Pt
  std::sort(leptons.begin(), leptons.end(), [](const Lepton& a, const Lepton &b){ return a.p4.Pt() > b.p4.Pt(); });

  // Store indices to leptons for each ID/Iso combination
  for(uint16_t idx = 0; idx < leptons.size(); idx++){
    for(const LepID::LepID& id: LepID::it){
      for(const LepIso::LepIso& iso: LepIso::it){
         
        uint16_t comb = LepIDIso(id, iso);
        const Lepton& m_lepton = leptons[idx];
        if(m_lepton.ID[id] && m_lepton.iso[iso])
          leptons_IDIso[comb].push_back(idx);

      }
    }
  }

  ///////////////////////////
  //       DILEPTONS       //
  ///////////////////////////

  #ifdef _TT_DEBUG_
    std::cout << "Dileptons" << std::endl;
  #endif

  for(uint16_t i1 = 0; i1 < leptons.size(); i1++){
    for(uint16_t i2 = i1 + 1; i2 < leptons.size(); i2++){
      const Lepton& l1 = leptons[i1];
      const Lepton& l2 = leptons[i2];

      DiLepton m_diLepton;

      m_diLepton.p4 = l1.p4 + l2.p4; 
      m_diLepton.idxs = std::make_pair(l1.idx, l2.idx); 
      m_diLepton.lidxs = std::make_pair(i1, i2); 
      m_diLepton.isElEl = l1.isEl && l2.isEl;
      m_diLepton.isElMu = l1.isEl && l2.isMu;
      m_diLepton.isMuEl = l1.isMu && l2.isEl;
      m_diLepton.isMuMu = l1.isMu && l2.isMu;
      m_diLepton.isOS = l1.charge != l2.charge;
      m_diLepton.isSF = m_diLepton.isElEl || m_diLepton.isMuMu;
 
      // Save the combination of IDs
      for(const LepID::LepID& id1: LepID::it){
        for(const LepID::LepID& id2: LepID::it){
          uint16_t idx = LepLepID(id1, id2);
          m_diLepton.ID[idx] = l1.ID[id1] && l2.ID[id2];
        }
      }

      // Save the combination of isolations
      for(const LepIso::LepIso& iso1: LepIso::it){
        for(const LepIso::LepIso& iso2: LepIso::it){
          uint16_t idx = LepLepIso(iso1, iso2);
          m_diLepton.iso[idx] = l1.iso[iso1] && l2.iso[iso2];
        }
      }
      
      m_diLepton.DR = VectorUtil::DeltaR(l1.p4, l2.p4);
      m_diLepton.DEta = TTWAnalysis::DeltaEta(l1.p4, l2.p4);
      m_diLepton.DPhi = VectorUtil::DeltaPhi(l1.p4, l2.p4);

      diLeptons.push_back(m_diLepton);
    }
  }

  // Save indices to DiLeptons for the combinations of IDs & Isolationss
  for(uint16_t i = 0; i < diLeptons.size(); i++){
    const DiLepton& m_diLepton = diLeptons[i];
    
    for(const LepID::LepID& id1: LepID::it){
      for(const LepID::LepID& id2: LepID::it){
        for(const LepIso::LepIso& iso1: LepIso::it){
          for(const LepIso::LepIso& iso2: LepIso::it){
            
            uint16_t idx_ids = LepLepID(id1, id2);
            uint16_t idx_isos = LepLepIso(iso1, iso2);
            uint16_t idx_comb = LepLepIDIso(id1, iso1, id2, iso2);

            if(m_diLepton.ID[idx_ids] && m_diLepton.iso[idx_isos])
              diLeptons_IDIso[idx_comb].push_back(i);
          
          }
        }
      }
    }
    
  }

  ///////////////////////////
  //       JETS            //
  ///////////////////////////

  #ifdef _TT_DEBUG_
    std::cout << "Jets" << std::endl;
  #endif

  const JetsProducer& jets = producers.get<JetsProducer>(m_jets_producer);

  // First find the jets passing kinematic cuts and save them as Jet objects

  uint16_t jetCounter(0);
  for(uint16_t ijet = 0; ijet < jets.p4.size(); ijet++){
    // Save the jets that pass the kinematic cuts
    if (std::abs(jets.p4[ijet].Eta()) < m_jetEtaCut && jets.p4[ijet].Pt() > m_jetPtCut){
      Jet m_jet;
      
      m_jet.p4 = jets.p4[ijet];
      m_jet.idx = ijet;
      m_jet.ID[JetID::L] = jets.passLooseID[ijet]; 
      m_jet.ID[JetID::T] = jets.passTightID[ijet]; 
      m_jet.ID[JetID::TLV] = jets.passTightLeptonVetoID[ijet];
      m_jet.CSVv2 = jets.getBTagDiscriminant(ijet, m_jetCSVv2Name);
      m_jet.BWP[BWP::L] = m_jet.CSVv2 > m_jetCSVv2L;
      m_jet.BWP[BWP::M] = m_jet.CSVv2 > m_jetCSVv2M;
      m_jet.BWP[BWP::T] = m_jet.CSVv2 > m_jetCSVv2T;
      
      // Save minimal DR(l,j) using selected leptons, for each Lepton ID/Iso
      for(const LepID::LepID& id: LepID::it){
        for(const LepIso::LepIso& iso: LepIso::it){
              
          uint16_t idx_comb = LepIDIso(id, iso);
          
          for(const uint16_t& lepIdx: leptons_IDIso[idx_comb]){
            const Lepton& m_lepton = leptons[lepIdx];
            float DR = (float) VectorUtil::DeltaR(jets.p4[ijet], m_lepton.p4);
            if( DR < m_jet.minDRjl_lepIDIso[idx_comb] )
              m_jet.minDRjl_lepIDIso[idx_comb] = DR;
          }
          
          // Save the indices to Jets passing the selected jetID and minDRjl > cut for this lepton ID/Iso
          if( m_jet.minDRjl_lepIDIso[idx_comb] > m_jetDRleptonCut && jetIDAccessor(jets, ijet, m_jetID) ){
            selJets_selID_DRCut[idx_comb].push_back(jetCounter);

            // Out of these, save the indices for different b-tagging working points
            for(const BWP::BWP& wp: BWP::it){
              uint16_t idx_comb_b = LepIDIsoJetBWP(id, iso, wp);
              if ((m_jet.BWP[wp]) && (std::abs(m_jet.p4.Eta()) < m_bJetEtaCut))
                selBJets_DRCut_BWP_PtOrdered[idx_comb_b].push_back(jetCounter);
            }
          }
        }
      }
      
      if(jetIDAccessor(jets, ijet, m_jetID)) // Save the indices to Jets passing the selected jet ID
        selJets_selID.push_back(jetCounter);
      
      selJets.push_back(m_jet);
      
      jetCounter++;
    }
  }

  // Sort the b-jets according to decreasing CSVv2 value
  selBJets_DRCut_BWP_CSVv2Ordered = selBJets_DRCut_BWP_PtOrdered;
  for(const LepID::LepID& id: LepID::it){
    for(const LepIso::LepIso& iso: LepIso::it){
      for(const BWP::BWP& wp: BWP::it){ 
        uint16_t idx_comb_b = LepIDIsoJetBWP(id, iso, wp);
        std::sort(selBJets_DRCut_BWP_CSVv2Ordered[idx_comb_b].begin(), selBJets_DRCut_BWP_CSVv2Ordered[idx_comb_b].end(), jetBTagDiscriminantSorter(jets, m_jetCSVv2Name, selJets));
      }
    }
  }
        
  ///////////////////////////
  //       DIJETS          //
  ///////////////////////////

  #ifdef _TT_DEBUG_
    std::cout << "Dijets" << std::endl;
  #endif

  // Next, construct DiJets out of selected jets with selected ID (not accounting for minDRjl here)

  uint16_t diJetCounter(0);

  for(uint16_t j1 = 0; j1 < selJets_selID.size(); j1++){
    for(uint16_t j2 = j1 + 1; j2 < selJets_selID.size(); j2++){
      const uint16_t jidx1 = selJets_selID[j1];
      const Jet& jet1 = selJets[jidx1];
      const uint16_t jidx2 = selJets_selID[j2];
      const Jet& jet2 = selJets[jidx2];

      DiJet m_diJet; 
      m_diJet.p4 = jet1.p4 + jet2.p4;
      m_diJet.idxs = std::make_pair(jet1.idx, jet2.idx);
      m_diJet.jidxs = std::make_pair(jidx1, jidx2);
      
      m_diJet.DR = VectorUtil::DeltaR(jet1.p4, jet2.p4);
      m_diJet.DEta = DeltaEta(jet1.p4, jet2.p4);
      m_diJet.DPhi = VectorUtil::DeltaPhi(jet1.p4, jet2.p4);
     
      for(const BWP::BWP& wp1: BWP::it){
        for(const BWP::BWP& wp2: BWP::it){
          uint16_t comb = JetJetBWP(wp1, wp2);
          m_diJet.BWP[comb] = jet1.BWP[wp1] && jet2.BWP[wp2];
        }
      }
      
      for(const LepID::LepID& id: LepID::it){
        for(const LepIso::LepIso& iso: LepIso::it){
          uint16_t combIDIso = LepIDIso(id, iso);

          m_diJet.minDRjl_lepIDIso[combIDIso] = std::min(jet1.minDRjl_lepIDIso[combIDIso], jet2.minDRjl_lepIDIso[combIDIso]);
          
          // Save the DiJets which have minDRjl>cut, for each leptonIDIso
          if(m_diJet.minDRjl_lepIDIso[combIDIso] > m_jetDRleptonCut){
            diJets_DRCut[combIDIso].push_back(diJetCounter);

            // Out of these, save di-b-jets for each combination of b-tagging working points
            for(const BWP::BWP& wp1: BWP::it){
              for(const BWP::BWP& wp2: BWP::it){
                uint16_t combB = JetJetBWP(wp1, wp2);
                uint16_t combAll = LepIDIsoJetJetBWP(id, iso, wp1, wp2);
                if ((m_diJet.BWP[combB])
                        && (std::abs(jet1.p4.Eta()) < m_bJetEtaCut)
                        && (std::abs(jet2.p4.Eta()) < m_bJetEtaCut))
                  diBJets_DRCut_BWP_PtOrdered[combAll].push_back(diJetCounter);
              }
            }
          
          }
        
        }
      }
      
      diJets.push_back(m_diJet); 
      diJetCounter++;
    }
  }

  // Order selected di-b-jets according to decreasing CSVv2 discriminant
  diBJets_DRCut_BWP_CSVv2Ordered = diBJets_DRCut_BWP_PtOrdered;
  for(const LepID::LepID& id: LepID::it){
    for(const LepIso::LepIso& iso: LepIso::it){
      for(const BWP::BWP& wp1: BWP::it){ 
        for(const BWP::BWP& wp2: BWP::it){ 
          uint16_t idx_comb_b = LepIDIsoJetJetBWP(id, iso, wp1, wp2);
          std::sort(diBJets_DRCut_BWP_CSVv2Ordered[idx_comb_b].begin(), diBJets_DRCut_BWP_CSVv2Ordered[idx_comb_b].end(), diJetBTagDiscriminantSorter(jets, m_jetCSVv2Name, diJets));
        }
      }
    }
  }
  
  ///////////////////////////
  //    EVENT VARIABLES    //
  ///////////////////////////
  
  #ifdef _TT_DEBUG_
    std::cout << "Dileptons-dijets" << std::endl;
  #endif

  // leptons-(b-)jets

  uint16_t diLepDiJetCounter(0);

  for(uint16_t dilep = 0; dilep < diLeptons.size(); dilep++){
    const DiLepton& m_diLepton = diLeptons[dilep];
    
    for(uint16_t dijet = 0; dijet < diJets.size(); dijet++){
      const DiJet& m_diJet =  diJets[dijet];
      
      DiLepDiJet m_diLepDiJet(m_diLepton, dilep, m_diJet, dijet);

      m_diLepDiJet.minDRjl = std::min( {
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.maxDRjl = std::max( {
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaR(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.minDEtajl = std::min( {
          DeltaEta(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          DeltaEta(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          DeltaEta(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          DeltaEta(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.maxDEtajl = std::max( {
          DeltaEta(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          DeltaEta(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          DeltaEta(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          DeltaEta(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.minDPhijl = std::min( {
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );
      m_diLepDiJet.maxDPhijl = std::max( {
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.first].p4, selJets[m_diJet.jidxs.second].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.first].p4),
          (float) VectorUtil::DeltaPhi(leptons[m_diLepton.lidxs.second].p4, selJets[m_diJet.jidxs.second].p4)
          } );

      diLepDiJets.push_back(m_diLepDiJet);

      for(const LepID::LepID& id1: LepID::it){
        for(const LepID::LepID& id2: LepID::it){
          for(const LepIso::LepIso& iso1: LepIso::it){
            for(const LepIso::LepIso& iso2: LepIso::it){
              
              uint16_t combID = LepLepID(id1, id2);
              LepID::LepID minID = std::min(id1, id2);
              
              uint16_t combIso = LepLepIso(iso1, iso2);
              LepIso::LepIso minIso = std::min(iso1, iso2);

              uint16_t minCombIDIso = LepIDIso(minID, minIso);
              uint16_t diLepCombIDIso = LepLepIDIso(id1, iso1, id2, iso2);
             
              // Store objects for each combined lepton ID/Iso, with jets having minDRjl>cut for leptons corresponding to the loosest combination of the aforementioned ID/Iso
              if(m_diLepton.ID[combID] && m_diLepton.iso[combIso] && m_diJet.minDRjl_lepIDIso[minCombIDIso] > m_jetDRleptonCut){
                diLepDiJets_DRCut[diLepCombIDIso].push_back(diLepDiJetCounter);
                
                // Out of these, store combinations of b-tagging working points
                for(const BWP::BWP& wp1: BWP::it){
                  for(const BWP::BWP& wp2: BWP::it){
                    uint16_t combB = JetJetBWP(wp1, wp2);
                    uint16_t combAll = LepLepIDIsoJetJetBWP(id1, iso1, id2, iso2, wp1, wp2);
                    if ((m_diJet.BWP[combB])
                            && (std::abs(jets.p4[m_diJet.idxs.first].Eta()) < m_bJetEtaCut)
                            && (std::abs(jets.p4[m_diJet.idxs.second].Eta()) < m_bJetEtaCut))
                      diLepDiBJets_DRCut_BWP_PtOrdered[combAll].push_back(diLepDiJetCounter);
                  }
                } // end b-jet loops

              } // end minDRjl>cut

            }
          } // end lepton iso loops
        }
      } // end lepton ID loops

      diLepDiJetCounter++;
    } // end dijet loop
  } // end dilepton loop

  // Order selected di-lepton-di-b-jets according to decreasing CSVv2 discriminant
  diLepDiBJets_DRCut_BWP_CSVv2Ordered = diLepDiBJets_DRCut_BWP_PtOrdered;
  
  for(const LepID::LepID& id1: LepID::it){
    for(const LepID::LepID& id2: LepID::it){
      
      for(const LepIso::LepIso& iso1: LepIso::it){
        for(const LepIso::LepIso& iso2: LepIso::it){
          
          for(const BWP::BWP& wp1: BWP::it){ 
            for(const BWP::BWP& wp2: BWP::it){ 
              
              uint16_t idx_comb_all = LepLepIDIsoJetJetBWP(id1, iso1, id2, iso2, wp1, wp2);
              std::sort(diLepDiBJets_DRCut_BWP_CSVv2Ordered[idx_comb_all].begin(), diLepDiBJets_DRCut_BWP_CSVv2Ordered[idx_comb_all].end(), diJetBTagDiscriminantSorter(jets, m_jetCSVv2Name, diLepDiJets));
            
            }
          }
        
        }
      }
    
    }
  }
      
  // leptons-(b-)jets-MET

  #ifdef _TT_DEBUG_
    std::cout << "Dileptons-Dijets-MET" << std::endl;
  #endif

  const METProducer &met = producers.get<METProducer>(m_met_producer);
  
  for(uint16_t i = 0; i < diLepDiJets.size(); i++){
    // Using regular MET
    DiLepDiJetMet m_diLepDiJetMet(diLepDiJets[i], i, met.p4);
    
    m_diLepDiJetMet.minDR_l_Met = std::min(
        (float) VectorUtil::DeltaR(leptons[m_diLepDiJetMet.diLepton->lidxs.first].p4, met.p4),
        (float) VectorUtil::DeltaR(leptons[m_diLepDiJetMet.diLepton->lidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.maxDR_l_Met = std::max(
        (float) VectorUtil::DeltaR(leptons[m_diLepDiJetMet.diLepton->lidxs.first].p4, met.p4),
        (float) VectorUtil::DeltaR(leptons[m_diLepDiJetMet.diLepton->lidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.minDEta_l_Met = std::min(
        DeltaEta(leptons[m_diLepDiJetMet.diLepton->lidxs.first].p4, met.p4),
        DeltaEta(leptons[m_diLepDiJetMet.diLepton->lidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.maxDEta_l_Met = std::max(
        DeltaEta(leptons[m_diLepDiJetMet.diLepton->lidxs.first].p4, met.p4),
        DeltaEta(leptons[m_diLepDiJetMet.diLepton->lidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.minDPhi_l_Met = std::min(
        (float) VectorUtil::DeltaPhi(leptons[m_diLepDiJetMet.diLepton->lidxs.first].p4, met.p4),
        (float) VectorUtil::DeltaPhi(leptons[m_diLepDiJetMet.diLepton->lidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.maxDPhi_l_Met = std::max(
        (float) VectorUtil::DeltaPhi(leptons[m_diLepDiJetMet.diLepton->lidxs.first].p4, met.p4),
        (float) VectorUtil::DeltaPhi(leptons[m_diLepDiJetMet.diLepton->lidxs.second].p4, met.p4)
        );

    m_diLepDiJetMet.minDR_j_Met = std::min(
        (float) VectorUtil::DeltaR(selJets[m_diLepDiJetMet.diJet->jidxs.first].p4, met.p4),
        (float) VectorUtil::DeltaR(selJets[m_diLepDiJetMet.diJet->jidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.maxDR_j_Met = std::max(
        (float) VectorUtil::DeltaR(selJets[m_diLepDiJetMet.diJet->jidxs.first].p4, met.p4),
        (float) VectorUtil::DeltaR(selJets[m_diLepDiJetMet.diJet->jidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.minDEta_j_Met = std::min(
        DeltaEta(selJets[m_diLepDiJetMet.diJet->jidxs.first].p4, met.p4),
        DeltaEta(selJets[m_diLepDiJetMet.diJet->jidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.maxDEta_j_Met = std::max(
        DeltaEta(selJets[m_diLepDiJetMet.diJet->jidxs.first].p4, met.p4),
        DeltaEta(selJets[m_diLepDiJetMet.diJet->jidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.minDPhi_j_Met = std::min(
        (float) VectorUtil::DeltaPhi(selJets[m_diLepDiJetMet.diJet->jidxs.first].p4, met.p4),
        (float) VectorUtil::DeltaPhi(selJets[m_diLepDiJetMet.diJet->jidxs.second].p4, met.p4)
        );
    m_diLepDiJetMet.maxDPhi_j_Met = std::max(
        (float) VectorUtil::DeltaPhi(selJets[m_diLepDiJetMet.diJet->jidxs.first].p4, met.p4),
        (float) VectorUtil::DeltaPhi(selJets[m_diLepDiJetMet.diJet->jidxs.second].p4, met.p4)
        );

    diLepDiJetsMet.push_back(m_diLepDiJetMet);

    for(const LepID::LepID& id1: LepID::it){
      for(const LepID::LepID& id2: LepID::it){
        for(const LepIso::LepIso& iso1: LepIso::it){
          for(const LepIso::LepIso& iso2: LepIso::it){
            
            uint16_t combID = LepLepID(id1, id2);
            LepID::LepID minID = std::min(id1, id2);
            
            uint16_t combIso = LepLepIso(iso1, iso2);
            LepIso::LepIso minIso = std::min(iso1, iso2);

            uint16_t minCombIDIso = LepIDIso(minID, minIso);
            uint16_t diLepCombIDIso = LepLepIDIso(id1, iso1, id2, iso2);
           
            // Store objects for each combined lepton ID/Iso, with jets having minDRjl>cut for leptons corresponding to the loosest combination of the aforementioned ID/Iso
            
            // First regular MET
            if(m_diLepDiJetMet.diLepton->ID[combID] && m_diLepDiJetMet.diLepton->iso[combIso] && m_diLepDiJetMet.diJet->minDRjl_lepIDIso[minCombIDIso] > m_jetDRleptonCut){
              diLepDiJetsMet_DRCut[diLepCombIDIso].push_back(i);
              
              // Out of these, store combinations of b-tagging working points
              for(const BWP::BWP& wp1: BWP::it){
                for(const BWP::BWP& wp2: BWP::it){
                  uint16_t combB = JetJetBWP(wp1, wp2);
                  uint16_t combAll = LepLepIDIsoJetJetBWP(id1, iso1, id2, iso2, wp1, wp2);
                  if ((m_diLepDiJetMet.diJet->BWP[combB])
                          && (std::abs(jets.p4[m_diLepDiJetMet.diJet->idxs.first].Eta()) < m_bJetEtaCut)
                          && (std::abs(jets.p4[m_diLepDiJetMet.diJet->idxs.second].Eta()) < m_bJetEtaCut))
                    diLepDiBJetsMet_DRCut_BWP_PtOrdered[combAll].push_back(i);
                }
              } // end b-jet loops

            } // end minDRjl>cut

          }
        } // end lepton iso loops
      }
    } // end lepton ID loops
     
  } // end diLepDiJet loop
  
  // Store objects according to CSVv2
  // First regular MET
  diLepDiBJetsMet_DRCut_BWP_CSVv2Ordered = diLepDiBJetsMet_DRCut_BWP_PtOrdered; 
  for(const LepID::LepID& id1: LepID::it){
    for(const LepID::LepID& id2: LepID::it){
      
      for(const LepIso::LepIso& iso1: LepIso::it){
        for(const LepIso::LepIso& iso2: LepIso::it){
          
          for(const BWP::BWP& wp1: BWP::it){ 
            for(const BWP::BWP& wp2: BWP::it){ 
              
              uint16_t idx_comb_all = LepLepIDIsoJetJetBWP(id1, iso1, id2, iso2, wp1, wp2);
              std::sort(diLepDiBJetsMet_DRCut_BWP_CSVv2Ordered[idx_comb_all].begin(), diLepDiBJetsMet_DRCut_BWP_CSVv2Ordered[idx_comb_all].end(), diJetBTagDiscriminantSorter(jets, m_jetCSVv2Name, diLepDiJetsMet));
            
            }
          }
        
        }
      }
    
    }
  }
  
  ///////////////////////////
  //       TRIGGER         //
  ///////////////////////////

  #ifdef _TT_DEBUG_
    std::cout << "Trigger" << std::endl;
  #endif

  if (producers.exists("hlt")) {

      const HLTProducer& hlt = producers.get<HLTProducer>("hlt");

      if (hlt.paths.empty()) {
#if TT_HLT_DEBUG
          std::cout << "No HLT path triggered for this event. Skipping HLT matching." << std::endl;
#endif
          goto after_hlt_matching;
      }

#if TT_HLT_DEBUG
      std::cout << "HLT path triggered for this event:" << std::endl;
      for (const std::string& path: hlt.paths) {
          std::cout << "\t" << path << std::endl;
      }
#endif

      /*
       * Try to match `lepton` with an online object, using a deltaR and a deltaPt cut
       * Returns the index inside the HLTProducer collection, or -1 if no match is found.
       */
      auto matchOfflineLepton = [&](Lepton& lepton) {

          if (lepton.hlt_already_tried_matching)
              return lepton.hlt_idx;

#if TT_HLT_DEBUG
          std::cout << "Trying to match offline lepton: " << std::endl;
          std::cout << "\tMuon? " << lepton.isMu << " ; Pt: " << lepton.p4.Pt() << " ; Eta: " << lepton.p4.Eta() << " ; Phi: " << lepton.p4.Phi() << " ; E: " << lepton.p4.E() << std::endl;
#endif

          float min_dr = std::numeric_limits<float>::max();

          int16_t index = -1;
          for (size_t hlt_object = 0; hlt_object < hlt.object_p4.size(); hlt_object++) {

              float dr = VectorUtil::DeltaR(lepton.p4, hlt.object_p4[hlt_object]);
              float dpt_over_pt = std::abs(lepton.p4.Pt() - hlt.object_p4[hlt_object].Pt()) / lepton.p4.Pt();

              if (dr > m_hltDRCut)
                  continue;

              if (dpt_over_pt > m_hltDPtCut)
                  continue;

              if (dr < min_dr) {
                  min_dr = dr;
                  index = hlt_object;
              }
          }

#if TT_HLT_DEBUG
          if (index != -1) {
              std::cout << "\033[32mMatched with online object:\033[00m" << std::endl;
              std::cout << "\tPDG Id: " << hlt.object_pdg_id[index] << " ; Pt: " << hlt.object_p4[index].Pt() << " ; Eta: " << hlt.object_p4[index].Eta() << " ; Phi: " << hlt.object_p4[index].Phi() << " ; E: " << hlt.object_p4[index].E() << std::endl;
              std::cout << "\tΔR: " << min_dr << " ; ΔPt / Pt: " << std::abs(lepton.p4.Pt() - hlt.object_p4[index].Pt()) / lepton.p4.Pt() << std::endl;
          } else {
              std::cout << "\033[31mNo match found\033[00m" << std::endl;
          }
#endif

          lepton.hlt_idx = index;
          lepton.hlt_already_tried_matching = true;
          lepton.hlt_DR_matched_object = min_dr;
          lepton.hlt_DPt_matched_object = std::abs(lepton.p4.Pt() - hlt.object_p4[index].Pt()) / lepton.p4.Pt();

          return index;
      };

      // Iterate over all dilepton pairs
      for (auto& m_diLepton: diLeptons) {
          // For each lepton of this pair, find the online object
          m_diLepton.hlt_idxs = std::make_pair(
                  matchOfflineLepton(leptons[m_diLepton.lidxs.first]),
                  matchOfflineLepton(leptons[m_diLepton.lidxs.second])
         );
      }

  }

after_hlt_matching:

    ///////////////////////////
    //       GEN INFO        //
    ///////////////////////////

    #ifdef _TT_DEBUG_
      std::cout << "Generator" << std::endl;
    #endif

    if (event.isRealData())
        return;

    const GenParticlesProducer& gen_particles = producers.get<GenParticlesProducer>("gen_particles");

    // 'Pruned' particles are from the hard process
    // 'Packed' particles are stable particles

    std::function<bool(size_t, size_t)> pruned_decays_from = [&pruned_decays_from, &gen_particles](size_t particle_index, size_t mother_index) -> bool {
        // Iterator over all pruned particles to find if the particle `particle_index` has `mother_index` in its decay history
        if (gen_particles.pruned_mothers_index[particle_index].empty())
            return false;

        size_t index = gen_particles.pruned_mothers_index[particle_index][0];

        if (index == mother_index) {
            return true;
        }

        if (pruned_decays_from(index, mother_index))
            return true;

        return false;
    };

#if TT_GEN_DEBUG
    std::function<void(size_t)> print_mother_chain = [&gen_particles, &print_mother_chain](size_t p) {

        if (gen_particles.pruned_mothers_index[p].empty()) {
            std::cout << std::endl;
            return;
        }

        size_t index = gen_particles.pruned_mothers_index[p][0];
            std::cout << " <- #" << index << "(" << gen_particles.pruned_pdg_id[index] << ")";
            print_mother_chain(index);
    };
#endif

    // We need to initialize everything to -1, since 0 is a valid entry in the tt_genParticles array
    gen_t = -1;
    gen_t_beforeFSR = -1;
    gen_tbar = -1;
    gen_tbar_beforeFSR = -1;
    
    gen_b = -1;
    gen_b_beforeFSR = -1;
    gen_bbar = -1;
    gen_bbar_beforeFSR = -1;
    
    gen_jet1_t = -1;
    gen_jet1_t_beforeFSR = -1;
    gen_jet1_tbar = -1;
    gen_jet1_tbar_beforeFSR = -1;
    gen_jet2_t = -1;
    gen_jet2_t_beforeFSR = -1;
    gen_jet2_tbar = -1;
    gen_jet2_tbar_beforeFSR = -1;
    
    gen_lepton_t = -1;
    gen_lepton_t_beforeFSR = -1;
    gen_neutrino_t = -1;
    gen_neutrino_t_beforeFSR = -1;
    gen_lepton_tbar = -1;
    gen_lepton_tbar_beforeFSR = -1;
    gen_neutrino_tbar = -1;
    gen_neutrino_tbar_beforeFSR = -1;
    
    gen_matched_lepton_t = -1;
    gen_matched_lepton_tbar = -1;

    size_t gen_index(0);
    for (size_t i = 0; i < gen_particles.pruned_pdg_id.size(); i++) {

        int16_t pdg_id = gen_particles.pruned_pdg_id[i];
        uint16_t a_pdg_id = std::abs(pdg_id);

        // We only care of particles with PDG id <= 16 (16 is neutrino tau)
        if (a_pdg_id > 16)
            continue;

        GenStatusFlags flags(gen_particles.pruned_status_flags[i]);

        if (! flags.isLastCopy() && ! flags.isFirstCopy())
            continue;

        if (! flags.fromHardProcess())
            continue;

#if TT_GEN_DEBUG
        std::cout << "---" << std::endl;
        std::cout << "Gen particle #" << i << ": PDG id: " << gen_particles.pruned_pdg_id[i];
        print_mother_chain(i);
        flags.dump();
#endif

        if (pdg_id == 6) {
            FILL_GEN_COLL(t);
            continue;
        } else if (pdg_id == -6) {
            FILL_GEN_COLL(tbar);
            continue;
        }

        if (gen_t == -1 || gen_tbar == -1) {
            // Don't bother if we don't have found the tops
            continue;
        }

        bool from_t_decay = pruned_decays_from(i, genParticles[gen_t].pruned_idx);
        bool from_tbar_decay = pruned_decays_from(i, genParticles[gen_tbar].pruned_idx);

        // Only keep particles coming from the tops decay
        if (! from_t_decay && ! from_tbar_decay)
            continue;

        if (pdg_id == 5) {
            // Maybe it's a b coming from the W decay
            if (!flags.isFirstCopy() && flags.isLastCopy() && gen_b == -1) {

                // This can be a B decaying from a W
                // However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
                // Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
                // If yes, then it's not the B coming directly from the top decay
                if ((gen_jet1_t_beforeFSR != -1 && std::abs(genParticles[gen_jet1_t_beforeFSR].pdg_id) == 5) ||
                    (gen_jet2_t_beforeFSR != -1 && std::abs(genParticles[gen_jet2_t_beforeFSR].pdg_id) == 5) ||
                    (gen_jet1_tbar_beforeFSR != -1 && std::abs(genParticles[gen_jet1_tbar_beforeFSR].pdg_id) == 5) ||
                    (gen_jet2_tbar_beforeFSR != -1 && std::abs(genParticles[gen_jet2_tbar_beforeFSR].pdg_id) == 5)) {

#if TT_GEN_DEBUG
                    std::cout << "A quark coming from W decay is a b" << std::endl;
#endif

                    if (! (gen_jet1_tbar_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet1_tbar_beforeFSR].pruned_idx)) &&
                        ! (gen_jet2_tbar_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet2_tbar_beforeFSR].pruned_idx)) &&
                        ! (gen_jet1_t_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet1_t_beforeFSR].pruned_idx)) &&
                        ! (gen_jet2_t_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet2_t_beforeFSR].pruned_idx))) {
#if TT_GEN_DEBUG
                        std::cout << "This after-FSR b quark is not coming from a W decay" << std::endl;
#endif
                        FILL_GEN_COLL(b);
                        continue;
                    }
#if TT_GEN_DEBUG
                    else {
                        std::cout << "This after-FSR b quark comes from a W decay" << std::endl;
                    }
#endif
                } else {
#if TT_GEN_DEBUG
                    std::cout << "Assigning gen_b" << std::endl;
#endif
                    FILL_GEN_COLL(b);
                    continue;
                }
            } else if (flags.isFirstCopy() && gen_b_beforeFSR == -1) {
                FILL_GEN_COLL(b);
                continue;
            } else {
#if TT_GEN_DEBUG
                std::cout << "This should not happen!" << std::endl;
#endif
            }
        } else if (pdg_id == -5) {
            if (!flags.isFirstCopy() && flags.isLastCopy() && gen_bbar == -1) {

                // This can be a B decaying from a W
                // However, we can't rely on the presence of the W in the decay chain, as it may be generator specific
                // Since it's the last copy (ie, after FSR), we can check if this B comes from the B assigned to the W decay (ie, gen_jet1_t_beforeFSR, gen_jet2_t_beforeFSR)
                // If yes, then it's not the B coming directly from the top decay
                if ((gen_jet1_t_beforeFSR != -1 && std::abs(genParticles[gen_jet1_t_beforeFSR].pdg_id) == 5) ||
                    (gen_jet2_t_beforeFSR != -1 && std::abs(genParticles[gen_jet2_t_beforeFSR].pdg_id) == 5) ||
                    (gen_jet1_tbar_beforeFSR != -1 && std::abs(genParticles[gen_jet1_tbar_beforeFSR].pdg_id) == 5) ||
                    (gen_jet2_tbar_beforeFSR != -1 && std::abs(genParticles[gen_jet2_tbar_beforeFSR].pdg_id) == 5)) {

#if TT_GEN_DEBUG
                    std::cout << "A quark coming from W decay is a bbar" << std::endl;
#endif

                    if (! (gen_jet1_tbar_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet1_tbar_beforeFSR].pruned_idx)) &&
                        ! (gen_jet2_tbar_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet2_tbar_beforeFSR].pruned_idx)) &&
                        ! (gen_jet1_t_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet1_t_beforeFSR].pruned_idx)) &&
                        ! (gen_jet2_t_beforeFSR != -1 && pruned_decays_from(i, genParticles[gen_jet2_t_beforeFSR].pruned_idx))) {
#if TT_GEN_DEBUG
                        std::cout << "This after-fsr b anti-quark is not coming from a W decay" << std::endl;
#endif
                        FILL_GEN_COLL(bbar);
                        continue;
                    }
#if TT_GEN_DEBUG
                    else {
                        std::cout << "This after-fsr b anti-quark comes from a W decay" << std::endl;
                    }
#endif
                } else {
#if TT_GEN_DEBUG
                    std::cout << "Assigning gen_bbar" << std::endl;
#endif
                    FILL_GEN_COLL(bbar);
                    continue;
                }
            } else if (flags.isFirstCopy() && gen_bbar_beforeFSR == -1) {
                FILL_GEN_COLL(bbar);
                continue;
            }
        }

        if ((gen_tbar == -1) || (gen_t == -1))
            continue;

        if (gen_t != -1 && from_t_decay) {
#if TT_GEN_DEBUG
        std::cout << "Coming from the top chain decay" << std::endl;
#endif
            if (a_pdg_id >= 1 && a_pdg_id <= 5) {
                FILL_GEN_COLL2(jet1_t, jet2_t, "Error: more than two quarks coming from top decay");
            } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
                FILL_GEN_COLL(lepton_t);
            } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
                FILL_GEN_COLL(neutrino_t);
            } else {
                std::cout << "Error: unknown particle coming from top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
            }
        } else if (gen_tbar != -1 && from_tbar_decay) {
#if TT_GEN_DEBUG
        std::cout << "Coming from the anti-top chain decay" << std::endl;
#endif
            if (a_pdg_id >= 1 && a_pdg_id <= 5) {
                FILL_GEN_COLL2(jet1_tbar, jet2_tbar, "Error: more than two quarks coming from anti-top decay");
            } else if (a_pdg_id == 11 || a_pdg_id == 13 || a_pdg_id == 15) {
                FILL_GEN_COLL(lepton_tbar);
            } else if (a_pdg_id == 12 || a_pdg_id == 14 || a_pdg_id == 16) {
                FILL_GEN_COLL(neutrino_tbar);
            } else {
                std::cout << "Error: unknown particle coming from anti-top decay - #" << i << " ; PDG Id: " << pdg_id << std::endl;
            }
        }
    }

    if (!gen_t || !gen_tbar) {
#if TT_GEN_DEBUG
        std::cout << "This is not a ttbar event" << std::endl;
#endif
        gen_ttbar_decay_type = NotTT;
        return;
    }

    if ((gen_jet1_t != -1) && (gen_jet2_t != -1) && (gen_jet1_tbar != -1) && (gen_jet2_tbar != -1)) {
#if TT_GEN_DEBUG
        std::cout << "Hadronic ttbar decay" << std::endl;
#endif
        gen_ttbar_decay_type = Hadronic;
    } else if (
            ((gen_lepton_t != -1) && (gen_lepton_tbar == -1)) ||
            ((gen_lepton_t == -1) && (gen_lepton_tbar != -1))
            ) {

#if TT_GEN_DEBUG
        std::cout << "Semileptonic ttbar decay" << std::endl;
#endif

        uint16_t lepton_pdg_id;
        if (gen_lepton_t != -1)
            lepton_pdg_id = std::abs(genParticles[gen_lepton_t].pdg_id);
        else
            lepton_pdg_id = std::abs(genParticles[gen_lepton_tbar].pdg_id);

        if (lepton_pdg_id == 11)
            gen_ttbar_decay_type = Semileptonic_e;
        else if (lepton_pdg_id == 13)
            gen_ttbar_decay_type = Semileptonic_mu;
        else
            gen_ttbar_decay_type = Semileptonic_tau;
    } else if (gen_lepton_t != -1 && gen_lepton_tbar != -1) {
        uint16_t lepton_t_pdg_id = std::abs(genParticles[gen_lepton_t].pdg_id);
        uint16_t lepton_tbar_pdg_id = std::abs(genParticles[gen_lepton_tbar].pdg_id);

#if TT_GEN_DEBUG
        std::cout << "Dileptonic ttbar decay" << std::endl;
#endif

        if (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 11)
            gen_ttbar_decay_type = Dileptonic_ee;
        else if (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 13)
            gen_ttbar_decay_type = Dileptonic_mumu;
        else if (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 15)
            gen_ttbar_decay_type = Dileptonic_tautau;
        else if (
                (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 13) ||
                (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 11)
                ) {
            gen_ttbar_decay_type = Dileptonic_mue;
        }
        else if (
                (lepton_t_pdg_id == 11 && lepton_tbar_pdg_id == 15) ||
                (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 11)
                ) {
            gen_ttbar_decay_type = Dileptonic_etau;
        }
        else if (
                (lepton_t_pdg_id == 13 && lepton_tbar_pdg_id == 15) ||
                (lepton_t_pdg_id == 15 && lepton_tbar_pdg_id == 13)
                ) {
            gen_ttbar_decay_type = Dileptonic_mutau;
        } else {
            std::cout << "Error: unknown dileptonic ttbar decay." << std::endl;
            gen_ttbar_decay_type = NotTT;
            return;
        }
    } else {
        std::cout << "Error: unknown ttbar decay." << std::endl;
        gen_ttbar_decay_type = UnknownTT;
    }

    gen_ttbar_p4 = genParticles[gen_t].p4 + genParticles[gen_tbar].p4;
    if (gen_t_beforeFSR != -1 && gen_tbar_beforeFSR != -1)
        gen_ttbar_beforeFSR_p4 = genParticles[gen_t_beforeFSR].p4 + genParticles[gen_tbar_beforeFSR].p4;

    gen_t_tbar_deltaR = VectorUtil::DeltaR(genParticles[gen_t].p4, genParticles[gen_tbar].p4);
    gen_t_tbar_deltaEta = DeltaEta(genParticles[gen_t].p4, genParticles[gen_tbar].p4);
    gen_t_tbar_deltaPhi = VectorUtil::DeltaPhi(genParticles[gen_t].p4, genParticles[gen_tbar].p4);

    gen_b_bbar_deltaR = VectorUtil::DeltaR(genParticles[gen_b].p4, genParticles[gen_bbar].p4);

    if (gen_ttbar_decay_type > Hadronic) {

        float min_dr_lepton_t = std::numeric_limits<float>::max();
        float min_dr_lepton_tbar = std::numeric_limits<float>::max();

        size_t lepton_index = 0;
        for (auto& lepton: leptons) {

            if (gen_lepton_t != -1) {
                float dr = VectorUtil::DeltaR(genParticles[gen_lepton_t].p4, lepton.p4);
                if (dr < min_dr_lepton_t &&
                        (std::abs(lepton.pdg_id()) == std::abs(genParticles[gen_lepton_t].pdg_id))) {
                    min_dr_lepton_t = dr;
                    gen_matched_lepton_t = lepton_index;
                }
                gen_lepton_t_deltaR.push_back(dr);
            }

            if (gen_lepton_tbar != -1) {
                float dr = VectorUtil::DeltaR(genParticles[gen_lepton_tbar].p4, lepton.p4);
                if (dr < min_dr_lepton_tbar &&
                        (std::abs(lepton.pdg_id()) == std::abs(genParticles[gen_lepton_tbar].pdg_id))) {
                    min_dr_lepton_tbar = dr;
                    gen_matched_lepton_tbar = lepton_index;
                }
                gen_lepton_tbar_deltaR.push_back(dr);
            }

            lepton_index++;
        }
    }

    // Match b quarks to jets

    const float MIN_DR_JETS = 0.8;
    for (const auto& id: LepID::it) {
      for (const auto& iso: LepIso::it) {
          uint16_t IdWP = LepIDIso(id, iso);

          float min_dr_b = MIN_DR_JETS;
          float min_dr_bbar = MIN_DR_JETS;
          float min_dr_b_beforeFSR = MIN_DR_JETS;
          float min_dr_bbar_beforeFSR = MIN_DR_JETS;
          size_t jet_index = 0;

          int16_t local_gen_matched_b = -1;
          int16_t local_gen_matched_bbar = -1;
          int16_t local_gen_matched_b_beforeFSR = -1;
          int16_t local_gen_matched_bbar_beforeFSR = -1;
          for (auto& jet: selJets_selID_DRCut[IdWP]) {
              float dr = VectorUtil::DeltaR(genParticles[gen_b].p4, jets.p4[jet]);
              if (dr < min_dr_b) {
                  min_dr_b = dr;
                  local_gen_matched_b = jet_index;
              }
              gen_b_deltaR[IdWP].push_back(dr);

              dr = VectorUtil::DeltaR(genParticles[gen_b_beforeFSR].p4, jets.p4[jet]);
              if (dr < min_dr_b_beforeFSR) {
                  min_dr_b_beforeFSR = dr;
                  local_gen_matched_b_beforeFSR = jet_index;
              }
              gen_b_beforeFSR_deltaR[IdWP].push_back(dr);

              dr = VectorUtil::DeltaR(genParticles[gen_bbar].p4, jets.p4[jet]);
              if (dr < min_dr_bbar) {
                  min_dr_bbar = dr;
                  local_gen_matched_bbar = jet_index;
              }
              gen_bbar_deltaR[IdWP].push_back(dr);

              dr = VectorUtil::DeltaR(genParticles[gen_bbar_beforeFSR].p4, jets.p4[jet]);
              if (dr < min_dr_bbar_beforeFSR) {
                  min_dr_bbar_beforeFSR = dr;
                  local_gen_matched_bbar_beforeFSR = jet_index;
              }
              gen_bbar_beforeFSR_deltaR[IdWP].push_back(dr);

              jet_index++;
          }

          gen_matched_b[IdWP] = local_gen_matched_b;
          gen_matched_bbar[IdWP] = local_gen_matched_bbar;
          gen_matched_b_beforeFSR[IdWP] = local_gen_matched_b_beforeFSR;
          gen_matched_bbar_beforeFSR[IdWP] = local_gen_matched_bbar_beforeFSR;
      }
    }

    if (gen_b > -1 && gen_lepton_t > -1) {
        gen_b_lepton_t_deltaR = VectorUtil::DeltaR(genParticles[gen_b].p4, genParticles[gen_lepton_t].p4);
    }

    if (gen_bbar > -1 && gen_lepton_tbar > -1) {
        gen_bbar_lepton_tbar_deltaR = VectorUtil::DeltaR(genParticles[gen_bbar].p4, genParticles[gen_lepton_tbar].p4);
    }

    #ifdef _TT_DEBUG_
      std::cout << "End event." << std::endl;
    #endif

}

void TTAnalyzer::registerCategories(CategoryManager& manager, const edm::ParameterSet& config) {
  manager.new_category<TTWAnalysis::ElElCategory>("elel", "Category with leading leptons as two electrons", config);
  manager.new_category<TTWAnalysis::ElMuCategory>("elmu", "Category with leading leptons as electron, muon", config);
  manager.new_category<TTWAnalysis::MuElCategory>("muel", "Category with leading leptons as muon, electron", config);
  manager.new_category<TTWAnalysis::MuMuCategory>("mumu", "Category with leading leptons as two muons", config);
}

#include <FWCore/PluginManager/interface/PluginFactory.h>
DEFINE_EDM_PLUGIN(ExTreeMakerAnalyzerFactory, TTAnalyzer, "tt_analyzer");
