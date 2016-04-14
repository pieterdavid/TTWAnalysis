#pragma once

#include <string>
#include <utility>
#include <vector>
#include <limits>

#include <cp3_llbb/Framework/interface/MuonsProducer.h>
#include <cp3_llbb/Framework/interface/JetsProducer.h>
#include <cp3_llbb/Framework/interface/Analyzer.h>

#include <cp3_llbb/TTAnalysis/interface/Types.h>
#include <cp3_llbb/TTAnalysis/interface/Tools.h>

class TTAnalyzer: public Framework::Analyzer {
    public:
        TTAnalyzer(const std::string& name, const ROOT::TreeGroup& tree_, const edm::ParameterSet& config):
            Analyzer(name, tree_, config),

            // Not untracked as these parameters are mandatory
            m_electrons_producer(config.getParameter<std::string>("electronsProducer")),
            m_muons_producer(config.getParameter<std::string>("muonsProducer")),
            m_jets_producer(config.getParameter<std::string>("jetsProducer")),
            m_met_producer(config.getParameter<std::string>("metProducer")),
            
            m_electronPtCut( config.getUntrackedParameter<double>("electronPtCut", 20) ),
            m_electronEtaCut( config.getUntrackedParameter<double>("electronEtaCut", 2.5) ),
            m_electronVetoIDName( config.getUntrackedParameter<std::string>("electronVetoIDName") ),
            m_electronLooseIDName( config.getUntrackedParameter<std::string>("electronLooseIDName") ),
            m_electronMediumIDName( config.getUntrackedParameter<std::string>("electronMediumIDName") ),
            m_electronTightIDName( config.getUntrackedParameter<std::string>("electronTightIDName") ),
            
            m_muonPtCut( config.getUntrackedParameter<double>("muonPtCut", 20) ),
            m_muonEtaCut( config.getUntrackedParameter<double>("muonEtaCut", 2.4) ),
            m_muonLooseIsoCut( config.getUntrackedParameter<double>("muonLooseIsoCut", 0.2) ),
            m_muonTightIsoCut( config.getUntrackedParameter<double>("muonTightIsoCut", 0.12) ),
            
            m_jetPtCut( config.getUntrackedParameter<double>("jetPtCut", 30) ),
            m_jetEtaCut( config.getUntrackedParameter<double>("jetEtaCut", 2.5) ),
            m_bJetEtaCut( config.getUntrackedParameter<double>("bJetEtaCut", 2.4) ),
            m_jetPUID( config.getUntrackedParameter<double>("jetPUID", std::numeric_limits<float>::min()) ),
            m_jetDRleptonCut( config.getUntrackedParameter<double>("jetDRleptonCut", 0.3) ),
            m_jetID( config.getUntrackedParameter<std::string>("jetID", "loose") ),
            m_jetCSVv2Name( config.getUntrackedParameter<std::string>("jetCSVv2Name", "pfCombinedInclusiveSecondaryVertexV2BJetTags") ),
            m_jetCSVv2L( config.getUntrackedParameter<double>("jetCSVv2L", 0.605) ),
            m_jetCSVv2M( config.getUntrackedParameter<double>("jetCSVv2M", 0.89) ),
            m_jetCSVv2T( config.getUntrackedParameter<double>("jetCSVv2T", 0.97) ),
            
            m_hltDRCut( config.getUntrackedParameter<double>("hltDRCut", std::numeric_limits<float>::max()) ),
            m_hltDPtCut( config.getUntrackedParameter<double>("hltDPtCut", std::numeric_limits<float>::max()) )
        {
        }

        virtual void analyze(const edm::Event&, const edm::EventSetup&, const ProducersManager&, const AnalyzersManager&, const CategoryManager&) override;
        virtual void registerCategories(CategoryManager& manager, const edm::ParameterSet&) override;

        BRANCH(electrons_IDIso, std::vector<std::vector<uint16_t>>);
        BRANCH(muons_IDIso, std::vector<std::vector<uint16_t>>);

        BRANCH(leptons, std::vector<TTAnalysis::Lepton>);
        BRANCH(leptons_IDIso, std::vector<std::vector<uint16_t>>);

        BRANCH(diLeptons, std::vector<TTAnalysis::DiLepton>);
        BRANCH(diLeptons_IDIso, std::vector<std::vector<uint16_t>>);

        BRANCH(selJets, std::vector<TTAnalysis::Jet>);
        BRANCH(selJets_selID, std::vector<uint16_t>);
        // ex.: selectedJets_..._DRCut[X][0] is the highest Pt selected jet with minDRjl>0.3 taking into account ID/Iso-X Leptons
        BRANCH(selJets_selID_DRCut, std::vector<std::vector<uint16_t>>);
        // ex.: selectedBJets_..._PtOrdered[X][0] is the highest Pt selected jet with minDRjl>0.3 taking into account ID/Iso/Btag-X combination
        BRANCH(selBJets_DRCut_BWP_PtOrdered, std::vector<std::vector<uint16_t>>);
        BRANCH(selBJets_DRCut_BWP_CSVv2Ordered, std::vector<std::vector<uint16_t>>);

        BRANCH(diJets, std::vector<TTAnalysis::DiJet>);
        // ex.: diJets_DRCut[X][0] is first diJet with minDRjl>0.3 taking into account ID/Iso-X Leptons
        BRANCH(diJets_DRCut, std::vector<std::vector<uint16_t>>); 
        // ex.: diBJets_..._CSVv2Ordered[X][0] is the b-jet pair with highest CSVv2 values and with minDRjl>0.3 taking into account the leptonID/Iso/Btag-X combination
        BRANCH(diBJets_DRCut_BWP_PtOrdered, std::vector<std::vector<uint16_t>>);
        BRANCH(diBJets_DRCut_BWP_CSVv2Ordered, std::vector<std::vector<uint16_t>>);

        // For all the following: indices are combinations of LeptonID/LeptonIso/(B-tagging working point)

        BRANCH(diLepDiJets, std::vector<TTAnalysis::DiLepDiJet>);
        
        BRANCH(diLepDiJets_DRCut, std::vector<std::vector<uint16_t>>); // di-leptons of combined ID/Iso with di-jets built out of jets having minDRjl>cut taking into account lepton ID/Iso corresponding to the loosest combination of the two leptons of the object
        BRANCH(diLepDiBJets_DRCut_BWP_PtOrdered, std::vector<std::vector<uint16_t>>);
        BRANCH(diLepDiBJets_DRCut_BWP_CSVv2Ordered, std::vector<std::vector<uint16_t>>);

        BRANCH(diLepDiJetsMet, std::vector<TTAnalysis::DiLepDiJetMet>);
        
        BRANCH(diLepDiJetsMet_DRCut, std::vector<std::vector<uint16_t>>); 
        BRANCH(diLepDiBJetsMet_DRCut_BWP_PtOrdered, std::vector<std::vector<uint16_t>>);
        BRANCH(diLepDiBJetsMet_DRCut_BWP_CSVv2Ordered, std::vector<std::vector<uint16_t>>);

        // Gen matching. All indexes are from the `genParticles` collection
        BRANCH(genParticles, std::vector<TTAnalysis::GenParticle>);
        BRANCH(gen_t, int16_t); // Index of the top quark
        BRANCH(gen_t_beforeFSR, int16_t); // Index of the top quark, before any FSR
        BRANCH(gen_tbar, int16_t); // Index of the anti-top quark
        BRANCH(gen_tbar_beforeFSR, int16_t); // Index of the anti-top quark, before any FSR
        BRANCH(gen_t_tbar_deltaR, float); // DeltaR between the top and the anti-top quark
        BRANCH(gen_t_tbar_deltaEta, float); // DeltaEta between the top and the anti-top quark
        BRANCH(gen_t_tbar_deltaPhi, float); // DeltaPhi between the top and the anti-top quark

        BRANCH(gen_b, int16_t); // Index of the b quark coming from the top decay
        BRANCH(gen_b_beforeFSR, int16_t); // Index of the b quark coming from the top decay, before any FSR
        BRANCH(gen_bbar, int16_t); // Index of the anti-b quark coming from the anti-top decay
        BRANCH(gen_bbar_beforeFSR, int16_t); // Index of the anti-b quark coming from the anti-top decay, before any FSR
        BRANCH(gen_b_bbar_deltaR, float); // DeltaR between the b and the anti-b quark

        BRANCH(gen_jet1_t, int16_t); // Index of the first jet from the top decay chain
        BRANCH(gen_jet1_t_beforeFSR, int16_t); // Index of the first jet from the top decay chain, before any FSR
        BRANCH(gen_jet2_t, int16_t); // Index of the second jet from the top decay chain
        BRANCH(gen_jet2_t_beforeFSR, int16_t); // Index of the second jet from the top decay chain, before any FSR

        BRANCH(gen_jet1_tbar, int16_t); // Index of the first jet from the anti-top decay chain
        BRANCH(gen_jet1_tbar_beforeFSR, int16_t); // Index of the first jet from the anti-top decay chain, before any FSR
        BRANCH(gen_jet2_tbar, int16_t); // Index of the second jet from the anti-top decay chain
        BRANCH(gen_jet2_tbar_beforeFSR, int16_t); // Index of the second jet from the anti-top decay chain, before any FSR

        BRANCH(gen_lepton_t, int16_t); // Index of the lepton from the top decay chain
        BRANCH(gen_lepton_t_beforeFSR, int16_t); // Index of the lepton from the top decay chain, before any FSR
        BRANCH(gen_neutrino_t, int16_t); // Index of the neutrino from the top decay chain
        BRANCH(gen_neutrino_t_beforeFSR, int16_t); // Index of the neutrino from the top decay chain, before any FSR

        BRANCH(gen_lepton_tbar, int16_t); // Index of the lepton from the anti-top decay chain
        BRANCH(gen_lepton_tbar_beforeFSR, int16_t); // Index of the lepton from the anti-top decay chain, before any FSR
        BRANCH(gen_neutrino_tbar, int16_t); // Index of the neutrino from the anti-top decay chain
        BRANCH(gen_neutrino_tbar_beforeFSR, int16_t); // Index of the neutrino from the anti-top decay chain, before any FSR

        BRANCH(gen_ttbar_decay_type, char); // Type of ttbar decay. Can take any values from TTDecayType enum

        BRANCH(gen_ttbar_beforeFSR_p4, LorentzVector);
        BRANCH(gen_ttbar_p4, LorentzVector);

        // Matching for the dileptonic case

        BRANCH(gen_b_lepton_t_deltaR, float); // DeltaR between the b quark and the lepton coming from the top decay chain
        BRANCH(gen_bbar_lepton_tbar_deltaR, float); // DeltaR between the b quark and the lepton coming from the top decay chain

        // These two vectors are indexed wrt LepLepId enum
        BRANCH(gen_b_deltaR, std::vector<std::vector<float>>); // DeltaR between the gen b coming from the top decay and each selected jets. Indexed as `selectedJets_tightID_DRcut` array
        BRANCH(gen_bbar_deltaR, std::vector<std::vector<float>>); // DeltaR between the gen bbar coming from the anti-top decay chain and each selected jets. Indexed as `selectedJets_tightID_DRcut` array

        // These two vectors are indexed wrt LepLepId enum
        BRANCH(gen_b_beforeFSR_deltaR, std::vector<std::vector<float>>); // DeltaR between the gen b coming from the top decay and each selected jets. Indexed as `selectedJets_tightID_DRcut` array
        BRANCH(gen_bbar_beforeFSR_deltaR, std::vector<std::vector<float>>); // DeltaR between the gen bbar coming from the anti-top decay chain and each selected jets. Indexed as `selectedJets_tightID_DRcut` array

        BRANCH(gen_lepton_t_deltaR, std::vector<float>); // DeltaR between the gen lepton coming from the top decay chain and each selected lepton. Indexed as `leptons` array
        BRANCH(gen_lepton_tbar_deltaR, std::vector<float>); // DeltaR between the gen lepton coming from the anti-top decay chain and each selected lepton. Indexed as `leptons` array

        // These two vectors are indexed wrt LepLepId enum
        BRANCH(gen_matched_b, std::vector<int8_t>); // Index inside the `selectedJets_tightID_DRcut` collection of the jet with the smallest deltaR with the gen b coming from the top decay
        BRANCH(gen_matched_bbar, std::vector<int8_t>); // Index inside the `selectedJets_tightID_DRcut` collection of the jet with the smallest deltaR with the gen bbar coming from the anti-top decay

        // These two vectors are indexed wrt LepLepId enum
        BRANCH(gen_matched_b_beforeFSR, std::vector<int8_t>); // Index inside the `selectedJets_tightID_DRcut` collection of the jet with the smallest deltaR with the gen b coming from the top decay
        BRANCH(gen_matched_bbar_beforeFSR, std::vector<int8_t>); // Index inside the `selectedJets_tightID_DRcut` collection of the jet with the smallest deltaR with the gen bbar coming from the anti-top decay

        BRANCH(gen_matched_lepton_t, int16_t); // Index inside the `leptons` collection of the lepton with the smallest deltaR with the gen lepton coming from the top decay chain
        BRANCH(gen_matched_lepton_tbar, int16_t); // Index inside the `leptons` collection of the lepton with the smallest deltaR with the gen lepton coming from the anti-top decay chain

    private:

        // Producers name
        const std::string m_electrons_producer;
        const std::string m_muons_producer;
        const std::string m_jets_producer;
        const std::string m_met_producer;

        const float m_electronPtCut, m_electronEtaCut;
        const std::string m_electronVetoIDName;
        const std::string m_electronLooseIDName;
        const std::string m_electronMediumIDName;
        const std::string m_electronTightIDName;

        const float m_muonPtCut, m_muonEtaCut, m_muonLooseIsoCut, m_muonTightIsoCut;

        const float m_jetPtCut, m_jetEtaCut, m_bJetEtaCut, m_jetPUID, m_jetDRleptonCut;
        const std::string m_jetID, m_jetCSVv2Name;
        const float m_jetCSVv2L, m_jetCSVv2M, m_jetCSVv2T;

        const float m_hltDRCut, m_hltDPtCut;

        static inline bool muonIDAccessor(const MuonsProducer& muons, const uint16_t index, const std::string& muonID){
            if(index >= muons.p4.size())
              throw edm::Exception(edm::errors::StdException, "Invalid muon index passed to ID accessor");

            if(muonID == "loose")
                return muons.isLoose[index];
            
            if(muonID == "medium")
                return muons.isMedium[index];
            
            if(muonID == "tight")
                return muons.isTight[index];
            
            throw edm::Exception(edm::errors::NotFound, "Unknown muonID passed to analyzer");
        }
        
        static inline bool jetIDAccessor(const JetsProducer& jets, const uint16_t index, const std::string& jetID){
            if(index >= jets.p4.size())
              throw edm::Exception(edm::errors::StdException, "Invalid jet index passed to ID accessor");
            
            if(jetID == "loose")
                return jets.passLooseID[index];
            
            if(jetID == "tight")
                return jets.passTightID[index];
            
            if(jetID == "tightLeptonVeto")
                return jets.passTightLeptonVetoID[index];
            
            throw edm::Exception(edm::errors::NotFound, "Unknown jetID passed to analyzer");
        }
};

