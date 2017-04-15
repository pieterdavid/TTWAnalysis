# Analyzers for the ttW analysis #

Plugins and configuration to save information related to electrons, muons and jets, and di-lepton+di-jet candidates made from miniAOD inputs to small trees.

See the [cp3-llbb framework](https://github.com/cp3-llbb/Framework) project for an overview of framework, and the [HH](https://github.com/cp3-llbb/HHAnalysis) and [ZA](https://github.com/cp3-llbb/ZAAnalysis) repositories for similar analyzers.

## Overview of the code ##

The structure of the ttW configuration differs a bit from the other analyses.
The main idea behind the differences is to make the code more modular, such that it is more transparent and easier to adapt to new datasets, selections and final states.

The typical monolithic analyzer was split into a `TTWAnalyzer`, which is only responsible for applying selections and constructing candidates, a `TTWTruthAnalyzer` for the branches related to MC truth, and `TTWAnalysis::DictTool<objecttype>` plugins (where objecttype is `pat::Electron`, `pat::Muon` or `pat::Jet` for now) that each fill a number of branches (e.g. kinematic quantities, identification variables, trigger matching...).

The configuration of the path and analyzers can be found in [`Configuration.py`](python/Configuration.py).

### Output tree structure ###

The only classes used in the output tree are `TTWAnalysis::GenParticle` (similar to `TTAnalysis::GenParticle`), `LorentzVector` (`ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float>>`), and (`std::vector` of) float, integer and boolean types.

Because of this, the code that reads the trees only needs to load a small class dictionary that does not need frequent updates.
It does, however, need to re-construct the 'objects' out of the different branches --- this is done by following the common convention of a prefix ending with an underscore (see [ttWtools](https://gitlab.cern.ch/piedavid/ttWTools) for more details on how this can be automatically handled).
The main advantage of this is that a branch (e.g. a new ID variable) can be added with a single line in a `DictTool` (and no changes are need to be done to the downstream code until it is used there).

### `TTWAnalyzer` ###

Moving the truth matching to [`TTWTruthAnalyzer`](plugins/TTWTruthAnalyzer.cc) and HLT matching to [`HLTDictTools`](plugins/HLTDictTools.cc) already makes the [`TTWAnalyzer`](plugins/TTWAnalyzer.cc) quite a bit smaller.
In addition, the definition of object selection working points is left to the configuration: a dictionary of working point names and [generic cut expression](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePhysicsCutParser) strings is taken, and the lists of objects and the indices of those that are selected, are available through public methods, e.g. [`const edm::PtrVector<pat::Electron>& TTWAnalyzer::getPtrList<pat::Electron>() const`](interface/TTWAnalyzer.h#L212) and [`const TTWAnalysis::indexlist_t& TTWAnalyzer::selectedElectrons(const std::string& wpName) const`](interface/TTWAnalyzer.h#L127) (these could be moved to an interface class).

### `Producers` ###

Since the object information is added to the tree by a `DictTool`, the producers do not do anything that depends on the physics object type, and can be replaced by a trivial generic version, [`TTWAnalysis::CandidatesProducer<ObjectType>`](interface/CandidatesProducer.h) (plugin definitions in [plugins/CandidateProducers.cc](plugins/CandidatesProducers.cc)).

To avoid duplicate work and ensure a consistent definition of variables that are used to derive other quantities, `DictTool` plugins can also be added to the producers, which will add the return values as `userFloat` to the objects (this is used to store the impact parameter (significance) for the leptons, which is used directly and as an input to the lepton MVA).
This very basic handling of dependencies between the `DictTool` plugins is a limitation of the current implementation.

### `DictTool<objecttype>` ###

The `DictTool` plugins fill most of the branches, see [ `BasicDictTools`](plugins/BasicDictTools.cc) for some related to the candidates, and the `*DictTools.cc` files in [`plugins`](plugins) for more.

They are added to the path through a [`TTWAnalysis::DelegatingAnalyzer`](plugins/DelegatingAnalyzer.cc), which simply calls a number of [`TTWAnalysis::AnalyzerHelper`](interface/AnalyzerHelper.h) plugins.
One of these, which adds the return values of a list of `DictTool`s for the corresponding object to `std::vector` branches, is added for each object type (`reco::Vertex`, `pat::Electron`, `pat::Muon`, `pat::Jet`, `TTWAnalysis::Lepton`, `TTWAnalysis::DiLepton`, `TTWAnalysis::DiJet`, `TTWAnalysis::DiLeptonDiJet` and `TTWAnalysis::DiLeptonDiJetMet` --- the candidate classes in the `TTWAnalysis` namespace only keep track of the references, and are not written to the output tree); the implementation is generic (templated on the object type and the return type of `TTWAnalyzer::getPtrList<objecttype>` --- the second argument could be avoided through type erasure, e.g. with `boost::any_range` when the Boost version in CMSSW is updated).

### Configuration ###

The configuration in [`python/Configuration.py`](python/Configuration.py) is kept in one file, and intentionally kept in a few rather deeply nested declarations.
This hopefully makes it easy to keep most of the inputs and settings closely together, and to find back how the values for a specific branch are calculated.

### Scale factors ###

From the `CMSSW_8_0_X` version onwards, scale factors are not stored in the trees.
The rationale behind this is that they are one of the things that can change frequently (and at any point in the analysis), and they can simply be looked up from a `json` file based on quantities that are stored.
This also saves space in the trees.
See [`scalefactors.py`](https://gitlab.cern.ch/piedavid/ttWTools/tree/dev/python/scalefactors.py) and [`ScaleFactors.h`](https://gitlab.cern.ch/piedavid/ttWTools/tree/dev/common/include/ScaleFactors.h) in [ttWtools](https://gitlab.cern.ch/piedavid/ttWTools) for more information on how this is dealt with later on.
