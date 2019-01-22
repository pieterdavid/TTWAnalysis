#!/usr/bin/env python2

import argparse
parser = argparse.ArgumentParser(description="Generate synchronization JSON from llbb Framework + TTWAnalysis trees")
parser.add_argument("inputfile", type=str, nargs=1, help="Input file with tree")
parser.add_argument("-t", "--treename", type=str, default="framework/t", help="Name of the tree inside the file")
parser.add_argument("-n", "--numevents", type=int, default=-1, help="Number of events to process (default: all)")
parser.add_argument("-o", "--output", type=str, default="test.json", help="Output filename")
args = parser.parse_args()

import ROOT
tup = ROOT.TChain("t")
for infn in args.inputfile:
    tup.Add(infn)

import functools
def getCont(tree, index, prefix, varName):
    return getattr(tree, "{0}_{1}".format(prefix, varName))[index]

nEvts = args.numevents
if nEvts < 0:
    nEvts = tup.GetEntries()
events = []
for i in xrange(nEvts):
    tup.GetEntry(i)

    evt = {"id" : (tup.event_run, tup.event_lumi, tup.event_event)}

    nLooseLep = tup.ttW_lepton_idx.size()
    evt["muons"] = []
    evt["electrons"] = []
    for iL in xrange(nLooseLep):
        li = tup.ttW_lepton_idx[iL]
        isMu = tup.ttW_lepton_isMu[iL]
        getL = functools.partial(getCont, tup, li, "muon" if isMu else "electron")
        p4 = getL("p4")
        lepton = {"p4": (p4.Pt(), p4.Eta(), p4.Phi(), p4.E())}
        lepton["IP"] = {
                "dxy": getL("dxy"), "dz": getL("dz"),
                "sip3d": abs(getL("dca"))
                }
        lepton["ISO"] = {
                "rel": getL("relativeIsoR03_withEA"),
                "miniCh" : getL("miniRelIsoCharged"),
                "miniNeu": getL("miniRelIsoNeutral")
                }
        lepton["JET"] = {
                "trackMult": int(getL("jetNDaugMVASel")),
                "ptRel": getL("jetPtRelv2"),
                "ptRatio": getL("jetPtRatio"),
                "Btag": getL("jetBTagCSV")
                }
        lepton["MVA"] = {
                "ttVtZq": getL("LeptonMVAGhent")
                }
        if isMu:
            lepton["ID"] = { "segmentCompatibility": getL("segmentCompatibility") }
            evt["muons"].append(lepton)
        else:
            lepton["ID"] = { "MVAIdForLepMVA": getL("MVAIdForLepMVA") }
            evt["electrons"].append(lepton)

    evt["jets"] = []
    for iJ in xrange(tup.jet_p4.size()):
        getJ = functools.partial(getCont, tup, iJ, "jet")
        if getJ("IDLoose"):
            p4 = getJ("p4")
            jet = {"p4": (p4.Pt(), p4.Eta(), p4.Phi(), p4.E())}
            ## TODO add more
            evt["jets"].append(jet)

    events.append(evt)

import json
with open(args.output, "w") as outF:
    json.dump({"Events": events}, outF)
