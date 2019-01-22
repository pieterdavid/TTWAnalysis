#!/usr/bin/env python2
import ROOT
tup = ROOT.TChain("t")
import sys
tup.Add(sys.argv[1])

def printEventHdr(runID, lumiID, eventID):
    print "\nrun={0:d}, lumi={1:d}, event={2:d}".format(runID, lumiID, eventID)

def printMomentum(p4):
    print " - ({0:6g}, {1:6g}, {2:6g}, {3:6g})".format(p4.Pt(), p4.Eta(), p4.Phi(), p4.E())

def printIso(relIso, minIsoCh, minIsoNeu):
    print "   ISO: rel={0:6g}, miniCharged={1:6g}, miniNeutral={2:6g}".format(relIso, minIsoCh, minIsoNeu)

def printIP(dxy, dz, sip3d):
    print "   IP: dxy={0:6g}, dz={1:6g}, sip3d={2:6g}".format(dxy, dz, sip3d)

def printClosestJet(trackMult, ptrel, ptratio, btag):
    print "   Closest jet: trackMult={0:d}, pTRel={1:6g}, ptRatio={2:6g}, deepCSV={3:6g}".format(trackMult, ptrel, ptratio, btag)

iEntry = 0
nPass = 0
while nPass < 10:
    tup.GetEntry(iEntry)
    nLooseLep = tup.ttW_lepton_idx.size()
    nJet = tup.jet_p4.size() ##
    if nLooseLep >= 2:
        nPass += 1
        printEventHdr(tup.event_run, tup.event_lumi, tup.event_event)
        print("Muons")
        for iL in range(nLooseLep):
            li = tup.ttW_lepton_idx[iL]
            if tup.ttW_lepton_isMu[iL]:
                printMomentum(tup.muon_p4[li])
        print("Electrons")
        for iL in range(nLooseLep):
            li = tup.ttW_lepton_idx[iL]
            if tup.ttW_lepton_isEl[iL]:
                printMomentum(tup.electron_p4[li])
        print("Jets")
        for iJ in range(tup.jet_p4.size()):
            if tup.jet_IDLoose[iJ]:
                printMomentum(tup.jet_p4[iJ])
        print("relIso, miniIsoCharged, miniIsoNeutral")
        for iL in range(nLooseLep):
            li = tup.ttW_lepton_idx[iL]
            if tup.ttW_lepton_isEl[iL]:
                printIso(tup.electron_relativeIsoR03_withEA[li], tup.electron_miniRelIsoCharged[li], tup.electron_miniRelIsoNeutral[li])
            else:
                printIso(tup.muon_relativeIsoR03_withEA[li], tup.muon_miniRelIsoCharged[li], tup.muon_miniRelIsoNeutral[li])
        print("dxy, dz, sip3d")
        for iL in range(nLooseLep):
            li = tup.ttW_lepton_idx[iL]
            if tup.ttW_lepton_isEl[iL]:
                printIP(tup.electron_dxy[li], tup.electron_dz[li], abs(tup.electron_dca[li]))
            else:
                printIP(tup.muon_dxy[li], tup.muon_dz[li], abs(tup.muon_dca[li]))
        print("trackMultClosestJet, pTRel, ptRatio, deepCsvClosestJet")
        for iL in range(nLooseLep):
            li = tup.ttW_lepton_idx[iL]
            if tup.ttW_lepton_isEl[iL]:
                printClosestJet(int(tup.electron_jetNDaugMVASel[li]), tup.electron_jetPtRelv2[li], tup.electron_jetPtRatio[li], tup.electron_jetBTagCSV[li])
            else:
                printClosestJet(int(tup.muon_jetNDaugMVASel[li]), tup.muon_jetPtRelv2[li], tup.muon_jetPtRatio[li], tup.muon_jetBTagCSV[li])
        print("muon segmentCompatibily")
        for iL in range(nLooseLep):
            li = tup.ttW_lepton_idx[iL]
            if tup.ttW_lepton_isMu[iL]:
                print("   ({0:6g})".format(tup.muon_segmentCompatibility[li]))
        print("Electron MVA ID Summer16")
        for iL in range(nLooseLep):
            li = tup.ttW_lepton_idx[iL]
            if tup.ttW_lepton_isEl[iL]:
                print("   ({0:6g})".format(tup.electron_MVAIdForLepMVA[li]))
        print("LeptonMVAGhent16 output")
        for iL in range(nLooseLep):
            li = tup.ttW_lepton_idx[iL]
            if tup.ttW_lepton_isEl[iL]:
                print("   {0:6g}".format(tup.electron_LeptonMVAGhent[li]))
            else:
                print("   {0:6g}".format(tup.muon_LeptonMVAGhent[li]))

    iEntry += 1
