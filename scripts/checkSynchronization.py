#!/usr/bin/env python2

import functools
from itertools import izip
import json
import math
from pprint import pprint

def evtIdLess(evtIdA, evtIdB):
    runA,lmbA,evtA = evtIdA
    runB,lmbB,evtB = evtIdB
    if runA != runB:
        return runA < runB
    if lmbA != lmbB:
        return lmbA < lmbB
    return evtA < evtB

class Event(object):
    """ Event class: unique ID and collections """
    def __init__(self, evtID, collections=None):
        self.evtID = tuple(evtID)
        self.collections = collections if collections else dict()

    def __str__(self):
        return "Event #{0:d}:{1:d}:{2:d}".format(*self.evtID)

class PhysObject(object):
    """ Reconstructec object, comparison is done by DeltaR """
    def __init__(self, p4, stages=None):
        self.p4 = tuple(p4)
        self.stages = stages if stages else dict()
    def deltaR(self, other):
        dEta = other.p4[1]-self.p4[1]
        dPhi = other.p4[2]-self.p4[2]
        return math.sqrt(dEta**2+dPhi**2)
    def __str__(self):
        return "({0:6g}, {1:6g}, {2:6g}, {3:6g})".format(*self.p4)

def parse(syncJson):
    """ Parse JSON into objects """
    events = []
    for evtJson in syncJson["Events"]:
        evt = Event(evtJson["id"])
        for ky,val in evtJson.iteritems():
            if ky != "id":
                if isinstance(val, list):
                    objs = []
                    for obJSON in val:
                        obj = PhysObject(obJSON["p4"])
                        for oky,oval in obJSON.iteritems():
                            obj.stages[oky] = oval
                        objs.append(obj)
                    evt.collections[ky] = objs
        events.append(evt)
    return sorted(events, key=lambda evt : evt.evtID)

def isclose(a,b):
    return ( 2*abs(a-b)/(a+b) < 1.e-6 if a+b != 0. else abs(a-b) < 1.e-6 )

def fillGaps(coll1, coll2, match=(lambda x,y:x==y), comp=(lambda x,y:x<y), empty=None):
    n1,n2 = len(coll1), len(coll2)
    nc1, nc2 = [], []
    i1 = 0
    i2 = 0
    while i1 < n1 and i2 < n2:
        if match(coll1[i1], coll2[i2]):
            nc1.append(coll1[i1])
            nc2.append(coll2[i2])
            i1 += 1
            i2 += 1
        while i1 < n1 and comp(coll1[i1], coll2[i2]): ## 1-item before 2-item -> figure out where 2-item belongs
            nc1.append(coll1[i1])
            i1 += 1
            nc2.append(empty)
        while i2 < n2 and comp(coll2[i2], coll1[i1]): ## 2-item before 1-item -> figure out where 1-item belongs
            nc1.append(empty)
            nc2.append(coll2[i2])
            i2 += 1
    return nc1, nc2

def compareCollection(refColl, testColl, logger=None, pre=""):
    nDiffs = 0
    rcb, tcb = fillGaps(refColl, testColl, match=(lambda oR,oT : oR.deltaR(oT) < .01), comp=(lambda oR,oT : (oR.deltaR(oT) > .01) and (oR.p4[0] - oT.p4[0] > 1.e-6)))
    for refObj, testObj in izip(rcb, tcb):
        if testObj is None:
            if logger: logger.error("{0}{1:s} not found in test collection".format(pre, refObj))
            nDiffs += 1
        elif refObj is None:
            if logger: logger.error("{0}{1:s} not found in reference collection".format(pre, refObj))
            nDiffs += 1
        else:
            for stNm, refSt in refObj.stages.iteritems():
                if stNm not in testObj.stages:
                    nDiffs += 1
                    if logger: logger.error("No stage '{0}' in test object {1}{2:s}".format(stNm, pre, testObj))
                elif stNm == "p4":
                    testP4 = testObj.stages["p4"]
                    for ix,xNm in enumerate(("pt", "eta", "phi", "e")):
                        if not isclose(refSt[ix], testP4[ix]):
                            nDiffs += 1
                            if logger: logger.error("Different momentum component {0}: {rV} vs. {tV} for {pre}{rO} and {pre}{tO}".format(xNm, rV=refSt[ix], tV=testP4[ix], pre=pre, rO=refObj, tO=testObj))
                else:
                    testSt = testObj.stages[stNm]
                    for ky,refVal in refSt.iteritems():
                        if ky not in testSt:
                            nDiffs += 1
                            if logger: logger.error("No entry '{0}' for stage '{1}' in test object {2}{3:s}".format(ky, stNm, pre, testObj))
                        else:
                            testVal = testSt[ky]
                            if not isclose(refVal, testVal):
                                nDiffs += 1
                                if logger: logger.error("Different values for '{0}' for stage '{1}': {rV} vs. {tV} for {pre}{rO} and {pre}{tO}".format(ky, stNm, rV=refVal, tV=testVal, pre=pre, rO=refObj, tO=testObj))
                    for ky in testSt:
                        if ky not in refSt:
                            nDiffs += 1
                            if logger: logger.error("No entry '{0}' for stage '{1}' in ref object {2}{3:s}".format(ky, stNm, pre, refObj))
            for stNm in testObj.stages.iterkeys():
                if stNm != "p4" and stNm not in refObj.stages:
                    nDiffs += 1
                    if logger: logger.error("No stage '{0}' in ref object {1}{2:s}".format(stNm, pre, refObj))

    return nDiffs

def compareEvents(refEvents, testEvents, logger=None):
    """ Compare (return True if equal, False otherwise), and print a summary of differences to the logger """
    nDiffs = 0
    reb, teb = fillGaps(refEvents, testEvents, match=(lambda rEv,tEv: rEv.evtID==tEv.evtID), comp=(lambda rEv,tEv : evtIdLess(rEv.evtID, tEv.evtID)))
    for refEvt, testEvt in izip(reb, teb):
        if testEvt is None:
            if logger: logger.error("{0:s} not found in test collection".format(refEvt))
            nDiffs += 1
        elif refEvt is None:
            if logger: logger.error("{0:s} not found in reference collection".format(refEvt))
            nDiffs += 1
        else:
            evDiffs = 0
            for collName, refColl in refEvt.collections.iteritems():
                if collName != "id":
                    if collName not in testEvt.collections:
                        evDiffs += 1
                        if logger: logger.error("Collection {0} not found in test event".format(collName))
                    elif not isinstance(refColl, list):
                        if logger: logger.warning("Found non-id non-collection entry in reference: {0}".format(collName))
                    else:
                        testColl = testEvt.collections[collName]
                        if len(refColl) != len(testColl):
                            if logger: logger.error("Different length for {0}: {1:d} (ref) versus {2:d} (test)".format(collName, len(refColl), len(testColl)))
                        collDiffs = compareCollection(refColl, testColl, logger=logger, pre="{0} entry with p4=".format(collName))
                        if collDiffs:
                            if logger: logger.error("Found {0:d} differences in collection {0}".format(collDiffs, collName))
                            evDiffs += collDiffs
            for collName,testColl in testEvt.collections.iteritems():
                if isinstance(testColl, list) and collName not in refEvt.collections:
                    evDiffs += 1
                    if logger: logger.error("Collection {0} not found in reference event".format(collName))
            if logger: logger.info("{0:d} differences in {1:s}".format(evDiffs, refEvt))
            nDiffs += evDiffs
    if logger: logger.info("{0:d} differences in total".format(nDiffs))
    return nDiffs

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Test equality of synchronization JSON outputs")
    parser.add_argument("reference", type=str, nargs=1, help="Reference output")
    parser.add_argument("test", type=str, nargs=1, help="Test output")
    args = parser.parse_args()

    with open(args.reference[0]) as refF:
        ref = parse(json.load(refF))
    with open(args.test[0]) as testF:
        test = parse(json.load(testF))

    import logging
    diffLogger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    diffLogger.info("Comparing {0} with {1}".format(args.reference[0], args.test[0]))

    nDiffs = compareEvents(ref, test, logger=diffLogger)

    import sys
    sys.exit(0 if nDiffs == 0 else 1)
