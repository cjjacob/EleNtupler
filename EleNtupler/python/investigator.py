import os
from ROOT import *
from efficiency import *

def tree_loop(tree,max_events=-1):
  count = 0
  for ev in tree:
    if not count < max_events:
      break
    count += 1
    if ev.nEle > 1:
      Zlist = filter(GoodZMass, MakeZs(ev))
      if len(Zlist) > 0:
        for Z in Zlist:
          filts = ev.eleFiredHLTFilters
          if abs(ev.eleZ[Z.eles()[0]] - ev.eleZ[Z.eles()[1]]) > 0.2:
            print "Z with mass %5.3f made from electrons %d and %d" % (Z.M(), Z.eles()[0], Z.eles()[1])
            print "Vertex z:   %8.5f" % ev.vtz
            print "BeamSpot z: %8.5f" % ev.bsz
            for e in range(ev.nEle): #Z.eles():
              filtInfo = ev.eleFiredHLTFilters[e]
              print "-Ele %d z: %8.5f, Gsf z: %8.5f, Q: %2.0f, DZ: %d" % (e, ev.eleZ[e], ev.gsfTrackZ[e], ev.eleCharge[e], (filts[e]>>2)&1)
            PairwiseMass(ev)
            print "Z Eles Dz:  %8.5f" % abs(ev.eleZ[Z.eles()[0]] - ev.eleZ[Z.eles()[1]])
            print "Z e-Gsf Dz: %8.5f" % abs(ev.gsfTrackZ[Z.eles()[0]] - ev.gsfTrackZ[Z.eles()[1]])
            print ""

def PairwiseMass(event):
  for e0 in range(event.nEle-1):
    for e1 in range(e0+1,event.nEle):
      ve0 = TLorentzVector()
      ve0.SetPtEtaPhiM(event.elePt[e0], event.eleEta[e0], event.elePhi[e0], 0.000511)
      ve1 = TLorentzVector()
      ve1.SetPtEtaPhiM(event.elePt[e1], event.eleEta[e1], event.elePhi[e1], 0.000511)
      zv = ve0 + ve1
      print " -- Mass of pair %d and %d: %5.3f, with charge %2.0f" % (e0, e1, zv.M(), event.eleCharge[e0]+event.eleCharge[e1])

def main():
  tfile = TFile("electronNtupler.root","read")
  tree = tfile.Get("/EleNtupler/EventTree")
  tree_loop(tree,500)
  return

if __name__ == "__main__":
  try:
    main()
  except KeyboardInterrupt:
    print "\nKeyboardInterrupt\n"
