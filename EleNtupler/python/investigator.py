import os
from ROOT import *
from efficiency import *
from math import sqrt, pi

def tree_loop_high_dz(tree,max_events=-1):
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
          e0, e1 = Z.eles()
          if abs(ev.eleZ[e0] - ev.eleZ[e1]) > 0.2:
            print "Z with mass %5.3f made from electrons %d and %d.    (%4d)" % (Z.M(), e0, e1, count)
            print "Vertex z:   %8.5f" % ev.vtz
            print "BeamSpot z: %8.5f" % ev.bsz
            for e in range(ev.nEle): #Z.eles():
              print "-Ele %d z: %8.5f, Gsf z: %8.5f, Q: %2.0f, DZ: %d" % (e, ev.eleZ[e], ev.gsfTrackZ[e], ev.eleCharge[e], (filts[e]>>2)&1)
              print "  Ele dz: %8.5f, (ele z - vtx z): %8.5f" % (ev.eleDz[e], ev.eleZ[e]-ev.vtz)
            PairwiseMass(ev)
            print "Z Eles Dz:  %8.5f" % abs(ev.eleZ[e0] - ev.eleZ[e1])
            print "Z e-Gsf Dz: %8.5f" % abs(ev.gsfTrackZ[e0] - ev.gsfTrackZ[e1])
            print "Z D dz:     %8.5f" % abs(ev.eleDz[e0] - ev.eleDz[e1])
            shifted0 = ev.vtz - ev.bsz
            shifted1 = -shifted0
            dz0 = [ev.eleZ[e] - shifted0 for e in Z.eles()]
            dz1 = [ev.eleZ[e] - shifted1 for e in Z.eles()]
            print "Shifted0 dz's: %8.5f, %8.5f; diff: %8.5f" % (dz0[0], dz0[1], abs(dz0[0]-dz0[1]))
            print "Shifted1 dz's: %8.5f, %8.5f; diff: %8.5f" % (dz1[0], dz1[1], abs(dz1[0]-dz1[1]))
            print ""

def PairwiseMass(event):
  for e0 in range(event.nEle-1):
    for e1 in range(e0+1,event.nEle):
      if event.eleCharge[e0]*event.eleCharge[e1] < 0.:
        ve0 = TLorentzVector()
        ve0.SetPtEtaPhiM(event.elePt[e0], event.eleEta[e0], event.elePhi[e0], 0.000511)
        ve1 = TLorentzVector()
        ve1.SetPtEtaPhiM(event.elePt[e1], event.eleEta[e1], event.elePhi[e1], 0.000511)
        zv = ve0 + ve1
        vsc0 = TLorentzVector()
        vsc0.SetPtEtaPhiM(event.elePt[e0], event.eleSCEta[e0], event.eleSCPhi[e0], 0.000511)
        vsc1 = TLorentzVector()
        vsc1.SetPtEtaPhiM(event.elePt[e1], event.eleSCEta[e1], event.eleSCPhi[e1], 0.000511)
        scv = vsc0 + vsc1
        dz0 = (event.eleFiredHLTFilters[e0]>>2) & 1
        dz1 = (event.eleFiredHLTFilters[e1]>>2) & 1
        dxy = sqrt( (event.eleX[e0]-event.eleX[e1])*(event.eleX[e0]-event.eleX[e1]) + (event.eleY[e0]-event.eleY[e1])*(event.eleY[e0]-event.eleY[e1]) ) 
        Dz = abs(event.eleZ[e0]-event.eleZ[e1])
        DR = DeltaR(event.eleEta[e0],event.elePhi[e0],event.eleEta[e1],event.elePhi[e1])
        print " :::: Mass of pair %d and %d: %7.3f (%7.3f SC), DZ: %d & %d, Delta z: %6.5f, Delta xy: %6.5f --- Delta R: %6.4f" % (e0, e1, zv.M(), scv.M(), dz0, dz1, Dz, dxy, DR)

def DeltaPhi(phi0, phi1):
  dphi = phi0 - phi1
  if dphi < -pi:
    dphi += 2*pi
  if dphi > pi:
    dphi = 2*pi - dphi
  return dphi

def DeltaR(eta0, phi0, eta1, phi1):
  Dphi = DeltaPhi(phi0,phi1)
  Deta = eta0 - eta1
  return sqrt(Deta*Deta + Dphi*Dphi)

def tree_loop_low_dz(tree,max_events=-1):
  count = 0
  for ev in tree:
    if not count < max_events:
      break
    count += 1
    if ev.nEle > 1:
      Zlist = filter(GoodZMass, MakeZs(ev))
      filts = ev.eleFiredHLTFilters
      any_leg2 = 0
      for e in range(ev.nEle):
        if (filts[e] & 2) > 0:
          any_leg2 += 1
      if any_leg2 > 1:
        if len(Zlist) > 0:
          for Z in Zlist:
            e0, e1 = Z.eles()
            dz = abs(ev.eleZ[e0] - ev.eleZ[e1])
            z_dz_pass = ((filts[e0]&4)>0) or ((filts[e1]&4)>0)
            if dz < 0.2 and not z_dz_pass:
              print " *** Z failed"
            print "Vertex z:   %8.5f" % ev.vtz
            print "BeamSpot z: %8.5f" % ev.bsz
            print "Z from %d and %d. %2d vertices in event      (%4d)" % (e0, e1, ev.nVtx, count)
            for e in range(ev.nEle):
              print " -Ele %d z: %8.5f, Gsf z: %8.5f, Q: %2.0f, L1: %d, L2: %d, DZ: %d. --- (pt, eta, phi): (%5.2f, %6.3f, %6.3f)" % (e, ev.eleZ[e], ev.gsfTrackZ[e], ev.eleCharge[e], (filts[e]>>0)&1, (filts[e]>>1)&1, (filts[e]>>2)&1, ev.elePt[e], ev.eleEta[e], ev.elePhi[e])
              print " |-  DR with TrigObj: %6.4f. TO (pt,eta,phi): (%5.2f, %6.3f, %6.3f). SC (eta, phi): (%6.3f, %6.3f)" % (ev.eleMatchedObjDR[e], ev.eleMatchedObjPt[e], ev.eleMatchedObjEta[e], ev.eleMatchedObjPhi[e], ev.eleSCEta[e], ev.eleSCPhi[e])
            print ""
            PairwiseMass(ev)
            print "\n"
  return

def tree_loop_filters_only(tree,max_events=-1):
  count = 0
  for ev in tree:
    if not count < max_events:
      break
    count += 1
    any_leg2 = 0
    for e in range(ev.nEle):
      if (ev.eleFiredHLTFilters[e]&2) > 0:
        any_leg2 += 1
    if any_leg2 > 0:
      hltpass = False
      # print "%d electron(s). HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ " % ev.nEle,
      if (ev.eleHLTs & 32) > 0:
        # print "passed. ",
        hltpass = True
      else:
        # print "failed. ",
        pass
      filts = ev.eleFiredHLTFilters
      if hltpass and any_leg2 < 2:
        print "Less than 2 leg2 passes. %d passes on %d eles. HLT: 1    (%4d)" % (any_leg2, ev.nEle, count)
        print "Run %6d, Event %10ld, LumiSection %4d" % (ev.run, ev.event, ev.lumis)
        for e in range(ev.nEle):
          print " -Ele %d z: %8.5f, vtz: %8.5f, Q: %2.0f, L1: %d, L2: %d, DZ: %d. --- (pt, eta, phi): (%5.2f, %6.3f, %6.3f)" % (e, ev.eleZ[e], ev.vtz, ev.eleCharge[e], (filts[e]>>0)&1, (filts[e]>>1)&1, (filts[e]>>2)&1, ev.elePt[e], ev.eleEta[e], ev.elePhi[e])
        print ""
      if any_leg2 > 2:
        print "More than 2 leg2 passes. %d passes on %d eles. HLT: %d    (%4d)" % (any_leg2, ev.nEle, (ev.eleHLTs>>5)&1, count)
        print "Run %6d, Event %10ld, LumiSection %4d" % (ev.run, ev.event, ev.lumis)
        for e in range(ev.nEle):
          print " -Ele %d z: %8.5f, vtz: %8.5f, Q: %2.0f, L1: %d, L2: %d, DZ: %d. --- (pt, eta, phi): (%5.2f, %6.3f, %6.3f)" % (e, ev.eleZ[e], ev.vtz, ev.eleCharge[e], (filts[e]>>0)&1, (filts[e]>>1)&1, (filts[e]>>2)&1, ev.elePt[e], ev.eleEta[e], ev.elePhi[e])
        if (ev.event == 257980647):
          PairwiseMass(ev)
        print ""
      print "%d electron(s) passed leg 2 filter.    (%4d)" % (any_leg2, count)
      for e in range(ev.nEle):
        filt = ev.eleFiredHLTFilters[e]
        print "--Ele %d, Leg1: %d, Leg2: %d, DZ: %d" % (e, (filt)&1, (filt>>1)&1, (filt>>2)&1)
      print ""
  return

def main():
  tfile = TFile("electronNtupler.root","read")
  tree = tfile.Get("/EleNtupler/EventTree")
#  tree_loop_high_dz(tree,5000)
  print " -*-*- ---------------------------------------- -*-*-\n"
  tree_loop_low_dz(tree,5000)
  print " -*-*- ---------------------------------------- -*-*-\n"
#  tree_loop_filters_only(tree,1000)
  return

if __name__ == "__main__":
  try:
    main()
  except KeyboardInterrupt:
    print "\nKeyboardInterrupt\n"
