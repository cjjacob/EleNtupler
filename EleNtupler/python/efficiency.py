from ROOT import *
from array import array
from random import randint
import os

Leg1Filter = "hltEle23Ele12CaloIdLTrackIdlIsoVLTrackIsoLeg1Filter"
Leg1FilterBit = 0
Leg2Filter = "hltEle23Ele12CaloIdLTrackIdlIsoVLTrackIsoLeg2Filter"
Leg2FilterBit = 1
DZFilter = "hltEle23Ele12CaloIdLTrackIdlIsoVLDZFilter"
DZFilterBit = 2
Ele27Filter = "hltSingleEle27WPTightGsfTrackIsoFilter"
Ele27FilterBit = 3;

Trigger    = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"
TriggerBit = 5

ptBins    = array('d', [0., 10., 11., 12., 13., 14., 21., 22., 23., 24., 25., 30., 40.])
nPtBins   = len(ptBins)-1
etaBins   = array('d', [-2.4, -2.0, -1.6, -1.44, -1., -0.5, 0., 0.5, 1., 1.44, 1.6, 2.0, 2.4])
nEtaBins  = len(etaBins)-1
DZBins    = array('d', [0., 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26])
nDZBins   = len(DZBins)-1
massBins  = array('d', [61., 62., 63., 64., 65., 66., 67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80., 81., 82., 83., 84., 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 96., 97., 98., 99., 100., 101., 102., 103., 104., 105., 106., 107., 108., 109., 110., 111., 112., 113., 114., 115., 116., 117., 118., 119., 120., 121.])
nMassBins = len(massBins)-1

allDZBins = array('d',[-0.26, -0.24, -0.22, -0.20, -0.18, -0.16, -0.14, -0.12, -0.10, -0.08, -0.06, -0.04, -0.03, -0.02, -0.01, 0., 0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26])
nAllDZBins = len(allDZBins)-1

ZMass = 91.19 # PDG value
eleMass = 0.000511

def Test(ttree,maxEv=10):
  count = 0
  for event in ttree:
    if count >= maxEv:
      break
#    print "%ld: %d e:" % (event.event, event.nEle),
#    for pt in event.elePt:
#        print "%f " % pt,
#    print ""
#    for eid in event.eleIDbit:
#      print "Ele ID bits: medium %d, tight %d;  " % (((eid>>2) & 1), ((eid>>3) & 1 )),
#    print ""
    if event.nEle > 1:
      Zlist = filter(GoodZMass, MakeZs(event))
      if len(Zlist) > 0:
        print "Vertex z: %f, " % event.vtz
        print "Possible Z bosons: "
        for Z in Zlist:
          print "%f GeV made by eles %d and %d" % (Z.M(), Z.eles()[0], Z.eles()[1])
          for e in Z.eles():
            filtInfo = event.eleFiredHLTFilters[e]
            print "( e%d z: %f, dz: %f; pt: %f, eta: %f, phi: %f, q: %f)" % (e, event.eleZ[e], event.eleDz[e], event.elePt[e], event.eleEta[e], event.elePhi[e], event.eleCharge[e])
            print "Leg1Filter: %d, Leg2Filter: %d, DZFilter: %d\n" % (filtInfo&1, (filtInfo>>1)&1, (filtInfo>>2)&1)

    count += 1

class ZwithEles:
  def __init__(self, e0, e1, mass):
    self.e0 = e0
    self.e1 = e1
    self.mass = mass

  def M(self):
    return self.mass

  def eles(self):
    out = [self.e0]
    out.append(self.e1)
    return out

def CompByZMass(a,b):
  if abs(a.M()-ZMass) < abs(b.M()-ZMass):
    return -1
  elif abs(a.M()-ZMass) == abs(b.M()-ZMass):
    return 0
  else:
    return 1

def GoodZMass(z):
  if abs(z.M() - ZMass) <= 20.:
    return True
  else:
    return False

def MakeZs(event, ptCut=0.0):
  Zlist = []
  if event.nEle > 1:
    for e0 in range(event.nEle-1):
      for e1 in range(e0+1,event.nEle):
        if (event.eleCharge[e0] * event.eleCharge[e1]) < 0.:
          if event.elePt[e0] > ptCut and event.elePt[e1] > ptCut:
            ev0 = TLorentzVector()
            ev0.SetPtEtaPhiM(event.elePt[e0], event.eleEta[e0], event.elePhi[e0], eleMass)
            ev1 = TLorentzVector()
            ev1.SetPtEtaPhiM(event.elePt[e1], event.eleEta[e1], event.elePhi[e1], eleMass)
            Z = ev0 + ev1
            Zlist.append( ZwithEles(e0,e1,Z.M()) )
    Zlist.sort(CompByZMass)
  return Zlist

def DoDZEffFromFile(tfile,output):
  if os.path.isfile(output):
    return
  else:
    infile = TFile(tfile,"read")
    tree = infile.Get("/EleNtupler/EventTree")

    outfile = TFile(output, "recreate")
    ROOT.gFile = outfile

    hMedEleP = TH2F("EtaEtaMedElePass",  "MedElePass;#eta_{1};#eta_{2}", nEtaBins, etaBins, nEtaBins, etaBins)
    hMedEleT = TH2F("EtaEtaMedEleTotal", "MedEleTot;#eta_{1};#eta_{2}",  nEtaBins, etaBins, nEtaBins, etaBins)

    hDZMedP = TH1F("DZMedElePass",  "DZMedElePass;abs(z_{e_{0}}-z_{e_{1}}) (cm)",  nDZBins, DZBins)
    hDZMedT = TH1F("DZMedEleTotal", "DZMedEleTotal;abs(z_{e_{0}}-z_{e_{1}}) (cm)", nDZBins, DZBins)

    hDZMedGsfP = TH1F("DZMedEleGsfPass",  "DZMedEleGsfPass;abs(z_{e_{0}}-z_{e_{1}}) (cm)",  nDZBins, DZBins)
    hDZMedGsfT = TH1F("DZMedEleGsfTotal", "DZMedEleGsfTotal;abs(z_{e_{0}}-z_{e_{1}}) (cm)", nDZBins, DZBins)

    hTightEleP = TH2F("EtaEtaTightElePass",  "TightElePass;#eta_{1};#eta_{2}", nEtaBins, etaBins, nEtaBins, etaBins)
    hTightEleT = TH2F("EtaEtaTightEleTotal", "TightEleTot;#eta_{1};#eta_{2}",  nEtaBins, etaBins, nEtaBins, etaBins)

    hDZTightP = TH1F("DZTightElePass",  "DZTightElePass;abs(z_{e_{0}}-z_{e_{1}}) (cm)",  nDZBins, DZBins)
    hDZTightT = TH1F("DZTightEleTotal", "DZTightEleTotal;abs(z_{e_{0}}-z_{e_{1}}) (cm)", nDZBins, DZBins)

    for ev in tree:
      if ev.nEle > 1:
        Zlist = filter(GoodZMass, MakeZs(ev,25.0))
        if len(Zlist) > 0:
          for Z in Zlist:
            e0, e1 = Z.eles()
            if abs(ev.eleZ[e0] - ev.eleZ[e1]) <= 0.26: #opening up range ## pass definition of DZ filter
              if randint(0,1) == 0: # randomly swap the electrons so e0 isn't always the leading pt ele
                e0, e1 = e1, e0
              filt0 = ev.eleFiredHLTFilters[e0]
              filt1 = ev.eleFiredHLTFilters[e1]
              # the electrons each have to pass one of the filters
              if (filt0 & 1 > 0 and filt1 & 2 > 0) or (filt0 & 2 > 0 and filt1 & 1 > 0):
                # the non-DZ trigger is passed
                id0 = ev.eleIDbit[e0]
                id1 = ev.eleIDbit[e1]
                if (id0>>2)&1 == 1 and (id1>>2)&1 == 1:
                  #medium id eles
                  hMedEleT.Fill(ev.eleEta[e0],ev.eleEta[e1])
                  hDZMedT.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))
                  hDZMedGsfT.Fill(abs(ev.gsfTrackZ[e0] - ev.gsfTrackZ[e1]))
                  if filt0 & 4 > 0 or filt1 & 4 > 0:
                    hMedEleP.Fill(ev.eleEta[e0],ev.eleEta[e1])
                    hDZMedP.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))
                    hDZMedGsfT.Fill(abs(ev.gsfTrackZ[e0] - ev.gsfTrackZ[e1]))
                if (id0>>3)&1 == 1 and (id1>>3)&1 == 1:
                  #tight id eles
                  hTightEleT.Fill(ev.eleEta[e0],ev.eleEta[e1])
                  hDZTightT.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))
                  if filt0 & 4 > 0 or filt1 & 4 > 0:
                    hTightEleP.Fill(ev.eleEta[e0],ev.eleEta[e1])
                    hDZTightP.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))

    hMedEleP.Write()
    hMedEleT.Write()

    hDZMedP.Write()
    hDZMedT.Write()

    hDZMedGsfP.Write()
    hDZMedGsfT.Write()

    hTightEleP.Write()
    hTightEleT.Write()

    hDZTightP.Write()
    hDZTightT.Write()

    hMedEleDZEff = TEfficiency(hMedEleP, hMedEleT)
    hTightEleDZEff = TEfficiency(hTightEleP, hTightEleT)
    hMedEleDZEff.Write()
    hTightEleDZEff.Write()

    hDZEffMed = TEfficiency(hDZMedP, hDZMedT)
    hDZEffTight = TEfficiency(hDZTightP, hDZTightT)
    hDZEffMed.Write()
    hDZEffTight.Write()

    outfile.Close()
    return

def DoDZEff(ttree,output="DZHistos.root",mode="recreate"):
  outfile = TFile(output, mode)
  ROOT.gFile = outfile

  hLooseEleP = TH2F("EtaEtaLooseElePass",  "LooseElePass;#eta_{1};#eta_{2}", nEtaBins, etaBins, nEtaBins, etaBins)
  hLooseEleT = TH2F("EtaEtaLooseEleTotal", "LooseEleTot;#eta_{1};#eta_{2}",  nEtaBins, etaBins, nEtaBins, etaBins)

  hDZLooseP = TH1F("DZLooseElePass",  "DZLooseElePass;abs(z_{e_{0}}-z_{e_{1}}) (cm)",  nDZBins, DZBins)
  hDZLooseT = TH1F("DZLooseEleTotal", "DZLooseEleTotal;abs(z_{e_{0}}-z_{e_{1}}) (cm)", nDZBins, DZBins)

  hMedEleP = TH2F("EtaEtaMedElePass",  "MedElePass;#eta_{1};#eta_{2}", nEtaBins, etaBins, nEtaBins, etaBins)
  hMedEleT = TH2F("EtaEtaMedEleTotal", "MedEleTot;#eta_{1};#eta_{2}",  nEtaBins, etaBins, nEtaBins, etaBins)

  hDZMedP = TH1F("DZMedElePass",  "DZMedElePass;abs(z_{e_{0}}-z_{e_{1}}) (cm)",  nDZBins, DZBins)
  hDZMedT = TH1F("DZMedEleTotal", "DZMedEleTotal;abs(z_{e_{0}}-z_{e_{1}}) (cm)", nDZBins, DZBins)

  hDZMedGsfP = TH1F("DZMedEleGsfPass",  "DZMedEleGsfPass;abs(z_{e_{0}}-z_{e_{1}}) (cm)",  nDZBins, DZBins)
  hDZMedGsfT = TH1F("DZMedEleGsfTotal", "DZMedEleGsfTotal;abs(z_{e_{0}}-z_{e_{1}}) (cm)", nDZBins, DZBins)

  hTightEleP = TH2F("EtaEtaTightElePass",  "TightElePass;#eta_{1};#eta_{2}", nEtaBins, etaBins, nEtaBins, etaBins)
  hTightEleT = TH2F("EtaEtaTightEleTotal", "TightEleTot;#eta_{1};#eta_{2}",  nEtaBins, etaBins, nEtaBins, etaBins)

  hDZTightP = TH1F("DZTightElePass",  "DZTightElePass;abs(z_{e_{0}}-z_{e_{1}}) (cm)",  nDZBins, DZBins)
  hDZTightT = TH1F("DZTightEleTotal", "DZTightEleTotal;abs(z_{e_{0}}-z_{e_{1}}) (cm)", nDZBins, DZBins)

  hZPtPzOppCorn = TH2F("ZPtPzOppCorn", "Z p_{T} p_{z} opposite #eta corners;p_{T} (GeV);p_{z} (GeV)", 30, 0., 300., 30, -600., 600.)
  hZPtPzSimCorn = TH2F("ZPtPzSimCorn", "Z p_{T} p_{z} similar  #eta corners;p_{T} (GeV);p_{z} (GeV)", 30, 0., 300., 30, -600., 600.)
  hZPtPzCentral = TH2F("ZPtPzCentral", "Z p_{T} p_{z} central #eta region;p_{T} (GeV);p_{z} (GeV)", 30, 0., 300., 30, -600., 600.)

  for ev in ttree:
    if ev.nEle > 1:
      for iZ in range(len(ev.ZM)):
        eta0 = ev.eleEta[ev.Ze0[iZ]]
        eta1 = ev.eleEta[ev.Ze1[iZ]]
        if (eta0 < -1. and eta1 > 1.) or (eta0 > 1. and eta1 < -1.):
          hZPtPzOppCorn.Fill(ev.ZPt[iZ], ev.ZPz[iZ])
        if (eta0 < -1. and eta1 < -1.) or (eta0 > 1. and eta1 > 1.):
          hZPtPzSimCorn.Fill(ev.ZPt[iZ], ev.ZPz[iZ])
        if (abs(eta0) < 1. and abs(eta1) < 1.):
          hZPtPzCentral.Fill(ev.ZPt[iZ], ev.ZPz[iZ])
      Zlist = filter(GoodZMass, MakeZs(ev,25.0))
      if len(Zlist) > 0:
        for Z in Zlist:
          e0, e1 = Z.eles()
          if abs(ev.eleZ[e0] - ev.eleZ[e1]) <= 0.26: #open up to fail region ## pass definition of DZ filter
            if randint(0,1) == 0: # randomly swap the electrons so e0 isn't always the leading pt ele
              e0, e1 = e1, e0
            filt0 = ev.eleFiredHLTFilters[e0]
            filt1 = ev.eleFiredHLTFilters[e1]
            # the electrons each have to pass one of the filters
            if (filt0 & 1 > 0 and filt1 & 2 > 0) or (filt0 & 2 > 0 and filt1 & 1 > 0):
              # the non-DZ trigger is passed
              id0 = ev.eleIDbit[e0]
              id1 = ev.eleIDbit[e1]
              if (id0>>1)&1 == 1 and (id1>>1)&1 == 1:
                #loose id eles
                hLooseEleT.Fill(ev.eleEta[e0],ev.eleEta[e1])
                hDZLooseT.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))
                if filt0 & 4 > 0 or filt1 & 4 > 0:
                  hLooseEleP.Fill(ev.eleEta[e0],ev.eleEta[e1])
                  hDZLooseP.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))
              if (id0>>2)&1 == 1 and (id1>>2)&1 == 1:
                #medium id eles
                hMedEleT.Fill(ev.eleEta[e0],ev.eleEta[e1])
                hDZMedT.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))
                hDZMedGsfT.Fill(abs(ev.gsfTrackZ[e0] - ev.gsfTrackZ[e1]))
                if filt0 & 4 > 0 or filt1 & 4 > 0:
                  hMedEleP.Fill(ev.eleEta[e0],ev.eleEta[e1])
                  hDZMedP.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))
                  hDZMedGsfP.Fill(abs(ev.gsfTrackZ[e0] - ev.gsfTrackZ[e1]))
              if (id0>>3)&1 == 1 and (id1>>3)&1 == 1:
                #tight id eles
                hTightEleT.Fill(ev.eleEta[e0],ev.eleEta[e1])
                hDZTightT.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))
                if filt0 & 4 > 0 or filt1 & 4 > 0:
                  hTightEleP.Fill(ev.eleEta[e0],ev.eleEta[e1])
                  hDZTightP.Fill(abs(ev.eleZ[e0] - ev.eleZ[e1]))

  hLooseEleP.Write()
  hLooseEleT.Write()

  hDZLooseP.Write()
  hDZLooseT.Write()

  hMedEleP.Write()
  hMedEleT.Write()

  hDZMedP.Write()
  hDZMedT.Write()

  hDZMedGsfP.Write()
  hDZMedGsfT.Write()

  hTightEleP.Write()
  hTightEleT.Write()

  hDZTightP.Write()
  hDZTightT.Write()

#  hMedEleDZEff = TEfficiency(hMedEleP, hMedEleT)
#  hTightEleDZEff = TEfficiency(hTightEleP, hTightEleT)
#  hMedEleDZEff.Write()
#  hTightEleDZEff.Write()

#  hDZEffMed = TEfficiency(hDZMedP, hDZMedT)
#  hDZEffTight = TEfficiency(hDZTightP, hDZTightT)
#  hDZEffMed.Write()
#  hDZEffTight.Write()

  hZPtPzOppCorn.Write()
  hZPtPzSimCorn.Write()
  hZPtPzCentral.Write()

  outfile.Close()

def DoFilters(ttree):
  outfile = TFile("FilterHistos.root", "recreate")
  ROOT.gFile = outfile

  # Leg 1 histos
  print "Making histos for %s" % Leg1Filter
  h1TrigEleP  = TH2F("Leg1TrigEleP",  "Leg1TrigEleP;p_{T} (GeV);#eta",  nPtBins, ptBins, nEtaBins, etaBins)
  h1MedEleP   = TH2F("Leg1MedEleP",   "Leg1MedEleP;p_{T} (GeV);#eta",   nPtBins, ptBins, nEtaBins, etaBins)
  h1TightEleP = TH2F("Leg1TightEleP", "Leg1TightEleP;p_{T} (GeV);#eta", nPtBins, ptBins, nEtaBins, etaBins)
  h1TrigEleF  = TH2F("Leg1TrigEleF",  "Leg1TrigEleF;p_{T} (GeV);#eta",  nPtBins, ptBins, nEtaBins, etaBins)
  h1MedEleF   = TH2F("Leg1MedEleF",   "Leg1MedEleF;p_{T} (GeV);#eta",   nPtBins, ptBins, nEtaBins, etaBins)
  h1TightEleF = TH2F("Leg1TightEleF", "Leg1TightEleF;p_{T} (GeV);#eta", nPtBins, ptBins, nEtaBins, etaBins)

  # Leg 2 histos
  print "Making histos for %s" % Leg2Filter
  h2TrigEleP  = TH2F("Leg2TrigEleP",  "Leg2TrigEleP;p_{T} (GeV);#eta",  nPtBins, ptBins, nEtaBins, etaBins)
  h2MedEleP   = TH2F("Leg2MedEleP",   "Leg2MedEleP;p_{T} (GeV);#eta",   nPtBins, ptBins, nEtaBins, etaBins)
  h2TightEleP = TH2F("Leg2TightEleP", "Leg2TightEleP;p_{T} (GeV);#eta", nPtBins, ptBins, nEtaBins, etaBins)
  h2TrigEleF  = TH2F("Leg2TrigEleF",  "Leg2TrigEleF;p_{T} (GeV);#eta",  nPtBins, ptBins, nEtaBins, etaBins)
  h2MedEleF   = TH2F("Leg2MedEleF",   "Leg2MedEleF;p_{T} (GeV);#eta",   nPtBins, ptBins, nEtaBins, etaBins)
  h2TightEleF = TH2F("Leg2TightEleF", "Leg2TightEleF;p_{T} (GeV);#eta", nPtBins, ptBins, nEtaBins, etaBins)

  # DZ histos
  print "Making histos for %s" % DZFilter
  hzTrigEleP  = TH1F("DZTrigEleP",  "DZTrigEleP;#Delta z (cm)",  nDZBins, DZBins)
  hzMedEleP   = TH1F("DZMedEleP",   "DZMedEleP;#Delta z (cm)",   nDZBins, DZBins)
  hzTightEleP = TH1F("DZTightEleP", "DZTightEleP;#Delta z (cm)", nDZBins, DZBins)
  hzTrigEleF  = TH1F("DZTrigEleF",  "DZTrigEleF;#Delta z (cm)",  nDZBins, DZBins)
  hzMedEleF   = TH1F("DZMedEleF",   "DZMedEleF;#Delta z (cm)",   nDZBins, DZBins)
  hzTightEleF = TH1F("DZTightEleF", "DZTightEleF;#Delta z (cm)", nDZBins, DZBins)

  for event in ttree:
    if event.nEle > 1:
      Zlist = filter(GoodZMass, MakeZs(event))
      for Z in Zlist:
        f1pass = False
        f2pass = False
        for e in Z.eles():
#          f1pass = False
#          f2pass = False
          if ((event.eleFiredHLTFilters[e]>>Leg1FilterBit) & 1) == 1: # leg 1 filter passed
            f1pass = True
            if ((event.eleIDbit[e]>>2) & 1) == 1:
              h1MedEleP.Fill(event.elePt[e], event.eleSCEta[e])
            if ((event.eleIDbit[e]>>3) & 1) == 1:
              h1TightEleP.Fill(event.elePt[e], event.eleSCEta[e])
            if ((event.eleFiredHLTFilters[e]>>Ele27FilterBit) & 1) == 1:
              h1TrigEleP.Fill(event.elePt[e], event.eleSCEta[e])
          else: # leg 1 filter failed
            if ((event.eleIDbit[e]>>2) & 1) == 1:
              h1MedEleF.Fill(event.elePt[e], event.eleSCEta[e])
            if ((event.eleIDbit[e]>>3) & 1) == 1:
              h1TightEleF.Fill(event.elePt[e], event.eleSCEta[e])
            if ((event.eleFiredHLTFilters[e]>>Ele27FilterBit) & 1) == 1:
              h1TrigEleF.Fill(event.elePt[e], event.eleSCEta[e])

          if event.nEle > 1: # leg 2 has a 2 electron requirement
            if ((event.eleFiredHLTFilters[e]>>Leg2FilterBit) & 1) == 1: # leg 2 filter passed
              f2pass = True
              if ((event.eleIDbit[e]>>2) & 1) == 1:
                h2MedEleP.Fill(event.elePt[e], event.eleSCEta[e])
              if ((event.eleIDbit[e]>>3) & 1) == 1:
                  h2TightEleP.Fill(event.elePt[e], event.eleSCEta[e])
              if ((event.eleFiredHLTFilters[e]>>Ele27FilterBit) & 1) == 1:
                h2TrigEleP.Fill(event.elePt[e], event.eleSCEta[e])
            else: # leg 2 filter failed
              if ((event.eleIDbit[e]>>2) & 1) == 1:
                h2MedEleF.Fill(event.elePt[e], event.eleSCEta[e])
              if ((event.eleIDbit[e]>>3) & 1) == 1:
                h2TightEleF.Fill(event.elePt[e], event.eleSCEta[e])
              if ((event.eleFiredHLTFilters[e]>>Ele27FilterBit) & 1) == 1:
                h2TrigEleF.Fill(event.elePt[e], event.eleSCEta[e])
      
        if f1pass and f2pass:    
          if ((event.eleFiredHLTFilters[e]>>DZFilterBit) & 1) == 1: # DZ filter passed
            if ((event.eleIDbit[e]>>2) & 1) == 1:
              hzMedEleP.Fill(abs(event.eleDz[e]))
            if ((event.eleIDbit[e]>>3) & 1) == 1:
              hzTightEleP.Fill(abs(event.eleDz[e]))
            if ((event.eleFiredHLTFilters[e]>>Ele27FilterBit) & 1) == 1:
              hzTrigEleP.Fill(abs(event.eleDz[e]))
            else: # DZ filter failed
              if ((event.eleIDbit[e]>>2) & 1) == 1:
                hzMedEleF.Fill(abs(event.eleDz[e]))
              if ((event.eleIDbit[e]>>3) & 1) == 1:
                hzTightEleF.Fill(abs(event.eleDz[e]))
              if ((event.eleFiredHLTFilters[e]>>Ele27FilterBit) & 1) == 1:
                hzTrigEleF.Fill(abs(event.eleDz[e]))    

  # finished loop over events
  h1MedEleP.Write()
  h1TightEleP.Write()
  h1TrigEleP.Write()
  h1MedEleF.Write()
  h1TightEleF.Write()
  h1TrigEleF.Write()
  h2MedEleP.Write()
  h2TightEleP.Write()
  h2TrigEleP.Write()
  h2MedEleF.Write()
  h2TightEleF.Write()
  h2TrigEleF.Write()
  hzMedEleP.Write()
  hzTightEleP.Write()
  hzTrigEleP.Write()
  hzMedEleF.Write()
  hzTightEleF.Write()
  hzTrigEleF.Write()
  # histogram much?
  outfile.Close()

def DoDZExploration(ttree):
  outfile = TFile("DZFilterCutExploration.root", "recreate")
  ROOT.gFile = outfile

  passing = TH2F("PassDZFilter", "PassDZFilter;#Delta z_{e_{0}} (cm);#Delta z_{e_{1}} (cm)", nAllDZBins, allDZBins, nAllDZBins, allDZBins)
  failing = TH2F("FailDZFilter", "FailDZFilter;#Delta z_{e_{0}} (cm);#Delta z_{e_{1}} (cm)", nAllDZBins, allDZBins, nAllDZBins, allDZBins)

  zcoordspass = TH2F("ZCoordsPass", "ZCoordsPass;z_{e_{0}} (cm);z_{e_{1}} (cm)", 200, -10.0, 10.0, 200, -10.0, 10.0)
  zcoordsfail = TH2F("ZCoordsFail", "ZCoordsFail;z_{e_{0}} (cm);z_{e_{1}} (cm)", 200, -10.0, 10.0, 200, -10.0, 10.0)

  dzfiltP = TH1F("DZpass", "DZpass;abs(z_{e_{0}}-z_{e_{1}}) (cm)", 25, 0., 0.25)
  dzfiltF = TH1F("DZfail", "DZfail;abs(z_{e_{0}}-z_{e_{1}}) (cm)", 25, 0., 0.25)

  for ev in ttree:
    if ev.nEles > 1:
      Zlist = filter(GoodZMass, MakeZs(ev,25.0)) # make |m_ll-m_Z|<=20GeV Zs out of pt>25GeV electrons
      if len(Zlist) > 0:
        for Z in Zlist:
          e0, e1 = Z.eles()
          if ev.eleIDbit[e0] & 4 > 0 and ev.eleIDbit[e1] & 4 > 0: # medium ID eles
            if randint(0,1) == 0: # randomly swap the electrons so e0 isn't always the leading pt ele
              e0, e1 = e1, e0
            filt0 = ev.eleFiredHLTFilters[e0]
            filt1 = ev.eleFiredHLTFilters[e1]
            if (filt0 & 1 > 0 and filt1 & 2 > 0) or (filt0 & 2 > 0 and filt1 & 1 > 0):
              if filt0 & 4 > 0 and filt1 & 4 > 0:
                passing.Fill(ev.eleDz[e0], ev.eleDz[e1])
                zcoordspass.Fill(ev.eleZ[e0], ev.eleZ[e1])
                dzfiltP.Fill(abs(ev.eleZ[e0]-ev.eleZ[e1]))
              else:
                failing.Fill(ev.eleDz[e0], ev.eleDz[e1])
                zcoordsfail.Fill(ev.eleZ[e0], ev.eleZ[e1])
                dzfiltF.Fill(abs(ev.eleZ[e0]-ev.eleZ[e1]))

  passing.Write()
  failing.Write()
  outfile.Close()

def main():
  tfile = TFile("electronNtupler.root","read")
  tree = tfile.Get("/EleNtupler/EventTree")
#  Test(tree,80)
  DoFilters(tree)
  DoDZEff(tree)
  DoDZExploration(tree)

if __name__ == "__main__":
  main()
