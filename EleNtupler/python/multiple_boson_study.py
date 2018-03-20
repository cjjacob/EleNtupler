from ROOT import *
from array import array
import os

massBins  = array('d', [float(m) for m in range(61,122)])
nMassBins = len(massBins)-1

elePtCut = 0.

ZMass = 91.19
eleMass = 0.000511

root_start = "root://eoscms.cern.ch/"

class eos_searcher:
  def __init__(self):
    self.eos_base  = "/eos/cms/store/group/phys_smp/cojacob/Eff/DZFilterEfficiency"
    self.root_base = "root://eoscms.cern.ch//eos/cms/store/group/phys_smp/cojacob/Eff/DZFilterEfficiency"
    self.root_start = "root://eoscms.cern.ch/"
    self.jobs = ["SingleElectron"]
    self.sub_dirs = dict()
    self.sub_dirs["SingleElectron"] = ["crab_2016rereco_RunB_SE", "crab_2016rereco_RunC_SE", "crab_2016rereco_RunD_SE", "crab_2016rereco_RunE_SE", "crab_2016rereco_RunF_se", "crab_2016rereco_RunG_SE", "crab_2016prompt_RunHv2_SE", "crab_2016prompt_RunHv3_SE"]
    self.here = "/afs/cern.ch/work/c/cojacob/public/CMSSW_8_0_27/src/EleNtupler/EleNtupler/python"

def find_root_files(directory):
  print "Looking for .root files in %s..." % directory,
  out = []
  for subdir, dirs, files in os.walk(directory):
#    print "Searching " + subdir + "..."
    if "failed" in subdir:
      continue
    for f in files:
      if f.endswith(".root"):
#        print f
        file_path = root_start + subdir + os.sep + f
        out.append(file_path)
#        print file_path
  print "Done!"
  return out

def run_on_eos():
  searcher = eos_searcher()
  for job in searcher.jobs:
    count = 0
    for subdir in searcher.sub_dirs[job]:
      loc = searcher.eos_base + os.sep + job + os.sep + job + os.sep + subdir + os.sep
      print " -- Analyzing .root files in %s" % loc
      files = find_root_files(loc)
      for f in files:
        print " -- -- Working on %s... " % f,
        tf = TFile(f,"read")
        tree = tf.Get("/EleNtupler/EventTree")
        outname = "Data" + os.sep + job + os.sep + "MultiBosonStrat_" + str(count).zfill(3) + ".root"
        if os.path.isfile(outname):
          print "already analyzed in %s." % outname
        else:
          print "find the results in %s soon... " % outname
          CompareBosonStrategies(tree,outname)
        count += 1

def InvSqMassDiff(mz):
    return 1.0/((mz-ZMass)*(mz-ZMass))

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

def MakeZs(event, ptCut=0.0):
  Zlist = []
  if event.nEle > 1:
    for e0 in range(event.nEle-1):
      for e1 in range(e0+1,event.nEle):
        if (event.eleCharge[e0] * event.eleCharge[e1]) < 0.: # opp charge
          if (event.eleIDbit[e0]&2)>0 and (event.eleIDbit[e1]&2)>0: # loose id
            if event.elePt[e0] > ptCut and event.elePt[e1] > ptCut: # above pt cut
              ev0 = TLorentzVector()
              ev0.SetPtEtaPhiM(event.elePt[e0], event.eleEta[e0], event.elePhi[e0], eleMass)
              ev1 = TLorentzVector()
              ev1.SetPtEtaPhiM(event.elePt[e1], event.eleEta[e1], event.elePhi[e1], eleMass)
              Z = ev0 + ev1
              Zlist.append( ZwithEles(e0,e1,Z.M()) )
#    Zlist.sort(CompByZMass)
  return Zlist

def CompareBosonStrategies(tree,outfile="BosonStrats.root",maxEv=-1):
  out = TFile(outfile,"recreate")

  HighestPtOnly  = TH1F("HighestPtOnly",  "HighestPtOnly;m_{Z} (GeV)",  nMassBins, massBins);
  AllEqualWeight = TH1F("AllEqualWeight", "AllEqualWeight;m_{Z} (GeV)", nMassBins, massBins);
  InvSqWeight    = TH1F("InvSqWeight",    "InvSqWeight;m_{Z} (GeV)",    nMassBins, massBins);
  BestMassOnly   = TH1F("BestMassOnly",   "BestMassOnly;m_{Z} (GeV)",   nMassBins, massBins);

  count = 0
  at_least_z = 0
  multi_zs = 0

  for event in tree:
    if maxEv != -1 and not count < maxEv:
      break
    if event.nEle > 1:
      Zlist = MakeZs(event)
      if len(Zlist) > 0:
        at_least_z += 1
        HighestPtOnly.Fill(Zlist[0].M())
        eq_weight = 1.0/len(Zlist)
        inv_sq_w = [InvSqMassDiff(Z.M()) for Z in Zlist]
        sum_inv_sq_w = reduce( lambda x,y: x+y, inv_sq_w )
        inv_sq_w = map( lambda w: w/sum_inv_sq_w, inv_sq_w )
        for iz, Z in enumerate(Zlist):
          AllEqualWeight.Fill(Z.M(), eq_weight)
          InvSqWeight.Fill(Z.M(), inv_sq_w[iz])
        Zlist.sort(CompByZMass)
        BestMassOnly.Fill(Zlist[0].M())
        if len(Zlist) > 1:
          multi_zs += 1
    count += 1

  ROOT.gFile = out
  
  HighestPtOnly.Write()
  AllEqualWeight.Write()
  InvSqWeight.Write()
  BestMassOnly.Write()

  out.Close()

  print "At least 1 Z happend in %d of %d events (%4.2f%%)" % (at_least_z, count, 100.*float(at_least_z)/count)
  print "Multiple Zs happened in %d of %d events with a Z (%4.2f%%)" % (multi_zs, at_least_z, 100.*float(multi_zs)/at_least_z)

def main():
#  print "Starting..."
#  tfile = TFile("electronNtupler.root","read")
#  tree = tfile.Get("/EleNtupler/EventTree")
#  CompareBosonStrategies(tree)
  run_on_eos()

if __name__ == "__main__":
  try:
    main()
  except KeyboardInterrupt:
    print "\nKeyboardInterrupt\n"
