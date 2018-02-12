import sys
import os
from subprocess import Popen, PIPE, call
import itertools
import ROOT

eos_base  = "/eos/cms/store/group/phys_smp/cojacob/Eff"
root_base = "root://eoscms.cern.ch//eos/cms/store/group/phys_smp/cojacob/Eff"

jobs = ["MC", "SingleElectron", "DoubleEG"]

sub_dirs = dict()
sub_dirs["MC"] = ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", "ZToEE_NNPDF30_13TeV-powheg_M_50_120", "ZToEE_NNPDF30_13TeV-powheg_M_120_200", "ZToEE_NNPDF30_13TeV-powheg_M_200_400", "ZToEE_NNPDF30_13TeV-powheg_M_400_800", "ZToEE_NNPDF30_13TeV-powheg_M_800_1400", "ZToEE_NNPDF30_13TeV-powheg_M_1400_2300","WWTo2L2Nu_13TeV-powheg"]
sub_dirs["SingleElectron"] = ["Run2016B", "Run2016C", "Run2016D", "Run2016E", "Run2016F", "Run2016G", "Run2016Hv2", "Run2016Hv3"]
sub_dirs["DoubleEG"] = ["Run2016B", "Run2016C", "Run2016D", "Run2016E", "Run2016F", "Run2016G", "Run2016Hv2", "Run2016Hv3"]

def find_root_files(directory):
  print "Looking for .root files in %s..." % directory
  out = []
  for subdir, dirs, files in os.walk(directory):
    print "Searching " + subdir + "..."
    if "failed" in subdir:
      continue
    for f in files:
      if f.endswith(".root"):
#        print f
        file_path = subdir + os.sep + f
        out.append(file_path)
  print "Done!"
  return out

def merge_root_files(root_files,outdir,name):
  print "Merging files..."
  command = ["hadd", "-f", eos_base + os.sep + outdir + os.sep + name + ".root"] + root_files
  print command
  process = Popen(command, stdout=PIPE, stderr=PIPE)
  stdout, stderr = process.communicate()
  print stdout
  print stderr
  print "Done!"
  return

def main():
  for job in jobs:
    for subdir in sub_dirs[job]:
      print " --- Combining .root files in ",
      loc = eos_base + os.sep + filt + os.sep + job + os.sep + subdir + os.sep
      print loc
      files = find_root_files(loc)
      merge_root_files(files, eos_base+os.sep+filt+os.sep+job, subdir)

if __name__ == "__main__":
  main()
