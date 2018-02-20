import sys
import os
from subprocess import Popen, PIPE, call
import itertools
import ROOT
import threading

import efficiency as ef

eos_base  = "/eos/cms/store/group/phys_smp/cojacob/Eff/DZFilterEfficiency"
root_base = "root://eoscms.cern.ch//eos/cms/store/group/phys_smp/cojacob/Eff/DZFilterEfficiency"
root_start = "root://eoscms.cern.ch/"

jobs = ["MC", "SingleElectron"] #, "DoubleEG"]

sub_dirs = dict()
sub_dirs["MC"] = ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", "ZToEE_NNPDF30_13TeV-powheg_M_50_120", "ZToEE_NNPDF30_13TeV-powheg_M_120_200", "ZToEE_NNPDF30_13TeV-powheg_M_200_400", "ZToEE_NNPDF30_13TeV-powheg_M_400_800", "ZToEE_NNPDF30_13TeV-powheg_M_800_1400", "ZToEE_NNPDF30_13TeV-powheg_M_1400_2300","WWTo2L2Nu_13TeV-powheg"]
sub_dirs["SingleElectron"] = ["Run2016B", "Run2016C", "Run2016D", "Run2016E", "Run2016F", "Run2016G", "Run2016Hv2", "Run2016Hv3"]
#sub_dirs["DoubleEG"] = ["Run2016B", "Run2016C", "Run2016D", "Run2016E", "Run2016F", "Run2016G", "Run2016Hv2", "Run2016Hv3"]

class eos_searcher:
  def __init__(self):
    self.eos_base = "/eos/cms/store/group/phys_smp/cojacob/Eff"
    self.root_base = "root://eoscms.cern.ch//eos/cms/store/group/phys_smp/cojacob/Eff"
    self.root_start = "root://eoscms.cern.ch/"
    self.jobs = ["MC", "SingleElectron"]#, "DoubleEG"]
    self.sub_dirs = dict()
    self.sub_dirs["MC"] = ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", "ZToEE_NNPDF30_13TeV-powheg_M_50_120", "ZToEE_NNPDF30_13TeV-powheg_M_120_200", "ZToEE_NNPDF30_13TeV-powheg_M_200_400", "ZToEE_NNPDF30_13TeV-powheg_M_400_800", "ZToEE_NNPDF30_13TeV-powheg_M_800_1400", "ZToEE_NNPDF30_13TeV-powheg_M_1400_2300","WWTo2L2Nu_13TeV-powheg"]
    self.sub_dirs["SingleElectron"] = ["crab_2016rereco_RunB_SE", "crab_2016rereco_RunC_SE", "crab_2016rereco_RunD_SE", "crab_2016rereco_RunE_SE", "crab_2016rereco_RunF_se", "crab_2016rereco_RunG_SE", "crab_2016prompt_RunHv2_SE", "crab_2016rereco_RunHv3_SE"]
    #self.sub_dirs["DoubleEG"] = ["crab_2016rereco_RunB_DE", "crab_2016rereco_RunC_DE", "crab_2016rereco_RunD_DE", "crab_2016rereco_RunE_DE", "crab_2016rereco_RunF_DE", "crab_2016rereco_RunG_DE", "crab_2016prompt_RunHv2_DE", "crab_2016prompt_RunHv3_DE"]
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
  print "Done!"
  return out

def merge_root_files(root_files,outdir,name):
  print "Merging files..."
  command = ["hadd", "-f", outdir + os.sep + name + ".root"] + root_files
  print command
  process = Popen(command, stdout=PIPE, stderr=PIPE)
  stdout, stderr = process.communicate()
  print stdout
  print stderr
  print "Done!"
  return

def pack_runner(pack_of_files,output,unique):
  count = 0
  if output[-5:] == ".root":
    output = output[:-5] + "_" + str(unique) + "_" + str(count).zfill(3) + ".root"
  else:
    output += "_" + str(unique) + "_" + str(count).zfill(3) + ".root"
  for f in pack_of_files:
    ef.DoDZEffFromFile(f,output)
    count += 1
  print "\nFinished analyzing ",
  print pack_of_files
  print "Find the output in /afs/cern.ch/work/c/cojacob/public/ForSemiray/Data/\n"
  return

def main():
#  threads = []
  searcher = eos_searcher()
  for job in searcher.jobs:
    # this process takes a while, so the poor man's solution is to run it in three terminals, one on each job
    if job == "DoubleEG" or job == "SingleElectron":
#    if job == "MC" or job == "DoubleEG":
#    if job == "MC" or job == "SingleElectron":
      continue
    count = 0
    for subdir in searcher.sub_dirs[job]:
      print " --- Analyzing .root files in ",
      loc = eos_base + os.sep + job + os.sep
      if job == "MC":
        loc += subdir + os.sep
      else:
        loc += job + os.sep + subdir + os.sep # because crab.
      print loc
      files = find_root_files(loc)
  #     step = max(1,len(files)//20)
  #     split_files = [files[i:i+step] for i in range(0, len(files), step)]
  #     out_name = "/afs/cern.ch/work/c/cojacob/public/ForSemiray/Data/" + job + os.sep + job
  #     for pack in split_files:
  #       #pack_runner(pack,out_name,count)
  #       t = threading.Thread(target=pack_runner, args=(pack,out_name,count,))
  #       threads.append(t)
  #       t.start()
  #       count += 1
  # for thr in threads:
  #   thr.join()
  #   print "\n\n -- All done! -- \n"
  # return
      for f in files:
        print "Analyzing %s... " % f,
        tf = ROOT.TFile(f,"read")
        tree = tf.Get("/EleNtupler/EventTree")
        output_filename = "Data" + os.sep + job + os.sep + job + str(count).zfill(3) + ".root"
        if os.path.isfile(output_filename):
          print "already analyzed in %s" % output_filename
        else:
          ef.DoDZEff(tree,output_filename)
          print "Results in %s" % output_filename
        count += 1
#      print files
#      merge_root_files(files, searcher.here+os.sep+"Data"+os.sep+job, subdir)

if __name__ == "__main__":
  try:
    main()
  except KeyboardInterrupt:
    print "\nKeyboardInterrupt!\n"
