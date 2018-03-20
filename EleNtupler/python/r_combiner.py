import os
from ROOT import *
from subprocess import Popen, PIPE, call
from time import sleep

top_dir  = "Data/"
sub_dirs = ["MC","SingleElectron"]#,"DoubleEG"]

for_semi = "/afs/cern.ch/work/c/cojacob/public/ForSemiray"
here_data = "/afs/cern.ch/work/c/cojacob/public/CMSSW_8_0_27/src/EleNtupler/EleNtupler/python/Data"

def MakeEfficiencyHistos(rfile):
  tfile = TFile(rfile,"read")

  eelp = TH2F()
  eelp = tfile.Get("EtaEtaLooseElePass")
  eelt = TH2F()
  eelt = tfile.Get("EtaEtaLooseEleTotal")
  eeleff = TEfficiency(eelp, eelt)
  eeleff.SetTitle("DZ filt eff #eta #eta, loose ID")

  dzlp = TH1F()
  dzlp = tfile.Get("DZLooseElePass")
  dzlt = TH1F()
  dzlt = tfile.Get("DZLooseEleTotal")
  dzleff = TEfficiency(dzlp, dzlt)
  dzleff.SetTitle("DZ filt eff #Delta z, loose ID")

  eemp = TH2F()
  eemp = tfile.Get("EtaEtaMedElePass")
  eemt = TH2F()
  eemt = tfile.Get("EtaEtaMedEleTotal")
  eemeff = TEfficiency(eemp, eemt)
  eemeff.SetTitle("DZ filt eff #eta #eta, med ID")

  dzmp = TH1F()
  dzmp = tfile.Get("DZMedElePass")
  dzmt = TH1F()
  dzmt = tfile.Get("DZMedEleTotal")
  dzmeff = TEfficiency(dzmp, dzmt)
  dzmeff.SetTitle("DZ filt eff #Delta z, loose ID")

  dzmgp = TH1F()
  dzmgp = tfile.Get("DZMedEleGsfPass")
  dzmgt = TH1F()
  dzmgt = tfile.Get("DZMedEleGsfTotal")
  dzmgeff = TEfficiency(dzmgp, dzmgt)
  dzmgeff.SetTitle("DZ filt eff #Delta z, loose ID")

  eetp = TH2F()
  eetp = tfile.Get("EtaEtaTightElePass")
  eett = TH2F()
  eett = tfile.Get("EtaEtaTightEleTotal")
  eeteff = TEfficiency(eetp, eett)
  eeteff.SetTitle("DZ filt eff #eta #eta, tight ID")

  dztp = TH1F()
  dztp = tfile.Get("DZTightElePass")
  dztt = TH1F()
  dztt = tfile.Get("DZTightEleTotal")
  dzteff = TEfficiency(dztp, dztt)
  dzteff.SetTitle("DZ filt eff #Delta z, loose ID")

  tfile.Close()
  outfile = rfile[:-5]
  outfile += "_Eff.root"
  print outfile
  out = TFile(outfile,"recreate")
  ROOT.gFile = out

  eeleff.Write()
  dzleff.Write()
  
  eemeff.Write()
  dzmeff.Write()
  dzmgeff.Write()

  eeteff.Write()
  dzteff.Write()

  out.Close()
#  tfile.Close()
  return

def Combine(directory,outdir,dz_name,boson_name,exclude=[]):
  dz_root_files = []
  boson_root_files = []
  for path, dirs, files in os.walk(directory):
    for exc in exclude:
      if exc in path:
        continue
    for f in files:
      if f.endswith(".root"):
        file_path = path + os.sep + f
        if "DZEff" in f:
          dz_root_files.append(file_path)
        if "BosonStrat" in f:
          boson_root_files.append(file_path)

  if not dz_name.endswith(".root"):
    dz_name += ".root"
  if not boson_name.endswith(".root"):
    boson_name += ".root"
  if not outdir.endswith(os.sep):
    outdir += os.sep

  dz_command = ["hadd", "-f", outdir+dz_name] + dz_root_files
  process = Popen(dz_command, stdout=PIPE, stderr=PIPE)
  stdout, stderr = process.communicate()
  print stdout
  print stderr

  boson_command = ["hadd", "-f", outdir+boson_name] + boson_root_files
  process = Popen(boson_command, stdout=PIPE, stderr=PIPE)
  stdout, stderr = process.communicate()
  print stdout
  print stderr
  return

def main():
  for sub in sub_dirs:
    Combine(top_dir+sub, here_data, "DZFilterEff"+sub, "MultiBosonStrat"+sub)

if __name__ == "__main__":
  main()
