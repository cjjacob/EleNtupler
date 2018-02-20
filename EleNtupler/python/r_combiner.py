import os
from ROOT import *
from subprocess import Popen, PIPE, call

top_dir  = "Data/"
sub_dirs = ["MC","SingleElectron"]#,"DoubleEG"]

for_semi = "/afs/cern.ch/work/c/cojacob/public/ForSemiray"
here_data = "/afs/cern.ch/work/c/cojacob/public/CMSSW_8_0_27/src/EleNtupler/EleNtupler/python/Data"

def MakeEfficiencyHistos(rfile):
  tfile = TFile(rfile,"update")

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
  dzmp = tfile.Get("DZLooseElePass")
  dzmt = TH1F()
  dzmt = tfile.Get("DZLooseEleTotal")
  dzmeff = TEfficiency(dzlp, dzlt)
  dzmeff.SetTitle("DZ filt eff #Delta z, loose ID")

  dzmgp = TH1F()
  dzmgp = tfile.Get("DZLooseElePass")
  dzmgt = TH1F()
  dzmgt = tfile.Get("DZLooseEleTotal")
  dzmgeff = TEfficiency(dzlp, dzlt)
  dzmgeff.SetTitle("DZ filt eff #Delta z, loose ID")

  eetp = TH2F()
  eetp = tfile.Get("EtaEtaTightElePass")
  eett = TH2F()
  eett = tfile.Get("EtaEtaTightEleTotal")
  eeteff = TEfficiency(eetp, eett)
  eeteff.SetTitle("DZ filt eff #eta #eta, tight ID")

  dztp = TH1F()
  dztp = tfile.Get("DZLooseElePass")
  dztt = TH1F()
  dztt = tfile.Get("DZLooseEleTotal")
  dzteff = TEfficiency(dzlp, dzlt)
  dzteff.SetTitle("DZ filt eff #Delta z, loose ID")

#  tfile.Close()
#  outfile = rfile[:-5]
#  outfile += "_Eff.root"
#  out = TFile(outfile,"recreate")
#  ROOT.gFile = out

  eeleff.Write()
  dzleff.Write()
  
  eemeff.Write()
  dzmeff.Write()
  dzmgeff.Write()

  eeteff.Write()
  dzteff.Write()

#  out.Close()
  tfile.Close()
  return

def Combine(directory,outdir,name,exclude=[]):
  root_files = []
  for path, dirs, files in os.walk(directory):
    for exc in exclude:
      if exc in path:
        continue
    for f in files:
      if f.endswith(".root"):
        file_path = path + os.sep + f
        root_files.append(file_path)

  if not name.endswith(".root"):
    name += ".root"
  if not outdir.endswith(os.sep):
    outdir += os.sep
  command = ["hadd", "-f", outdir+name] + root_files[:-1]
#  print command
  process = Popen(command, stdout=PIPE, stderr=PIPE)
  stdout, stderr = process.communicate()
  print stdout
  print stderr
  MakeEfficiencyHistos(outdir+name)
  return

def main():
  for sub in sub_dirs:
#    Combine(top_dir+sub, for_semi, "DZfilterEff"+sub)
    Combine(top_dir+sub, here_data, "DZfilterEff"+sub)

if __name__ == "__main__":
  main()
