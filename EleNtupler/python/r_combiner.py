import os
import ROOT
from subprocess import Popen, PIPE, call

top_dir  = "Data/"
sub_dirs = ["MC","SingleElectron","DoubleEG"]

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
  command = ["hadd", "-f", outdir+name] + root_files
  print command
  process = Popen(command, stdout=PIPE, stderr=PIPE)
  stdout, stderr = process.communicate()
  print stdout
  print stderr
  return

def main():
  for sub in sub_dirs:
    Combine(top_dir+sub, top_dir, sub)

if __name__ == "__main__":
  main()
