from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import sys

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

from time import sleep

config = config()

submitVersion = "DZFilterEfficiency"
  
mainOutputDir = '/store/group/phys_smp/cojacob/Eff/%s' % submitVersion

config.General.transferOutputs = True
config.General.transferLogs    = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = 'EleNtupler_cfg.py'
config.Data.allowNonValidInputDataset = False
config.JobType.sendExternalFolder     = True

config.Data.inputDBS = 'global'
config.Data.publication = False

#config.Data.publishDataName = 

config.Site.storageSite = 'T2_CH_CERN'

config.General.workArea = 'the_crab_%s' % submitVersion

def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

##### submit MC
config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'MC')
config.Data.splitting     = 'FileBased'
config.Data.unitsPerJob   = 8
config.JobType.pyCfgParams  = ['GlobalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v6']

print "--- Submitting ttbar madgraph ---"
config.General.requestName  = 'ttbar_madgraph'
config.Data.inputDataset    = '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
submit(config)

print "--- Submitting WWTo2L2Nu powheg ---"
config.General.requestName  = "WWToL2Nu_powheg_actual"
config.Data.inputDataset    = '/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
submit(config)

print "--- Submitting DYtoLL_amcatNLO ---"
config.General.requestName  = 'DYToLL_amcAtNLO'
config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
submit(config)

print "--- Submitting DYtoLL_madgraph ---"
config.General.requestName  = 'DYToLL_madgraph_Moriond17'
config.Data.inputDataset    = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
submit(config)

#config.General.requestName  = 'WJets_madgraph'
#config.Data.inputDataset    = '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/MINIAODSIM'
#submit(config)

print "--- Submitting DYtoEE_powheg_m50_120 ---"
config.General.requestName  = 'DYToEE_powheg_m50_120'
config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_50_120/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
submit(config)

print "--- Submitting DYtoEE_powheg_m120_200 ---"
config.General.requestName  = 'DYToEE_powheg_m120_200'
config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_120_200/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
submit(config)

print "--- Submitting DYtoEE_powheg_m200_400 ---"
config.General.requestName  = 'DYToEE_powheg_m200_400'
config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_200_400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
submit(config)

print "--- Submitting DYtoEE_powheg_m400_800 ---"
config.General.requestName  = 'DYToEE_powheg_m400_800'
config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_400_800/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
submit(config)

print "--- Submitting DYtoEE_powheg_m800_1400 ---"
config.General.requestName  = 'DYToEE_powheg_m800_1400'
config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_800_1400/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
submit(config)

print "--- Submitting DYtoEE_powheg_m1400_2300 ---"
config.General.requestName  = 'DYToEE_powheg_m1400_2300'
config.Data.inputDataset    = '/ZToEE_NNPDF30_13TeV-powheg_M_1400_2300/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
submit(config)

# #sys.exit(0)

##### now submit SingleElectron DATA
config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'SingleElectron')
config.Data.splitting     = 'LumiBased'
config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.unitsPerJob   = 100
config.JobType.pyCfgParams  = ['GlobalTag=80X_dataRun2_2016SeptRepro_v7']

print "--- Submitting Run2016B ---"
config.General.requestName  = '2016rereco_RunB_SE'
config.Data.inputDataset    = '/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD'
submit(config)

print "--- Submitting Run2016C ---"
config.General.requestName  = '2016rereco_RunC_SE'
config.Data.inputDataset    = '/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016D ---"
config.General.requestName  = '2016rereco_RunD_SE'
config.Data.inputDataset    = '/SingleElectron/Run2016D-03Feb2017-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016E ---"
config.General.requestName  = '2016rereco_RunE_SE'
config.Data.inputDataset    = '/SingleElectron/Run2016E-03Feb2017-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016F ---"
config.General.requestName  = '2016rereco_RunF_SE'
config.Data.inputDataset    = '/SingleElectron/Run2016F-03Feb2017-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016G ---"
config.General.requestName  = '2016rereco_RunG_SE'
config.Data.inputDataset    = '/SingleElectron/Run2016G-03Feb2017-v1/MINIAOD'
submit(config)

config.JobType.pyCfgParams = ['GlobalTag=80X_dataRun2_Prompt_v16']
print "--- Submitting Run2016H2 ---"
config.General.requestName  = '2016prompt_RunHv2_SE'
config.Data.inputDataset    = '/SingleElectron/Run2016H-03Feb2017_ver2-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016H3 ---"
config.General.requestName  = '2016prompt_RunHv3_SE'
config.Data.inputDataset    = '/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD'
submit(config)

##### now submit DoubleEG DATA
config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'DoubleEG')
config.Data.splitting     = 'LumiBased'
config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.unitsPerJob   = 100
config.JobType.pyCfgParams  = ['GlobalTag=80X_dataRun2_2016SeptRepro_v7']

print "--- Submitting Run2016B ---"
config.General.requestName  = '2016rereco_RunB_DE'
config.Data.inputDataset    = '/DoubleEG/Run2016B-03Feb2017_ver2-v2/MINIAOD'
submit(config)

print "--- Submitting Run2016C ---"
config.General.requestName  = '2016rereco_RunC_DE'
config.Data.inputDataset    = '/DoubleEG/Run2016C-03Feb2017-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016D ---"
config.General.requestName  = '2016rereco_RunD_DE'
config.Data.inputDataset    = '/DoubleEG/Run2016D-03Feb2017-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016E ---"
config.General.requestName  = '2016rereco_RunE_DE'
config.Data.inputDataset    = '/DoubleEG/Run2016E-03Feb2017-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016F ---"
config.General.requestName  = '2016rereco_RunF_DE'
config.Data.inputDataset    = '/DoubleEG/Run2016F-03Feb2017-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016G ---"
config.General.requestName  = '2016rereco_RunG_DE'
config.Data.inputDataset    = '/DoubleEG/Run2016G-03Feb2017-v1/MINIAOD'
submit(config)

config.JobType.pyCfgParams = ['GlobalTag=80X_dataRun2_Prompt_v16']
print "--- Submitting Run2016H2 ---"
config.General.requestName  = '2016prompt_RunHv2_DE'
config.Data.inputDataset    = '/DoubleEG/Run2016H-03Feb2017_ver2-v1/MINIAOD'
submit(config)

print "--- Submitting Run2016H3 ---"
config.General.requestName  = '2016prompt_RunHv3_DE'
config.Data.inputDataset    = '/DoubleEG/Run2016H-03Feb2017_ver3-v1/MINIAOD'
submit(config)
