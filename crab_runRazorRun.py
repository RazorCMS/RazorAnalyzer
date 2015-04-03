from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = False # turn on to get the crab logs
config.General.workArea = 'crab_projects'
config.General.requestName = 'sampleName'

config.section_('JobType')
config.JobType.psetName = 'pset_Razor_analysis.py' # empty pset for crab to produce a valid job report
config.JobType.pluginName = 'PrivateMC'
config.JobType.scriptExe = 'runRazorCrab.sh'
config.JobType.outputFiles = ['ntupleName.root']          # output file name from RazorRun
config.JobType.inputFiles = ['RazorRun','listfile.txt']

config.section_('Data')
config.Data.outLFN = '/store/group/phys_susy/razor/run2/RazorNtupleV1.6/Run1/Test/'
config.Data.unitsPerJob = 1  # each job will process one input file
config.Data.splitting = 'EventBased'
config.Data.primaryDataset = 'MinBias' # dummy name for output directory
config.Data.publication = False
config.Data.ignoreLocality = True #enable AAA
config.Data.totalUnits = 999 # number of files in the input list

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
