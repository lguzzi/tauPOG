DATASET         = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
FOLDER_NAME     = ''
OUTPUT_NAME     = ''
DESTINATION     = ''
PROJECT_NAME    = ''

file_list = os.popen('dasgoclient --query="file dataset={DST}"'.format( DST = DATASET)).readlines()
file_list = [ff.strip('\n') for ff in file_list]
import pdb; pdb.set_trace()