import sys, os
sys.path.append('{CMS}/src/RecoTauTag/RecoTau/test'.format(CMS = os.path.expandvars('$CMSSW_BASE')))

import FWCore.ParameterSet.Config as cms
from rerunTauRecoOnMiniAOD import process

readFiles = cms.untracked.vstring(
    'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/856887E2-9EC4-2740-B0ED-42E70BF1FFE6.root'   ## test on DY file
)
secFiles  = cms.untracked.vstring()

process.source = cms.Source( 'PoolSource', fileNames=readFiles, secondaryFileNames=secFiles)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32( 3000 ))

process.selectedPatTaus.cut = cms.string("pt > 18.")

## try to get rid of DM
#process.requireDecayMode.decayMode.cut = cms.double(-999999)
# single thread for debugging
# process.options.numberOfThreads = cms.untracked.uint32(1)

process.combinatoricRecoTaus.builders[0].signalConeSize = cms.string('max(min(0.1, 4.528/(pt()^0.8982)), 0.03)')    ## 0p95 quantile
# process.combinatoricRecoTaus.builders[0].verbosity = cms.int32(3)

process.output.fileName = cms.untracked.string('miniAOD_rerunTauRECO_DMFtest2.root')
