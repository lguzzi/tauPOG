submit_cfg = '''
import sys, os
sys.path.append('{CMS}/src/RecoTauTag/RecoTau/test'.format(CMS = os.path.expandvars('$CMSSW_BASE')))

import FWCore.ParameterSet.Config as cms
from rerunTauRecoOnMiniAOD import process

readFiles = cms.untracked.vstring(
    {FIL}
)
secFiles  = cms.untracked.vstring()

process.source = cms.Source( 'PoolSource', fileNames=readFiles, secondaryFileNames=secFiles)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32( 3000 ))

process.selectedPatTaus.cut = cms.string("pt > 18.")
process.combinatoricRecoTaus.builders[0].signalConeSize = cms.string('max(min(0.1, 4.528/(pt()^0.8982)), 0.03)')    ## 0p95 quantile
process.output.fileName = cms.untracked.string('{OUT}')
'''
import pdb; pdb.set_trace()