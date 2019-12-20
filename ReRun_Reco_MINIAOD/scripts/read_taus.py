'''
Loops on the events and operates the matching between reconstructed and generated taus.
It produces two flat ntuples:
    - one with an entry for each gen tau (useful for efficiencies)
    - one with an entry for each reconstructed tau (useful for fake studies)
'''
import ROOT
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes
from treeVariables import branches # here the ntuple branches are defined
from utils import isGenHadTau, finalDaughters, printer # utility functions
import sys
import argparse
import glob

parser = argparse.ArgumentParser(
		    description="Convert MiniAOD to flat ntuples!")
parser.add_argument("--type"        , default = 'TGUN')
parser.add_argument("--label"       , default = '_delme')
parser.add_argument("--collection"  , default = 'selectedPatTaus')
parser.add_argument("--file")
args = parser.parse_args()
sample = args.type
label  = args.label
file   = args.file if args.file is not None else '{}_miniAOD_rerunTauRECO{}.root'.format(sample, label) 
collection = args.collection

##########################################################################################
# initialise output files to save the flat ntuples
outfile_gen = ROOT.TFile('tau_gentau_tuple_{}{}.root'.format(sample, label), 'recreate')
ntuple_gen = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_gen = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       

outfile_jet = ROOT.TFile('tau_jet_tuple_{}{}.root'.format(sample, label), 'recreate')
ntuple_jet = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_jet = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       

##########################################################################################
# Get ahold of the events
print 'INFO: reading file', file 
events = Events(file)

maxevents = -1 # max events to process
totevents = events.size() # total number of events in the files

# PAT taus
#label_taus = ('slimmedTaus')
print 'INFO: reading', collection
label_taus  = (collection)
handle_taus = Handle('std::vector<pat::Tau>')
# PAT jets
label_jets = ('slimmedJets')
handle_jets = Handle('std::vector<pat::Jet>')
# gen particles
label_gen  = ('prunedGenParticles')
handle_gen = Handle('std::vector<reco::GenParticle>')
# vertices
handle_vtx = Handle('std::vector<reco::Vertex>')
label_vtx  = ('offlineSlimmedPrimaryVertices')

#deeptau_ids = ['byDeepTau2017v2VSjetraw'] +   ['by%sDeepTau2017v2VSjet' 
#                                                    %ii for ii in  ['VVVLoose', 'VVLoose', 'VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']
#                                                ]

##########################################################################################
# example histogram
histos = OrderedDict()
histos['pull_pt'] = ROOT.TH1F('pull_pt', 'pull_pt', 50, -2, 2)
histos['pull_pt'].GetXaxis().SetTitle('(p_{T}^{off} - p_{T}^{gen})/p_{T}^{gen}')
histos['pull_pt'].GetYaxis().SetTitle('counts')

##########################################################################################
# start looping on the events
for i, ev in enumerate(events):
    
    ######################################################################################
    # controls on the events being processed
    if maxevents>0 and i>maxevents:
        break
        
    if i%100==0:
        print '===> processing %d / %d event' %(i, totevents)
    
    ######################################################################################
    # access the taus
    ev.getByLabel(label_taus, handle_taus)
    taus = handle_taus.product()
    # loosely filter the reco taus 
    taus = [tau for tau in taus if tau.pt()>18.]
    #for tt in taus:
    #    for ii in deeptau_ids:
    #        setattr(tt, ii, tt.tauID(ii))
    ######################################################################################
    # access the jets
    ev.getByLabel(label_jets, handle_jets)
    jets = handle_jets.product()

    # loosely filter jets
    jets = [jet for jet in jets if jet.pt()>18. and abs(jet.eta())<2.5]

    ######################################################################################
    # access the vertices
    ev.getByLabel(label_vtx, handle_vtx)
    vertices = handle_vtx.product()

    ######################################################################################
    # access the gen taus
    ev.getByLabel (label_gen, handle_gen)
    gen_particles = handle_gen.product()
    
    # select only hadronically decaying taus
    gen_taus = [pp for pp in gen_particles if abs(pp.pdgId())==15 and \
                pp.status()==2 and isGenHadTau(pp)]

    # determine gen decaymode
    for gg in gen_taus:
        gg.decayMode = tauDecayModes.genDecayModeInt([d for d in finalDaughters(gg) \
                                                      if abs(d.pdgId()) not in [12, 14, 16]])

    ######################################################################################
    # match reco taus to gen taus
    for tt in taus    : tt.gen_tau  = None # first initialise the matching to None
    for gg in gen_taus: gg.reco_tau = None # first initialise the matching to None
    
    gen_taus_copy = gen_taus # we'll cyclically remove any gen taus that gets matched
    
    for tt in taus:
        tt.gamma_reco_min_pt = min([dd.pt() for dd in finalDaughters(tt) if abs(dd.pdgId()) == 22] + [-99])
        matches = [gg for gg in gen_taus_copy if deltaR(tt.p4(), gg.visp4)<0.3]
        if not len(matches):
            continue
        matches.sort(key = lambda gg : deltaR(tt.p4(), gg.visp4))
        bestmatch = matches[0]
        tt.gen_tau = bestmatch
        bestmatch.reco_tau = tt
        gen_taus_copy = [gg for gg in gen_taus_copy if gg != bestmatch]

    for gg in gen_taus:
        gg.signal_vis_dR     = max( [deltaR(gg.visp4, dd.p4()) for dd in finalDaughters(gg) if abs(dd.pdgId()) not in [12, 14, 16]])
        gg.signal_all_dR     = max( [deltaR(gg.p4(), dd.p4())  for dd in finalDaughters(gg) if abs(dd.pdgId()) not in [12, 14, 16]])
        gg.signal_vis_dR_pt1 = max( [deltaR(gg.visp4, dd.p4()) for dd in finalDaughters(gg) if dd.pt() > 1 and  abs(dd.pdgId()) not in [12, 14, 16]] + [-99])
        gg.signal_all_dR_pt1 = max( [deltaR(gg.p4(), dd.p4())  for dd in finalDaughters(gg) if dd.pt() > 1 and  abs(dd.pdgId()) not in [12, 14, 16]] + [-99])

        gg.gamma_gen_min_pt  = min([dd.pt() for dd in finalDaughters(gg) if abs(dd.pdgId()) == 22] + [-99])
    for tt in taus:
        tt.signal_dR     = max( [deltaR(tt.p4(), dd.p4()) for dd in finalDaughters(tt) if abs(dd.pdgId()) not in [12, 14, 16]] + [-99])

    ######################################################################################
    # match reco taus to reco jets
    for jj in jets : jj.tau = None # first initialise the matching to None

    taus_copy = taus # we'll cyclically remove any tau that gets matched

    for jj in jets:
        matches = [tt for tt in taus_copy if deltaR(jj.p4(), tt.p4())<0.3]
        if not len(matches):
            continue
        matches.sort(key = lambda tt : deltaR(jj.p4(), tt.p4()))
        bestmatch = matches[0]
        jj.tau = bestmatch
        taus_copy = [tt for tt in taus_copy if tt != bestmatch]

    ######################################################################################
    # fill histograms
    for gg in gen_taus:
        if hasattr(gg, 'reco_tau') and gg.reco_tau:
            pull = (gg.reco_tau.pt() - gg.vispt())/gg.vispt()
            histos['pull_pt'].Fill(pull)
    ######################################################################################
    # fill the ntuple: each gen tau makes an entry
    for gg in gen_taus:
        for k, v in tofill_gen.iteritems(): tofill_gen[k] = -99. # initialise before filling
        tofill_gen['run'               ] = ev.eventAuxiliary().run()
        tofill_gen['lumi'              ] = ev.eventAuxiliary().luminosityBlock()
        tofill_gen['event'             ] = ev.eventAuxiliary().event()
        tofill_gen['nvtx'              ] = vertices.size()
        tofill_gen['gamma_gen_min_pt'  ] = gg.gamma_gen_min_pt
        if hasattr(gg, 'reco_tau') and gg.reco_tau:
            tofill_gen['tau_reco_mass'     ] = gg.reco_tau.mass()
            tofill_gen['tau_reco_pt'       ] = gg.reco_tau.pt()
            tofill_gen['tau_reco_eta'      ] = gg.reco_tau.eta()
            tofill_gen['tau_reco_phi'      ] = gg.reco_tau.phi()
            tofill_gen['tau_reco_charge'   ] = gg.reco_tau.charge()
            tofill_gen['tau_reco_signal_dR'] = gg.reco_tau.signal_dR
            tofill_gen['tau_reco_decaymode'] = gg.reco_tau.decayMode()
            tofill_gen['tau_reco_iso'      ] = gg.reco_tau.tauID('byLooseIsolationMVArun2v1DBoldDMwLT')
            tofill_gen['gamma_reco_min_pt' ] = gg.reco_tau.gamma_reco_min_pt
            #for ii in deeptau_ids:
            #    tofill_gen['tau_reco_'+ii] = getattr(gg.reco_tau, ii, -99)
        tofill_gen['tau_gen_pt'        ] = gg.pt()
        tofill_gen['tau_gen_eta'       ] = gg.eta()
        tofill_gen['tau_gen_phi'       ] = gg.phi()
        tofill_gen['tau_gen_charge'    ] = gg.charge()
        tofill_gen['tau_gen_decaymode' ] = gg.decayMode
        tofill_gen['tau_gen_vis_mass'  ] = gg.vismass()
        tofill_gen['tau_gen_vis_pt'    ] = gg.vispt()
        tofill_gen['tau_gen_vis_eta'   ] = gg.viseta()
        tofill_gen['tau_gen_vis_phi'   ] = gg.visphi()
        tofill_gen['tau_gen_vis_signal_dR'] = gg.signal_vis_dR
        tofill_gen['tau_gen_all_signal_dR'] = gg.signal_all_dR
        tofill_gen['tau_gen_vis_signal_dR_pt1'] = gg.signal_vis_dR
        tofill_gen['tau_gen_all_signal_dR_pt1'] = gg.signal_all_dR
        ntuple_gen.Fill(array('f',tofill_gen.values()))
    # fill the ntuple: each jet makes an entry
    for jj in jets:
        for k, v in tofill_jet.iteritems(): tofill_jet[k] = -99. # initialise before filling
        tofill_jet['run'        ] = ev.eventAuxiliary().run()
        tofill_jet['lumi'       ] = ev.eventAuxiliary().luminosityBlock()
        tofill_jet['event'      ] = ev.eventAuxiliary().event()
        tofill_jet['nvtx'       ] = vertices.size()
        tofill_jet['jet_mass'   ] = jj.mass()
        tofill_jet['jet_pt'     ] = jj.pt()
        tofill_jet['jet_eta'    ] = jj.eta()
        tofill_jet['jet_phi'    ] = jj.phi()
        tofill_jet['jet_charge' ] = jj.charge()
        if hasattr(jj, 'tau') and jj.tau:
            #for ii in deeptau_ids:|
            #    tofill_jet['tau_reco_'+ii] = getattr(jj.tau, ii, -99)
            tofill_jet['tau_reco_mass'     ] = jj.tau.mass()
            tofill_jet['tau_reco_pt'       ] = jj.tau.pt()
            tofill_jet['tau_reco_eta'      ] = jj.tau.eta()
            tofill_jet['tau_reco_phi'      ] = jj.tau.phi()
            tofill_jet['tau_reco_charge'   ] = jj.tau.charge()
            tofill_jet['tau_reco_decaymode'] = jj.tau.decayMode()
            tofill_gen['gamma_reco_min_pt' ] = jj.tau.gamma_reco_min_pt
            tofill_gen['tau_reco_iso'      ] = jj.tau.tauID('byLooseIsolationMVArun2v1DBoldDMwLT')
            if hasattr(jj.tau, 'gen_tau') and jj.tau.gen_tau:
                tofill_gen['tau_gen_pt'        ] = jj.tau.gen_tau.pt()
                tofill_gen['tau_gen_eta'       ] = jj.tau.gen_tau.eta()
                tofill_gen['tau_gen_phi'       ] = jj.tau.gen_tau.phi()
                tofill_gen['tau_gen_charge'    ] = jj.tau.gen_tau.charge()
                tofill_gen['tau_gen_decaymode' ] = jj.tau.gen_tau.decayMode
                tofill_gen['tau_gen_vis_mass'  ] = jj.tau.gen_tau.vismass()
                tofill_gen['tau_gen_vis_pt'    ] = jj.tau.gen_tau.vispt()
                tofill_gen['tau_gen_vis_eta'   ] = jj.tau.gen_tau.viseta()
                tofill_gen['tau_gen_vis_phi'   ] = jj.tau.gen_tau.visphi()
        ntuple_jet.Fill(array('f',tofill_jet.values()))
    ######################################################################################
    # printout some info, if you want
    # printer(taus, gen_taus)

##########################################################################################
# write the ntuples and close the files
outfile_gen.cd()
ntuple_gen.Write()
outfile_gen.Close()

outfile_jet.cd()
ntuple_jet.Write()
outfile_jet.Close()
