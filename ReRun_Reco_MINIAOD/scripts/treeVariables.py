'''
In this file you can define the branch to be included in the flat ntuple.
You can find some basic quantities here, expand with more specific observables,
such as isolation etc...
'''
branches = [
    'run',
    'lumi',
    'event',
    'nvtx',
    'tau_reco_mass',
    'tau_reco_pt',
    'tau_reco_eta',
    'tau_reco_phi',
    'tau_reco_charge',
    'tau_reco_decaymode',
    'tau_gen_pt',
    'tau_gen_eta',
    'tau_gen_phi',
    'tau_gen_charge',
    'tau_gen_decaymode',
    'tau_reco_iso',
    'tau_gen_vis_mass',
    'tau_gen_vis_pt',
    'tau_gen_vis_signal_dR',
    'tau_gen_vis_signal_dR_pt1',
    'tau_gen_all_signal_dR_pt1',
    'tau_reco_signal_dR',
    'tau_gen_vis_eta',
    'tau_gen_vis_phi',
    'jet_mass',
    'jet_pt',
    'jet_eta',
    'jet_phi',
    'jet_charge',
    'gamma_gen_min_pt',
    'gamma_reco_min_pt',
]
#['tau_reco_byDeepTau2017v2p1VSjetraw'] +\
#['tau_reco_by%sDeepTau2017v2p1VSjet' %ii for ii in  ['VVVLoose', 'VVLoose', 'VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']]



