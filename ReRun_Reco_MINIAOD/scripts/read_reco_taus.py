import ROOT
import argparse
from DataFormats.FWLite import Events, Handle
from array import array
from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument('--file' , required = True)
parser.add_argument('--label', default  = 'test')
parser.add_argument('--max'  , default  = None)

args = parser.parse_args()

outfile = ROOT.TFile.Open('tau_reco_{LAB}.root'.format(LAB = args.label), 'RECREATE')

to_plot = OrderedDict(
    {
        'DM'        : -99   ,
        'newDMF'    : -99   ,
        'DMF'       : -99   ,
        'pt'        : -99   ,
    }
)

print 'INFO: reading file', args.file 
events = Events(args.file)

label_taus  = ('selectedPatTaus')
handle_taus = Handle('std::vector<pat::Tau>')

ntuple = ROOT.TNtuple('tree', 'tree', ':'.join(to_plot.keys()))

MAX = args.max
TOT = events.size()
for i, ev in enumerate(events): 
    if i%100==0:
        print '===> processing %d / %d event' %(i, TOT)

    ev.getByLabel(label_taus, handle_taus)
    taus = handle_taus.product()
    taus = [tau for tau in taus if tau.pt()>18.]

    for tt in taus:
        to_plot['DM'        ] = tt.decayMode()
        to_plot['newDMF'    ] = tt.tauID('decayModeFindingNewDMs')
        to_plot['DMF'       ] = tt.tauID('decayModeFinding')
        to_plot['pt'        ] = tt.pt()

    ntuple.Fill(array('f', to_plot.values()))

outfile.cd()
ntuple.Write()
outfile.Close()
    