#!/bin/env python

import glob
import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
r.gInterpreter.AutoLoad('TLorentzVector')
r.gSystem.Load('lib/libDelphes.so')
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
r.gStyle.SetOptStat(0)
r.gStyle.SetOptTitle(0)

def main() :
    verbose = False
    inputFileNames = {'sm_zb'    : glob.glob('data/sm_zb_mm_take1/*.root'),
                      'sm_zbbar' : glob.glob('data/sm_zbbar_mm_take1/*.root'),
                      'bm_zb'    : glob.glob('data/bm_zb_mm_take1/*.root'),
                      'bm_zbbar' : glob.glob('data/bm_zbbar_mm_take1/*.root'),
                      }
    histosPerSample = dict((s, buildHistos('_'+s)) for s in inputFileNames.keys())
    for s, filenames in inputFileNames.iteritems() :
        print s
        histos = histosPerSample[s]
        tree = buildInputChain(filenames)
        tree.GetEntry(0)
        branches_muon  = get_muon_branches(tree)
        branches_jets  = get_jet_branches(tree)
        branches_truth = get_true_part_branches(tree)
        nEntries = 10000 #100
        #nEntries = tree.GetEntriesFast()
        #nEntries = 50000*10
        for iEntry in xrange(nEntries) :
            tree.GetEntry(iEntry)
            # muons    = get_muons(branches_muon)
            # jets     = get_jets (branches_jets)
            # print "muons[%d]"%len(muons)," pt: ",["%.3f"%m.Pt() for m in muons]
            # print "jets[%d]"%len(jets)," pt: ",["%.3f"%j.Pt() for j in jets]
            # if len(muons) < 3 : continue
            truePart = get_true_particles(branches_truth)
            truePart = filter(lambda p : p.pid in [+13, -13,  +5,  -5], truePart)
            muons = sortedByPt(filter(isGoodMuon, truePart))
            bbar  = sortedByPt(filter(isGoodBparton, truePart))
            if verbose :
                printPart(truePart)
                printMuons(muons)
                printBbbar(bbar)
            if len(muons)<2 or len(bbar)<1 : continue
            mu0, mu1 = muons[0], muons[1]
            sameSign = mu0.charge==mu1.charge
            m_ll, m_z = (mu0 + mu1).M(), 91.2
            if sameSign or abs(m_ll-m_z) > 15.0 : continue
            b = bbar[0] # just take the highest pt one
            mup = mu0 if mu0.charge>0. else mu1
            mum = mu0 if mu0.charge<0. else mu1
            deta = deltaEta(mu0, mu1, b)
            if verbose :
                print m_ll
                print deta
            histos['deta' ].Fill(deta)
            histos['ptmup'].Fill(mup.Pt())
            histos['ptmum'].Fill(mum.Pt())
            histos['ptb'  ].Fill(b.Pt())
            h_detaDiff = histos['deta_pos'] if m_ll>m_z else  histos['deta_neg']
            h_detaDiff.Fill(deta)
        histosPerSample[s] = histos
    print "draw histos"
    histos = histosPerSample.itervalues().next().keys()
    samples = histosPerSample.keys()
    for h in histos :
        plotHistos(h, dict([(s,histosPerSample[s][h]) for s in samples]))


def isMuon(p) : return p.pid in [+13, -13] and p.status is 1
def isGoodMuon(p) : return isMuon(p) and p.Pt()>10.0 and abs(p.Eta())<2.5
def isBparton(p) : return p.pid in [ +5,  -5] and p.status is 3 # is this right?
def isGoodBparton(p) : return isBparton(p) and p.Pt()>10.0 and abs(p.Eta())<2.5
def sortedByPt(coll) : return sorted(coll, key=lambda p : p.Pt(), reverse=True)
def deltaEta(mu0, mu1, b) :
    mu = mu0 if mu0.charge*b.charge > 0 else mu1
    return mu.Eta() - b.Eta()
def openInputFiles(fnames_dict={}) :
    return dict([(k, r.TFile.Open(f)) for k, f in fnames_dict.iteritems()])
def getInputTrees(filesDict) :
    treename = 'Delphes'
    return dict([(k, f.Get(treename)) for k, f in filesDict.iteritems()])
def buildInputChain(filenames=[], treename='Delphes') :
    chain = r.TChain(treename)
    for f in filenames : chain.Add(f)
    return chain
def get_jet_branches(tree) :
    """
    These are the ones we need:
    Jet.PT
    Jet.Eta
    Jet.Phi
    Jet.Mass
    Jet.BTag
    Jet.Charge
    """
    return getattr(tree,'Jet')
class Jet(r.TLorentzVector):
    def __init__(self, pt, eta, phi, mass, btag, charge) :
        super(Jet, self).__init__()
        self.SetPtEtaPhiM(pt, eta, phi, mass)
        self.charge = charge
        self.btag = btag
def get_jets(jet_branches) :
    return [Jet(j.PT, j.Eta, j.Phi, j.Mass, j.BTag, j.Charge) for j in jet_branches]

def get_muon_branches(tree) :
    """
    These are the ones we need:
    Muon.PT
    Muon.Eta
    Muon.Phi
    Muon.Charge
    """
    return getattr(tree,'Muon')
class Muon(r.TLorentzVector):
    def __init__(self, pt, eta, phi, charge) :
        super(Muon, self).__init__()
        m_mu=0.1057
        self.SetPtEtaPhiM(pt, eta, phi, m_mu)
        self.charge = charge
def get_muons(muon_branches) :
    return [Muon(m.PT, m.Eta, m.Phi, m.Charge) for m in muon_branches]

def get_true_part_branches(tree) :
    """
    These are the ones we need:
    Particle.PID
    Particle.Status
    Particle.Charge
    Particle.Mass
    Particle.E
    Particle.Px
    Particle.Py
    Particle.Pz
    Particle.PT
    Particle.Eta
    Particle.Phi
    """
    return getattr(tree,'Particle')
class Particle(r.TLorentzVector):
    def __init__(self, pid, status, charge, pt, eta, phi, energy) :
        super(Particle, self).__init__()
        self.SetPtEtaPhiE(pt, eta, phi, energy)
        self.pid    = pid
        self.status = status
        self.charge = charge
def get_true_particles(true_part_branches) :
    return [Particle(p.PID, p.Status, p.Charge, p.PT, p.Eta, p.Phi, p.E)
            for p in true_part_branches]


def buildHistos(suffix='') :
    histos= { 'deta' : r.TH1F('h_deltaEta'+suffix, ';#Delta#eta', 64, -6.4, +6.4),
              'deta_pos' : r.TH1F('h_deltaEtaPos'+suffix, 'm_{ll}>m_{Z};#Delta#eta', 64, -6.4, +6.4),
              'deta_neg' : r.TH1F('h_deltaEtaNeg'+suffix, 'm_{ll}<m_{Z};#Delta#eta', 64, -6.4, +6.4),
              'ptmup': r.TH1F('h_pt_mu_p'+suffix, ';p_{T,#mu+}', 100, 0.0, 500.0),
              'ptmum': r.TH1F('h_pt_mu_m'+suffix, ';p_{T,#mu-}', 100, 0.0, 500.0),
              'ptb'  : r.TH1F('h_pt_b'+suffix,    ';p_{T,b}', 100, 0.0, 500.0),
              }
    for h in histos.values() : h.SetDirectory(0)
    return histos
def printPart(part) :
    print "part[%d]"%len(truePart)
    print "     pt : ",["%.3f"%p.Pt() for p in truePart]
    print "     pid: ",[p.pid for p in truePart]
def printMuons(muons) :
    print "muons[%d]"%len(muons)
    print "     pt : ",["%.3f"%m.Pt() for m in muons]
    print "     eta: ",["%.3f"%m.Eta() for m in muons]
    print "     phi: ",["%.3f"%m.Phi() for m in muons]
    print "  status: ",["%.3f"%m.status for m in muons]
def printBbbar(bbar) :
    print "bbar[%d]"%len(bbar)
    print "      pt: ",["%.3f"%b.Pt() for b in bbar]
    print "     eta: ",["%.3f"%b.Eta()  for b in bbar]
    print "     phi: ",["%.3f"%b.Phi()  for b in bbar]
    print "  status: ",[" %d"%b.status for b in bbar]
def computeAsymm(h) :
    nBins = h.GetNbinsX()
    bins = range(1, 1+nBins)
    binCenters = [h.GetBinCenter(b) for b in bins]
    posBins = [b for b,c in zip(bins, binCenters) if c>0]
    negBins = [b for b,c in zip(bins, binCenters) if c<0]
    nPosBins, nNegBins = len(posBins), len(negBins)
    assert nBins==(nPosBins+nNegBins),"N_b+ + N_b- = %d + %d !=%d"%(nNegBins, nPosBins, nBins)
    assert nPosBins==nNegBins,"N_b+ != N_b- (%d!=%d)"%(nNegBins, nPosBins)
    print posBins
    print negBins
    nP = sum([h.GetBinContent(b) for b in posBins])
    nN = sum([h.GetBinContent(b) for b in negBins])
    return float(nP-nN)/float(nP+nN) if nP or nN else 0.0
colors = {'sm_zb':r.kBlue,
          'sm_zbbar':r.kRed,
          'bm_zb':r.kOrange,
          'bm_zbbar':r.kMagenta}
def plotHistos(histoname='', histosDict={}, colors=colors) :
    c = r.TCanvas('c_'+histoname)
    c.cd()
    pm = histosDict.itervalues().next()
    pm.Draw('axis')
    pm.SetMaximum(1.1*max([h.GetMaximum() for h in histosDict.values()]))
    leg = r.TLegend(0.65, 0.8, 0.9, 0.9)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(1)
    for s,h in histosDict.iteritems() :
        print s
        print histoname
        h.SetLineColor(colors[s])
        h.SetMarkerColor(colors[s])
        h.SetLineWidth(2*h.GetLineWidth())
        h.Draw('same')
        label = s
        label += " (N=%d"%h.GetEntries()
        label += ", "
        label += "<v>=%.3f"%h.GetMean() if histoname!='deta' else "A=%.3f"%computeAsymm(h)
        label += ")"
        leg.AddEntry(h, label, 'l')
    leg.Draw()
    c.SaveAs(c.GetName()+'.png')

if __name__=='__main__' :
    main()
