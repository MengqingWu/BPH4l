#!/usr/bin/env python
## ---------------------------------------------
## Author: Mengqing Wu @ Nov 2016
## 0) You have a event selection
##    looping over events stored in a muon ntuple/tree produced by Heppy 
## 1) Histograms will be stored in 'outMuAna.root'
## 2) Events at each cut level appended to 'cutflow_out.txt'
## 3) The sequence is below:
##   -- Global Setup
##   -- init Histograms
##   -- define Cuts
##   -- Loop
## ---------------------------------------------

import os, sys, math, numpy
from itertools import combinations
from ROOT import *
#from python.Pair import *
from python.Particle import Particle
from python.Quad import Quad
from Cuts import *

#---- global setup:
dotest = False
mode = 'all'
indir = "../data"
sample = "MuOnia_2016ICHEP.root"
treename = 'tree'
chainname = [indir,sample,treename]
fchain = TChain("tree")

outdir='./'
fout=TFile(outdir+"outScanner.root","recreate")
outtxt = open('cutflow_scanner.txt', 'a')

status = fchain.Add('/'.join(chainname), 0)
if status<=0: raise RuntimeError, "[cutflow.py] ERROR! Invalid tree added to fchain, please check!"
else: print "== Successfully load root file =="


#---- init histos: 
cutflow = TH1F("h_cutflow", "cutflow; cut# ; / event ",10 , 0.,10. )

h_m2mu = TH1F("h_m2mu", "soft muons w/o overlap; m^{2#mu} [GeV]; / 25MeV", 480, 0., 12.)

#--- Cutflow txt output:
nhlt = nmuId = n4mu = 0
nMu4Cand_vtx = nMu4Cand_vtxMass = ngoodMu4 = ntightMu4 = 0

#--- Loop:
print fchain.GetEntriesFast(),' entries to process:'
outtxt.write('=> we have '+ str(fchain.GetEntriesFast())+ ' entries.\n')

maxEntries=5001 if dotest else fchain.GetEntriesFast()


for ientry in range(0, maxEntries):
        fchain.GetEntry(ientry)
        if ientry%50000==0: print "Entry ",ientry
        cutflow.Fill(0.)

        ##--> CUT: at least 1 neutral 4-mu combinations: FIXME!
        ##    (bad 4mu vtx w/ very default -99 chi2 prob filtered at ntuplisation level)
        if fchain.nmu4>0: cutflow.Fill(1.)
        else: continue
        
        ##--> CUT: pass the HLT
        if fchain.HLT_jpsi2mu+fchain.HLT_upsilon2mu>0: cutflow.Fill(2.)
        else: continue

        indexBadMuons = [fchain.mu_index[imu] for imu in range(fchain.nmu) if fchain.mu_isOverlap[imu]]
        indexTightMus = [fchain.mu_index[imu] for imu in range(fchain.nmu) if fchain.mu_tightMuonId[imu]]

        ##-> loop over all the 4mu combinations to get a di-mu collection:
        _all2mupairs=[]
        all2mupairs=[]
        #print "[Event level] we have %d quad-mu combinations, %d muons" % (fchain.nmu4, fchain.nmu)
        for imu4 in range(fchain.nmu4):
                i2mupairs = GetPairsIn4Mu(fchain, imu4)
                _all2mupairs+=[ i2mupairs[i] for i in i2mupairs]

        ##-> need to remove overlapped dimu (as you collect all dimus from all possible 4mu combinations)
        all2mupairs=list(_all2mupairs) # modification on all2mupairs will not place on _all2mupairs
        #print "[info] %d di-mu pairs" % (len(all2mupairs))
        for a2mu, b2mu in combinations(_all2mupairs, 2):
                if set(a2mu.index)==set(b2mu.index):
                        if a2mu in all2mupairs and b2mu in all2mupairs: all2mupairs.remove(b2mu)
                        #print "[info] we have duplicate 2mu pairs: a(%s) and b(%s)" % (a2mu.index, b2mu.index)
                else: pass
        #print "[info] %d di-mu pairs" % (len(all2mupairs))
        ##-> double check if overlapped dimu pairs still exist:
        for x,y in combinations(all2mupairs,2):
                if set(x.index)==set(y.index):
                        print "[Warning] we still have duplicate 2mu pairs: a(%s) and b(%s)" % (x.index, y.index)
                else: pass

        ##-> CUT: muon overlap removal
        good2mupairCands=list(all2mupairs)
        for ipair in all2mupairs:
                if filter( lambda x: x in indexBadMuons, ipair.index): good2mupairCands.remove(ipair)
                else: h_m2mu.Fill(ipair.mass())
                        
        ##-> CUT: at least 1 tight muon:
        good2mupairs=[]
        for icand in good2mupairCands:
                if len(filter( lambda x: x in indexTightMus, ipair.index))>=1: good2mupairs+=[ipair]
                else: continue
        
                
# outtxt.write(' raw events:     '+ str(maxEntries)+'\n')
# outtxt.write(' mu4 vtx cuts:   '+ str(nMu4Cand_vtx)+'\n')
# outtxt.write(' 1Jpsi+1subY:    '+ str(nMu4Cand_vtxMass)+'\n')
# outtxt.write(' no overlap mu:  '+ str(ngoodMu4)+'\n')
# outtxt.write(' >=2 tight mu:   '+ str(ntightMu4)+'\n')
# outtxt.write(' dimu HLT:       '+ str(nhlt)+'\n')
# outtxt.write(' mu ID+Kin :     '+ str(nmuId)+'\n')
# outtxt.write(' 0-charge 4mu:   '+ str(n4mu)+'\n')

outtxt.close()

#---- Finalize:
fout.cd()
fout.Write()
fout.Close()

