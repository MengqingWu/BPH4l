#!/usr/bin/env python
import os, sys, math, numpy
from  itertools import combinations, product
from ROOT import *
from python.Pair import *
from python.Particle import *

dotest = False
indir = "../data"
sample = "MuOnia_2016ICHEP.root"
treename = 'tree'
chainname = [indir,sample,treename]
fchain = TChain("tree")

outdir='./'
fout=TFile(outdir+"outMuAna.root","recreate")
outtxt = open('cutflow_out.txt', 'a')

status = fchain.Add('/'.join(chainname), 0)
if status<=0: raise RuntimeError, "[cutflow.py] ERROR! Invalid tree added to fchain, please check!"
else: print "== Successfully load root file =="


#---- init histos: 
cutflow = TH1F("h_cutflow", "cutflow; cut# ; / event ",10 , 0.,10. )
#cutflow_w = TH1F("h_cutflow_w", "weighted cutflow; cut# ; / yields ",10 , 0.,10. ) # only with MC
h_mu4 = TH1F("h_mu4", "# of 0-charge 4mu", 10, 0, 10)
h_mu = TH1F("h_mu", "# of mu passing ID", 10, 0, 10)
h_mu_pt = TH1F("h_mu_pt", "evts with 0-charged 4mu (id+kin); p_{T}^{#mu}; / 1GeV", 100, 0, 100)
h_mu_eta = TH1F("h_mu_eta", "evts with 0-charged 4mu (id+kin); |#eta^{#mu}|; / 0.05", 50, 0, 2.5)
h_mu_dxy = TH1F("h_mu_dxy", "evts with 0-charged 4mu (id+kin); |dxy^{#mu}|; / 0.01", 200, -1, 1)
h_mu_dz = TH1F("h_mu_dz", "evts with 0-charged 4mu (id+kin); |dz^{#mu}|; / 0.2", 200, -20, 20)

h_mu4_mass = TH1F("h_mu4_mass","invariant mass of best 4mu combination; m^{4#mu} [GeV]; / 200MeV", 120, 0., 24.)

#--- Cuts:
nhlt = nmuId = n4mu = 0

        
#--- Loop:
print fchain.GetEntriesFast(),' entries to process:'
outtxt.write('=> we have '+ str(fchain.GetEntriesFast())+ ' entries.\n')

maxEntries=5001 if dotest else fchain.GetEntriesFast()

for ientry in range(0, maxEntries):
        fchain.GetEntry(ientry)
        if ientry%50000==0: print "Entry ",ientry
        cutflow.Fill(0.)
        ##--> pass the HLT:
        if fchain.HLT_jpsi2mu+fchain.HLT_upsilon2mu>0:
                nhlt+=1; cutflow.Fill(1.)

                ##--> counting muons with a set of criteria
                goodMuBox = []
                mu4Box = []
                for ilep in range(fchain.nlep):
                        ##--> good muons: ID and Kin cuts
                        if abs(fchain.lep_pdgId[ilep])==13 and fchain.lep_preMuonId[ilep] and fchain.lep_softMuonId[ilep] and fchain.lep_pt[ilep]>=2.0 and abs(fchain.lep_eta[ilep])<=2.4:
                                goodMuBox.append(ilep)
                                ##--> get all 0-charged 4-mu combinations:
                                for l1,l2,l3,l4 in combinations(goodMuBox, 4):
                                        if fchain.lep_charge[l1]+fchain.lep_charge[l2]+fchain.lep_charge[l3]+fchain.lep_charge[l4]==0:
                                                mu4Box.append([l1,l2,l3,l4])
                # # di-mu pairing :
                # for quatreMu in mu4Box:
                #         mupos=[]
                #         muneg=[]
                #         comb =[]
                #         for imu in quatreMu:
                #                 if fchain.lep_charge[imu]>0: mupos.append(imu)
                #                 elif fchain.lep_charge[imu]<0: muneg.append(imu)
                #                 else: print "[Warning] This should be really weird: mu (%d) has charge = %d " % (imu, fchain.lep_charge[imu])
                #                 for mup, mun in product(mupos, muneg):
                #                         comb.append(Pair(fchain, mup, mun))
                
                ##--> apply cuts:
                if len(goodMuBox)>=4:
                        nmuId+=1; cutflow.Fill(2.)
                        h_mu.Fill(len(goodMuBox))
                        h_mu4.Fill(len(mu4Box))

                        ##--> pass 0-charge 4mu:
                        if len(mu4Box)>=1:
                                n4mu+=1; cutflow.Fill(3.)
                                for igoodMu in goodMuBox:
                                        h_mu_pt.Fill(fchain.lep_pt[igoodMu])
                                        h_mu_eta.Fill(abs(fchain.lep_eta[igoodMu]))
                                        h_mu_dxy.Fill(fchain.lep_dxy[igoodMu])
                                        h_mu_dz.Fill(fchain.lep_dz[igoodMu])
                                        
                                ##--> choose a 4mu comb based on dz:
                                mindz = 999
                                best4mu = None
                                for imu4 in mu4Box:
                                        sumdz  = 0
                                        for imu in imu4:
                                                sumdz += fchain.lep_dz[imu]
                                                #sumdxy += fchain.lep_dxy[imu]
                                        if sumdz<mindz: mindz=sumdz; best4mu = imu4
                                        #mindxy = min(sumdxy, mindxy)
                                if best4mu:
                                        ##--> fairly rare if no best4mu in this case (means you have all 4mu with sum(dz)>999!)
                                        best4muObj = [Particle(fchain, ibestmu) for ibestmu in best4mu]
                                        bestCombLV = best4muObj[0].p4()+best4muObj[1].p4()+best4muObj[2].p4()+best4muObj[3].p4()
                                        
                                        h_mu4_mass.Fill(bestCombLV.M())
                                        
                                
outtxt.write(' raw events:     '+ str(maxEntries)+'\n')
outtxt.write(' dimu HLT:       '+ str(nhlt)+'\n')
outtxt.write(' mu ID+Kin :     '+ str(nmuId)+'\n')
outtxt.write(' 0-charge 4mu:   '+ str(n4mu)+'\n')

outtxt.close()
#---- Finalize:
fout.cd()
fout.Write()
fout.Close()
