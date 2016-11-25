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
#h_mu4 = TH1F("h_mu4", "# of 0-charge 4mu", 10, 0, 10)
h_mu_pt = TH1F("h_mu_pt", "evts with 0-charged 4mu (id+kin); p_{T}^{#mu}; / 1GeV", 100, 0, 100)
h_mu_eta = TH1F("h_mu_eta", "evts with 0-charged 4mu (id+kin); |#eta^{#mu}|; / 0.05", 50, 0, 2.5)
h_mu_dxy = TH1F("h_mu_dxy", "evts with 0-charged 4mu (id+kin); |dxy^{#mu}|; / 0.01", 200, -1, 1)
h_mu_dz = TH1F("h_mu_dz", "evts with 0-charged 4mu (id+kin); |dz^{#mu}|; / 0.2", 200, -20, 20)


h_mu = TH1F("h_mu", "# of mu passing ID", 10, 0, 10)
h_mu2_1a_mass = TH1F("h_mu2_1a_mass","invariant mass of pair1a; m^{#mu#mu} [GeV]; / 200MeV", 100, 0., 20.)
h_mu2_1a_fitmass = TH1F("h_mu2_1a_fitmass","invariant mass of dimu '1a'; m^{#mu#mu} [GeV]; / 200MeV", 100, 0., 20.)
h_mu2_1a_fitmasserr = TH1F("h_mu2_1a_fitmasserr","invariant mass of dimu '1a'; m^{#mu#mu} [GeV]; / 2.5MeV", 100, 0, 0.25)

h_mu4_mass = TH1F("h_mu4_mass","invariant mass of best 4mu combination; m^{4#mu} [GeV]; / 200MeV", 120, 0., 24.)
h_mu4_pt = TH1F("h_mu4_pt","pT of best 4mu combination; p_{T}^{4#mu} [GeV]; / 500MeV", 180,0,90)
#--- Cuts:
nhlt = nmuId = n4mu = 0
nMu4Cand_vtx = nMu4Cand_vtxMass = ngoodMu4 = ntightMu4 = 0
jpsi_mass = 3.096916 # pdg mass in GeV
upsilon_mass = 9.46030 # pdg mass in GeV
massErr_sf = 1.105 # see AN2013_099_v17: derived from a fit to the mass pull distribution

min2muVtx=0.005
min4muVtx=0.05

quadMuVtxCut = (lambda tree,index: tree.mu4_quad_vtxProb[index]>min4muVtx )
# at least 1 pair of diMus both w/diMuVtxProb>0.005
diMuVtxCut = (lambda tree,index: (tree.mu4_1a_vtxProb[index]>min4muVtx and tree.mu4_1b_vtxProb[index]>min4muVtx) or\
              (tree.mu4_2a_vtxProb[index]>min4muVtx and tree.mu4_2b_vtxProb[index]>min4muVtx))

# pair = [(m, m_sigma),(m, m_sigma)], ordered in mass (high to low)
m_JpsiSubYCut = (lambda pair: abs(pair[0][0]-jpsi_mass)<3*pair[0][1] and\
                 abs(pair[1][0]-upsilon_mass)>3*pair[1][1])

def passJpsiSubYCut(tree, mu4_index):
        if fchain.mu4_1a_fitMassErr2[mu4_index]<0 or fchain.mu4_1b_fitMassErr2[mu4_index]<0 or fchain.mu4_2a_fitMassErr2[mu4_index]<0 or fchain.mu4_2b_fitMassErr2[mu4_index]<0:
                return False
        _1a=( fchain.mu4_1a_fitMass[mu4_index],  math.sqrt(fchain.mu4_1a_fitMassErr2[mu4_index])*massErr_sf )
        _1b=( fchain.mu4_1b_fitMass[mu4_index],  math.sqrt(fchain.mu4_1b_fitMassErr2[mu4_index])*massErr_sf )
        _2a=( fchain.mu4_2a_fitMass[mu4_index],  math.sqrt(fchain.mu4_2a_fitMassErr2[mu4_index])*massErr_sf )
        _2b=( fchain.mu4_2b_fitMass[mu4_index],  math.sqrt(fchain.mu4_2b_fitMassErr2[mu4_index])*massErr_sf )
        pair1=[_1a, _1b]
        pair2=[_2a, _2b]
        pair1.sort(key = lambda x: x[0], reverse = True)
        pair2.sort(key = lambda x: x[0], reverse = True)
        #print "[output] pair1 mass(%.4f, %.4f)" % (pair1[0][0], pair1[1][0])
        #print "[output] pair2 mass(%.4f, %.4f)" % (pair2[0][0], pair2[1][0])
        if m_JpsiSubYCut(pair1) or m_JpsiSubYCut(pair2):
                return True
        else:   return False


#--- Loop:
print fchain.GetEntriesFast(),' entries to process:'
outtxt.write('=> we have '+ str(fchain.GetEntriesFast())+ ' entries.\n')

maxEntries=5001 if dotest else fchain.GetEntriesFast()

for ientry in range(0, maxEntries):
        fchain.GetEntry(ientry)
        if ientry%50000==0: print "Entry ",ientry
        cutflow.Fill(0.)

        ##--> 4mu selections to select one 4mu:
        ## fchain.nmu4>0: skip those bad 4mu vertex w/ vtxProb==-99 (default)
        if fchain.nmu4>0:
                cutflow.Fill(1.)
                # good 4mu/2mu vertices selection:
                mu4Cand_vtx = [ imu4 for imu4 in range(fchain.nmu4) if quadMuVtxCut(fchain,imu4) and diMuVtxCut(fchain, imu4) ]
        else: continue

        if len(mu4Cand_vtx)>0:
                nMu4Cand_vtx+=1;cutflow.Fill(2.)
                #print "[output] evt = %d, nMu4Cand_vtx = %d"%(fchain.evt, len(mu4Cand_vtx))
                mu4Cand_vtxMass = [ imu4 for imu4 in mu4Cand_vtx if passJpsiSubYCut(fchain, imu4) ]
        else: continue

        if len(mu4Cand_vtxMass)>0:
                nMu4Cand_vtxMass+=1;cutflow.Fill(3.)
                ##--> algo to remove overlap mu before best mu4 decided:
                # indexBadMuons = [fchain.mu_index[imu] for imu in range(fchain.nmu) if fchain.mu_isOverlap]
                # goodMu4 = [ imu4 for imu4 in mu4Cand_vtxMass if\
                #             filter( lambda x: x in indexBadMuons,\
                #                     [fchain.mu4_quad_l1_index[imu4],fchain.mu4_quad_l2_index[imu4],\
                #                      fchain.mu4_quad_l3_index[imu4],fchain.mu4_quad_l3_index[imu4] ])]
                # if len(goodMu4)>0:
                #         ngoodMu4+=1
                #         print "[output] evt = %d, nMu4Cand_vtx = %d, nMu4Cand_vtxMass = %d, ngoodMu4 = %d"\
                #                 %(fchain.evt, len(mu4Cand_vtx), len(mu4Cand_vtxMass), len(goodMu4))
                #-->
                
                mu4Cand_dict = {imu4:fchain.mu4_quad_vtxProb[imu4] for imu4 in mu4Cand_vtxMass}
                bestMu4 = min(mu4Cand_dict, key=mu4Cand_dict.get)
        else: continue

        indexBadMuons = [fchain.mu_index[imu] for imu in range(fchain.nmu) if fchain.mu_isOverlap]
        indexTightMus = [fchain.mu_index[imu] for imu in range(fchain.nmu) if fchain.mu_tightMuonId]
        indexBest4Mus = [fchain.mu4_quad_l1_index[bestMu4],fchain.mu4_quad_l2_index[bestMu4],fchain.mu4_quad_l3_index[bestMu4],fchain.mu4_quad_l3_index[bestMu4]]
        ## muon overlap removal:
        if filter( lambda x: x in indexBadMuons, indexBest4Mus):
                ngoodMu4+=1; cutflow.Fill(4.)
        else: continue
        ## at least 2 tight muons:
        if len(filter( lambda x: x in indexTightMus, indexBest4Mus))>=2 :
                ntightMu4+=1; cutflow.Fill(5.)
        else: continue

        ##--> pass the HLT: 
        if fchain.HLT_jpsi2mu+fchain.HLT_upsilon2mu>0:
                nhlt+=1; cutflow.Fill(6.)
                h_mu4_mass.Fill(fchain.mu4_quad_mass[bestMu4])
                h_mu4_pt.Fill(fchain.mu4_quad_pt[bestMu4])
                h_mu2_1a_mass.Fill(fchain.mu4_1a_mass[bestMu4])
                h_mu2_1a_fitmass.Fill(fchain.mu4_1a_fitMass[bestMu4])
                h_mu2_1a_fitmasserr.Fill(fchain.mu4_1a_fitMassErr2[bestMu4])
        else: continue
                # ##--> counting muons with a set of criteria
                # goodMuBox = []
                # mu4Box = []
                # for imu in range(fchain.nmu):
                #         ##--> good muons: ID and Kin cuts
                #         if abs(fchain.mu_pdgId[imu])==13 and fchain.mu_preMuonId[imu] and fchain.mu_softMuonId[imu] and fchain.mu_pt[imu]>=2.0 and abs(fchain.mu_eta[imu])<=2.4:
                #                 goodMuBox.append(imu)
                #                 ##--> get all 0-charged 4-mu combinations:
                #                 for l1,l2,l3,l4 in combinations(goodMuBox, 4):
                #                         if fchain.mu_charge[l1]+fchain.mu_charge[l2]+fchain.mu_charge[l3]+fchain.mu_charge[l4]==0:
                #                                 mu4Box.append([l1,l2,l3,l4])
   
                # ##--> apply cuts:
                # if len(goodMuBox)>=4:
                #         nmuId+=1; cutflow.Fill(2.)
                #         h_mu.Fill(len(goodMuBox))
                #         h_mu4.Fill(len(mu4Box))

                #         ##--> pass 0-charge 4mu:
                #         if len(mu4Box)>=1:
                #                 n4mu+=1; cutflow.Fill(3.)
                #                 for igoodMu in goodMuBox:
                #                         h_mu_pt.Fill(fchain.mu_pt[igoodMu])
                #                         h_mu_eta.Fill(abs(fchain.mu_eta[igoodMu]))
                #                         h_mu_dxy.Fill(fchain.mu_dxy[igoodMu])
                #                         h_mu_dz.Fill(fchain.mu_dz[igoodMu])
                                        
                #                 ##--> choose a 4mu comb based on dz:
                #                 mindz = 999
                #                 best4mu = None
                #                 for imu4 in mu4Box:
                #                         sumdz  = 0
                #                         for imu in imu4:
                #                                 sumdz += fchain.mu_dz[imu]
                #                                 #sumdxy += fchain.mu_dxy[imu]
                #                         if sumdz<mindz: mindz=sumdz; best4mu = imu4
                #                         #mindxy = min(sumdxy, mindxy)
                #                 if best4mu:
                #                         ##--> fairly rare if no best4mu in this case (means you have all 4mu with sum(dz)>999!)
                #                         best4muObj = [Particle(fchain, ibestmu) for ibestmu in best4mu]
                #                         bestCombLV = best4muObj[0].p4()+best4muObj[1].p4()+best4muObj[2].p4()+best4muObj[3].p4()
                                        
                #                         h_mu4_mass.Fill(bestCombLV.M())
                                        
                                
outtxt.write(' raw events:     '+ str(maxEntries)+'\n')
outtxt.write(' mu4 vtx cuts:   '+ str(nMu4Cand_vtx)+'\n')
outtxt.write(' 1Jpsi+1subY:    '+ str(nMu4Cand_vtxMass)+'\n')
outtxt.write(' no overlap mu:  '+ str(ngoodMu4)+'\n')
outtxt.write(' >=2 tight mu:   '+ str(ntightMu4)+'\n')
outtxt.write(' dimu HLT:       '+ str(nhlt)+'\n')
#outtxt.write(' mu ID+Kin :     '+ str(nmuId)+'\n')
#outtxt.write(' 0-charge 4mu:   '+ str(n4mu)+'\n')

outtxt.close()
#---- Finalize:
fout.cd()
fout.Write()
fout.Close()
