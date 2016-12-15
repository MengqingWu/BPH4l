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
from itertools import combinations, product
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
fout=TFile(outdir+"outMuAna.root","recreate")
outtxt = open('cutflow_out.txt', 'a')

status = fchain.Add('/'.join(chainname), 0)
if status<=0: raise RuntimeError, "[cutflow.py] ERROR! Invalid tree added to fchain, please check!"
else: print "== Successfully load root file =="


#---- init histos: 
cutflow = TH1F("h_cutflow", "cutflow; cut# ; / event ",10 , 0.,10. )
#cutflow_w = TH1F("h_cutflow_w", "weighted cutflow; cut# ; / yields ",10 , 0.,10. ) # only with MC
#h_mu4 = TH1F("h_mu4", "# of 0-charge 4mu", 10, 0, 10)
#h_mu_pt = TH1F("h_mu_pt", "evts with 0-charged 4mu (id+kin); p_{T}^{#mu}; / 1GeV", 100, 0, 100)
#h_mu_eta = TH1F("h_mu_eta", "evts with 0-charged 4mu (id+kin); |#eta^{#mu}|; / 0.05", 50, 0, 2.5)
#h_mu_dxy = TH1F("h_mu_dxy", "evts with 0-charged 4mu (id+kin); |dxy^{#mu}|; / 0.01", 200, -1, 1)
#h_mu_dz = TH1F("h_mu_dz", "evts with 0-charged 4mu (id+kin); |dz^{#mu}|; / 0.2", 200, -20, 20)

h_mu = TH1F("h_mu", "# of mu passing ID", 10, 0, 10)
h_mu2_1a_mass = TH1F("h_mu2_1a_mass","invariant mass of pair1a; m^{#mu#mu} [GeV]; / 200MeV", 100, 0., 20.)
h_mu2_1a_fitmass = TH1F("h_mu2_1a_fitmass","invariant mass of dimu '1a'; m^{#mu#mu} [GeV]; / 200MeV", 100, 0., 20.)
h_mu2_1a_fitmasserr = TH1F("h_mu2_1a_fitmasserr","invariant mass of dimu '1a'; m^{#mu#mu} [GeV]; / 2.5MeV", 100, 0, 0.25)

h_mu4_ndiMuPass = TH1F("h_mu4_ndiMuPass", "# of di-mu pairing in the best 4mu; ; events", 3, -0.5, 2.5)
h_mu4_mass = TH1F("h_mu4_mass","invariant mass of best 4mu combination; m^{4#mu} [GeV]; / 200MeV", 120, 0., 24.)
h_mu4_pt   = TH1F("h_mu4_pt","pT of best 4mu combination; p_{T}^{4#mu} [GeV]; / 500MeV", 180,0,90)

h_jpsi_mass = TH1F("h_jpsi_mass","J#Psi mass; m^{J#Psi} [GeV]; / 25MeV", 200, 0, 5)
h_jpsi_pt   = TH1F("h_jpsi_pt","J#Psi pT; p_{T}^{J#Psi} [GeV]; / 500MeV", 100, 0, 50)

h_Y_mass = TH1F("h_Y_mass","Y mass; m^{#Upsilon} [GeV]; / 25MeV", 480, 0., 12.)
h_Y_pt   = TH1F("h_Y_pt","Y pT; p_{T}^{Y} [GeV]; / 500MeV", 100, 0, 50)

h2_mjpsi_mY   = TH2F("h2_mjpsi_mY",   " TLorentzMass;  m^{J#Psi} [GeV];  m^{#Upsilon} [GeV]", 200, 0, 5,  480, 0., 12.)
h2_mjpsi_m4mu = TH2F("h2_mjpsi_m4mu", " TLorentzMass;  m^{J#Psi} [GeV];  m^{4#mu} [GeV]"   , 200, 0, 5,  120, 0., 24.)
h2_mY_m4mu    = TH2F("h2_mY_m4mu ",   " TLorentzMass;  m^{#Upsilon} [GeV];  m^{4#mu} [GeV]", 480, 0, 12, 120, 0., 24.)

## --> check the pt and mass of 4mu when it has 2 di-mu combinations passing the final cuts: <---##
h_mu4_pt_2pair = TH1F("h_mu4_pt_2pair","pT of best 4mu combination (npass==2); p_{T}^{4#mu} [GeV]; / 500MeV", 180,0,90)
h_mu4_mass_2pair = TH1F("h_mu4_mass_2pair","invariant mass of best 4mu combination (npass==2); m^{4#mu} [GeV]; / 200MeV", 120, 0., 24.)

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
        
        ##--> 4mu selections to select one 4mu:
        ## fchain.nmu4>0: skip those bad 4mu vertex w/ vtxProb==-99 (default)
        if fchain.nmu4>0:
                cutflow.Fill(1.)
                # good 4mu/2mu vertices selection:
                #mu4Cand_vtx = [ imu4 for imu4 in range(fchain.nmu4) if quadMuVtxCut(fchain,imu4) and diMuVtxCut(fchain, imu4) ]
                mu4Cand_vtx = [ imu4 for imu4 in range(fchain.nmu4) if quadMuVtxCut(fchain,imu4) and diMuVtxCut(fchain, imu4) ]
        else: continue

        if len(mu4Cand_vtx)>0:
                nMu4Cand_vtx+=1;cutflow.Fill(2.)
                #print "[output] evt = %d, nMu4Cand_vtx = %d"%(fchain.evt, len(mu4Cand_vtx))
                mu4CandObj_vtxMass=[]
                for jmu4 in mu4Cand_vtx:
                        diMuPass = quadMuMassCut().run(fchain, jmu4, mode=mode)
                        #print "[debug] evt = %d, imu4 = %d, diMuPass type = '%s' while '%s'" % (fchain.evt, jmu4, type(diMuPass), type( passJpsiSubYCut(fchain, jmu4)))
                        if len(diMuPass):
                                jmu4CandObj=Quad(fchain, jmu4)
                                setattr(jmu4CandObj, 'ndiMuPass', len(diMuPass))
                                setattr(jmu4CandObj, 'diMuPass', diMuPass)
                                mu4CandObj_vtxMass.append(jmu4CandObj)
                        else: continue
        else: continue

        if len(mu4CandObj_vtxMass)>0:
                nMu4Cand_vtxMass+=1;cutflow.Fill(3.)

                #mu4Cand_dict = {imu4:fchain.mu4_quad_vtxProb[imu4] for imu4 in mu4CandObj_vtxMass}
                #bestMu4 = min(mu4Cand_dict, key=mu4Cand_dict.get) # this is a bug: you want the max not the min! @2016Dec1
                
                mu4CandObj_vtxMass.sort(key = lambda mu4: mu4.vtxProb, reverse = True)
                bestMu4=mu4CandObj_vtxMass[0].tIndex
                bestMu4Obj = mu4CandObj_vtxMass[0]
                
        else: continue
        
        indexBadMuons = [fchain.mu_index[imu] for imu in range(fchain.nmu) if fchain.mu_isOverlap[imu]]
        indexTightMus = [fchain.mu_index[imu] for imu in range(fchain.nmu) if fchain.mu_tightMuonId[imu]]
 
        indexBest4Mus = [fchain.mu4_quad_l1_index[bestMu4],
                         fchain.mu4_quad_l2_index[bestMu4],
                         fchain.mu4_quad_l3_index[bestMu4],
                         fchain.mu4_quad_l4_index[bestMu4]]
        ## muon overlap removal:
        if filter( lambda x: x in indexBadMuons, indexBest4Mus): continue
        else: ngoodMu4+=1; cutflow.Fill(4.)
        
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
                h_mu4_ndiMuPass.Fill(bestMu4Obj.ndiMuPass)
                
                for idiMuPass in bestMu4Obj.diMuPass:
                        #print "fill the jpsi/Y: %.3f, %.3f" % (idiMuPass[0].fitMass, idiMuPass[1].fitMass)

                        mY=mJpsi=-99
                        for idiMu in idiMuPass:
                                if 'Y' in idiMu.tag():
                                        h_Y_mass.Fill(idiMu.fitMass)
                                        h_Y_pt.Fill(idiMu.pt())
                                        mY=idiMu.mass()
                                elif 'Jpsi' in idiMu.tag():
                                        h_jpsi_mass.Fill(idiMu.fitMass)
                                        h_jpsi_pt.Fill(idiMu.pt())
                                        mJpsi=idiMu.mass()
                                else: pass
                                
                                h2_mjpsi_mY.Fill(mJpsi, mY)   
                                h2_mY_m4mu.Fill(mY, bestMu4Obj.mass())
                                h2_mjpsi_m4mu.Fill(mJpsi, bestMu4Obj.mass())
                        

                if bestMu4Obj.ndiMuPass==2:
                        h_mu4_pt_2pair.Fill(bestMu4Obj.pt())
                        h_mu4_mass_2pair.Fill(bestMu4Obj.mass())
                        
        else: continue
                                      
                                
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



                ##--> algo to remove overlap mu before best mu4 decided:
                # indexBadMuons = [fchain.mu_index[imu] for imu in range(fchain.nmu) if fchain.mu_isOverlap]
                # goodMu4 = [ imu4 for imu4 in mu4CandObj_vtxMass if\
                #             filter( lambda x: x in indexBadMuons,\
                #                     [fchain.mu4_quad_l1_index[imu4],fchain.mu4_quad_l2_index[imu4],\
                #                      fchain.mu4_quad_l3_index[imu4],fchain.mu4_quad_l3_index[imu4] ])]
                # if len(goodMu4)>0:
                #         ngoodMu4+=1
                #         print "[output] evt = %d, nMu4Cand_vtx = %d, nMu4Cand_vtxMass = %d, ngoodMu4 = %d"\
                #                 %(fchain.evt, len(mu4Cand_vtx), len(mu4CandObj_vtxMass), len(goodMu4))
                #-->



                


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
