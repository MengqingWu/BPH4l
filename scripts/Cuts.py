#!/usr/bin/env python
## ---------------------------------------------
## Author: Mengqing Wu @2016Dec02
## Cuts functions used for MuonAna.py
## ---------------------------------------------
import math
import ROOT
from python.Particle import Particle
from python.Quad import GetPairsIn4Mu

min2muVtx=0.005
min4muVtx=0.05

quadMuVtxCut = (lambda tree,index: tree.mu4_quad_vtxProb[index]>min4muVtx )
# at least 1 pair of diMus both w/diMuVtxProb>0.005
diMuVtxCut = (lambda tree,index: (tree.mu4_1a_vtxProb[index]>min4muVtx and tree.mu4_1b_vtxProb[index]>min4muVtx) or\
              (tree.mu4_2a_vtxProb[index]>min4muVtx and tree.mu4_2b_vtxProb[index]>min4muVtx))


#def passJpsiSubYCut(tree, mu4_i):
class quadMuMassCut(object):
    '''
    mode can be 'all', 'offY', 'offJpsi', 'bothOn' 
    '''
    def __init__(self):

        self.jpsi_mass = 3.096916 # pdg mass in GeV
        self.upsilon_mass = 9.46030 # pdg mass in GeV
        self.massErr_sf = 1.105 # see AN2013_099_v17: derived from a fit to the mass pull distribution

    #@staticmethod
    def run(self, tree, mu4_i, mode='all'):
        dimuPass=[]
        if tree.mu4_1a_fitMassErr2[mu4_i]<0 or tree.mu4_1b_fitMassErr2[mu4_i]<0 or tree.mu4_2a_fitMassErr2[mu4_i]<0 or tree.mu4_2b_fitMassErr2[mu4_i]<0:
            ''' should be really RARE! if so, bad evt to discard '''
            return dimuPass

        pair1, pair2 = self.getPairsIn4Mu(tree, mu4_i)
        for ipair in [pair1, pair2]:
            if mode == 'offY' or mode == 'all': 
                if self.m_JpsiSubYCut(ipair): dimuPass.append(ipair)
                
            if mode == 'offJpsi' or mode == 'all':
                if self.m_YSubJpsiCut(ipair): dimuPass.append(ipair)

            if mode == 'bothOn' or mode == 'all':
                if self.m_JpsiYCut(ipair): dimuPass.append(ipair)
                
        return dimuPass

    def m_JpsiYCut(self, pair):
        '''  mode == 'bothOn' for Jpsi(on-shell) + Upsilon(on-shell):
        pair = [diMu_Particle, diMu_Particle], ordered in fitted mass (high to low) 
        '''
        #print "[info] using 'off - Jpsi' "
        if abs(pair[0].fitMass - self.upsilon_mass) < 3*pair[0].fitMassErr and\
           abs(pair[1].fitMass - self.jpsi_mass) <  3*pair[1].fitMassErr:
            pair[0].SetTag('onY')
            pair[1].SetTag('onJpsi')
            return True
        else: return False
        
    def m_YSubJpsiCut(self, pair):
        '''  mode == 'offJpsi'
        no spectator in Jpsi(off-shell)+Upsilon(on-shell):
        pair = [diMu_Particle, diMu_Particle], ordered in fitted mass (high to low) 
        '''
        #print "[info] using 'off - Jpsi' "
        if abs(pair[0].fitMass - self.upsilon_mass) < 3*pair[0].fitMassErr and\
           pair[1].fitMass < self.jpsi_mass - 3*pair[1].fitMassErr:
            pair[0].SetTag('onY')
            pair[1].SetTag('offJpsi')
            return True
        else: return False
            
    def m_JpsiSubYCut(self, pair):
        '''  mode == 'offY'
        J/psi spectator in Jpsi(on-shell)+Upsilon(off-shell):
        spectator identified as a 3-sigma within nominal Jpsi pdg mass (and such evt is to removed)
        pair = [diMu_Particle, diMu_Particle], ordered in fitted mass (high to low)
        '''
        #print "[info] using 'off - Y' "
        njpsi=nY=0
        for dimu in pair:
            if abs(dimu.fitMass - self.jpsi_mass) < 3*dimu.fitMassErr:
                dimu.SetTag('onJpsi'); njpsi+=1
            elif dimu.fitMass < self.upsilon_mass - 3*dimu.fitMassErr:
                dimu.SetTag('offY'); nY+=1
            else: continue
        if njpsi==nY==1:
            return True
        else:
            #if njpsi==1: print '[debug]: onshell njpsi = %d; offshell nY = %d; m1, m2 = %.2f, %.2f' % (njpsi,nY,pair[0].fitMass, pair[1].fitMass )
            return False

    def getPairsIn4Mu(self, tree, mu4_i):
        # # _1a=( tree.mu4_1a_fitMass[mu4_i],  math.sqrt(tree.mu4_1a_fitMassErr2[mu4_i])*self.massErr_sf )
        # # _1b=( tree.mu4_1b_fitMass[mu4_i],  math.sqrt(tree.mu4_1b_fitMassErr2[mu4_i])*self.massErr_sf )
        # # _2a=( tree.mu4_2a_fitMass[mu4_i],  math.sqrt(tree.mu4_2a_fitMassErr2[mu4_i])*self.massErr_sf )
        # # _2b=( tree.mu4_2b_fitMass[mu4_i],  math.sqrt(tree.mu4_2b_fitMassErr2[mu4_i])*self.massErr_sf )
    
        # _1a = Particle(pt=tree.mu4_1a_pt[mu4_i], eta=tree.mu4_1a_eta[mu4_i], phi=tree.mu4_1a_phi[mu4_i], mass=tree.mu4_1a_mass[mu4_i])
        # _1b = Particle(pt=tree.mu4_1b_pt[mu4_i], eta=tree.mu4_1b_eta[mu4_i], phi=tree.mu4_1b_phi[mu4_i], mass=tree.mu4_1b_mass[mu4_i])
        # _2a = Particle(pt=tree.mu4_2a_pt[mu4_i], eta=tree.mu4_2a_eta[mu4_i], phi=tree.mu4_2a_phi[mu4_i], mass=tree.mu4_2a_mass[mu4_i])
        # _2b = Particle(pt=tree.mu4_2b_pt[mu4_i], eta=tree.mu4_2b_eta[mu4_i], phi=tree.mu4_2b_phi[mu4_i], mass=tree.mu4_2b_mass[mu4_i])
    
        # setattr(_1a,'fitMass',tree.mu4_1a_fitMass[mu4_i]); setattr(_1a,'fitMassErr',math.sqrt(tree.mu4_1a_fitMassErr2[mu4_i])*self.massErr_sf)
        # setattr(_1b,'fitMass',tree.mu4_1b_fitMass[mu4_i]); setattr(_1b,'fitMassErr',math.sqrt(tree.mu4_1b_fitMassErr2[mu4_i])*self.massErr_sf)
        # setattr(_2a,'fitMass',tree.mu4_2a_fitMass[mu4_i]); setattr(_2a,'fitMassErr',math.sqrt(tree.mu4_2a_fitMassErr2[mu4_i])*self.massErr_sf)
        # setattr(_2b,'fitMass',tree.mu4_2b_fitMass[mu4_i]); setattr(_2b,'fitMassErr',math.sqrt(tree.mu4_2b_fitMassErr2[mu4_i])*self.massErr_sf)

        diMus=GetPairsIn4Mu(tree, mu4_i, massErr_sf=self.massErr_sf )
        
        # Particle type stored in list
        pair1=[diMus['1a'], diMus['1b']] 
        pair2=[diMus['2a'], diMus['2b']]

        pair1.sort(key = lambda x: x.fitMass, reverse = True)
        pair2.sort(key = lambda x: x.fitMass, reverse = True)
        #print "[output] pair1 mass(%.4f, %.4f)" % (pair1[0][0], pair1[1][0])
        #print "[output] pair2 mass(%.4f, %.4f)" % (pair2[0][0], pair2[1][0])
    
        return pair1, pair2


    
