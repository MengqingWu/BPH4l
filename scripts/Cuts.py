#!/usr/bin/env python
## ---------------------------------------------
## Author: Mengqing Wu @2016Dec02
## Cuts functions used for MuonAna.py
## ---------------------------------------------
import math
import ROOT
from python.Particle import Particle

#def passJpsiSubYCut(tree, mu4_i):
class quadMuMassCut(object):
    '''
    mode can be 'both', 'offY', 'offJpsi' 
    '''
    def __init__(self):

        self.jpsi_mass = 3.096916 # pdg mass in GeV
        self.upsilon_mass = 9.46030 # pdg mass in GeV
        self.massErr_sf = 1.105 # see AN2013_099_v17: derived from a fit to the mass pull distribution

    #@staticmethod
    def run(self, tree, mu4_i, mode='both'):
        dimuPass=[]
        if tree.mu4_1a_fitMassErr2[mu4_i]<0 or tree.mu4_1b_fitMassErr2[mu4_i]<0 or tree.mu4_2a_fitMassErr2[mu4_i]<0 or tree.mu4_2b_fitMassErr2[mu4_i]<0:
            ''' should be really RARE! if so, bad evt to discard '''
            return dimuPass

        pair1, pair2 = self.getPairsIn4Mu(tree, mu4_i)
        for ipair in [pair1, pair2]:
            if mode=='offY':
                if self.m_JpsiSubYCut(ipair): dimuPass.append(ipair)
            elif mode=='offJpsi':
                if self.m_YSubJpsiCut(ipair): dimuPass.append(ipair)
            else: continue
        return dimuPass
    
    def m_YSubJpsiCut(self, pair):
        ''' 
        no spectator in Jpsi(off-shell)+Upsilon(on-shell):
        pair = [diMu_Particle, diMu_Particle], ordered in fitted mass (high to low) 
        '''
        if abs(pair[0].fitMass - self.upsilon_mass) < 3*pair[0].fitMassErr and\
           pair[1].fitMass < self.jpsi_mass - 3*pair[1].fitMassErr:
            pair[0].SetTag('onY')
            pair[1].SetTag('offJpsi')
            return True
        else: return False
            
    def m_JpsiSubYCut(self, pair):
        ''' 
        J/psi spectator in Jpsi(on-shell)+Upsilon(off-shell):
        spectator identified as a 3-sigma within nominal Jpsi pdg mass (and such evt is to removed)
        pair = [diMu_Particle, diMu_Particle], ordered in fitted mass (high to low)
        '''        
        njpsi=nY=0
        for dimu in pair:
            if abs(dimu.fitMass - self.jpsi_mass) < 3*dimu.fitMassErr:
                dimu.SetTag('onJpsi'); njpsi+=1
            elif dimu.fitMass < self.upsilon_mass - 3*dimu.fitMassErr:
                dimu.SetTag('offY');
            else: continue
        if njpsi==nY==1:
            return True
        else:
            print '[debug]: onshell njpsi = %d; offshell nY = %d' % (njpsi, nY)
            return False

    def getPairsIn4Mu(self, tree, mu4_i):
        # _1a=( tree.mu4_1a_fitMass[mu4_i],  math.sqrt(tree.mu4_1a_fitMassErr2[mu4_i])*self.massErr_sf )
        # _1b=( tree.mu4_1b_fitMass[mu4_i],  math.sqrt(tree.mu4_1b_fitMassErr2[mu4_i])*self.massErr_sf )
        # _2a=( tree.mu4_2a_fitMass[mu4_i],  math.sqrt(tree.mu4_2a_fitMassErr2[mu4_i])*self.massErr_sf )
        # _2b=( tree.mu4_2b_fitMass[mu4_i],  math.sqrt(tree.mu4_2b_fitMassErr2[mu4_i])*self.massErr_sf )
    
        _1a = Particle(pt=tree.mu4_1a_pt[mu4_i], eta=tree.mu4_1a_eta[mu4_i], phi=tree.mu4_1a_phi[mu4_i], mass=tree.mu4_1a_mass[mu4_i])
        _1b = Particle(pt=tree.mu4_1b_pt[mu4_i], eta=tree.mu4_1b_eta[mu4_i], phi=tree.mu4_1b_phi[mu4_i], mass=tree.mu4_1b_mass[mu4_i])
        _2a = Particle(pt=tree.mu4_2a_pt[mu4_i], eta=tree.mu4_2a_eta[mu4_i], phi=tree.mu4_2a_phi[mu4_i], mass=tree.mu4_2a_mass[mu4_i])
        _2b = Particle(pt=tree.mu4_2b_pt[mu4_i], eta=tree.mu4_2b_eta[mu4_i], phi=tree.mu4_2b_phi[mu4_i], mass=tree.mu4_2b_mass[mu4_i])
    
        setattr(_1a,'fitMass',tree.mu4_1a_fitMass[mu4_i]); setattr(_1a,'fitMassErr',math.sqrt(tree.mu4_1a_fitMassErr2[mu4_i])*self.massErr_sf)
        setattr(_1b,'fitMass',tree.mu4_1b_fitMass[mu4_i]); setattr(_1b,'fitMassErr',math.sqrt(tree.mu4_1b_fitMassErr2[mu4_i])*self.massErr_sf)
        setattr(_2a,'fitMass',tree.mu4_2a_fitMass[mu4_i]); setattr(_2a,'fitMassErr',math.sqrt(tree.mu4_2a_fitMassErr2[mu4_i])*self.massErr_sf)
        setattr(_2b,'fitMass',tree.mu4_2b_fitMass[mu4_i]); setattr(_2b,'fitMassErr',math.sqrt(tree.mu4_2b_fitMassErr2[mu4_i])*self.massErr_sf)

        # Particle type stored in list
        pair1=[_1a, _1b] 
        pair2=[_2a, _2b]
        pair1.sort(key = lambda x: x.fitMass, reverse = True)
        pair2.sort(key = lambda x: x.fitMass, reverse = True)
        #print "[output] pair1 mass(%.4f, %.4f)" % (pair1[0][0], pair1[1][0])
        #print "[output] pair2 mass(%.4f, %.4f)" % (pair2[0][0], pair2[1][0])
    
        return pair1, pair2


    
