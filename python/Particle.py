#!/usr/bin/env python
## ---------------------------------------------
## Author: Mengqing Wu @ 2016 Dec 02
## ---------------------------------------------

import ROOT

class Particle(object):

    def __init__(self, tag='', pt=0, eta=0, phi=0, mass=0, pdg=0):
        ''' 
        -> particle template for offline python analyzer 
        -> open for new attr. attached
        '''

        if pt*eta*phi*mass:
            self.LV = ROOT.TLorentzVector()
            self.LV.SetPtEtaPhiM(pt,eta,phi,mass)
        else: self.LV=None

        self.Tag=tag
        self.Pt=pt
        self.Phi=phi
        self.Eta=eta
        self.Mass=mass
        self.pdg = pdg

    def SetTag(self, newtag):
        if newtag:
            self.Tag = newtag
            return True
        else: return False
    
    def p4(self):    return self.LV
    def pdgId(self): return self.pdg
    def tag(self):   return self.Tag 

    def pt(self):    return self.Pt
    def eta(self):   return self.Eta
    def phi(self):   return self.Phi
    def mass(self):  return self.Mass
    
    def px(self):
        if LV: return self.LV.Px()
        else: return None
        
    def py(self):
        if LV: return self.LV.Py()
        else: return None
        
    def pz(self):
        if LV: return self.LV.Pz()
        else: return None
        
    def rapidity(self):
        if LV: return self.LV.Rapidity()
        else: return None
        
    def __str__(self):
        tmp = '{className} : {pdgId:>3}, pt={pt:5.1f}, eta={eta:5.2f}, phi={phi:5.2f}, mass={mass:5.2f}'
        return tmp.format( className = self.__class__.__name__,
                           pdgId = self.pdgId(),
                           pt = self.pt(),
                           eta = self.eta(),
                           phi = self.phi(),
                           mass = self.mass() )
    
