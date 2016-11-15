import math
import ROOT
from deltar import deltaR, deltaPhi

class Pair(object):
    # based on ntuple event:
    #  tchain is the tree,
    #  leg1 and leg2 refers to the object index.
    def __init__(self, tchain, leg1, leg2):
        self.l1=ROOT.TLorentzVector()
        self.l1.SetPtEtaPhiM(tchain.lep_pt[leg1],tchain.lep_eta[leg1],tchain.lep_phi[leg1],tchain.lep_mass[leg1])
        self.l2=ROOT.TLorentzVector()
        self.l2.SetPtEtaPhiM(tchain.lep_pt[leg2],tchain.lep_eta[leg2],tchain.lep_phi[leg2],tchain.lep_mass[leg2])

        # attention: only return the index
        self.leg1 = leg1 
        self.leg2 = leg2
        if abs(tchain.lep_pdgId[leg1]) == abs(tchain.lep_pdgId[leg2]): self.pdg = tchain.lep_pdgId[leg1]
        else: print "[ERROR!] the flavour of two lep in the Pair is different: lep1->%d and lep2->%d " % (tchain.lep_pdgId[leg1], tchain.lep_pdgId[leg2])

        self.LV = self.l1+self.l2
        #et1 = math.sqrt(self.l1.M()*self.l1.M()+self.l1.Pt()*self.l1.Pt())
        #et2 = math.sqrt(self.l1.M()*self.l1.M()+self.l2.Pt()*self.l2.Pt())

    def rawP4(self):
        return self.l1+self.l2

    def p4(self):
        return self.LV

    # invariant mass
    def m(self):
        return self.LV.M()

    def mass(self):
        return self.LV.M()

    def pdgId(self):
        return self.pdg
    
    def pt(self):
        return self.LV.Pt()

    def px(self):
        return self.LV.Px()

    def py(self):
        return self.LV.Py()

    def pz(self):
        return self.LV.Pz()

    def eta(self):
        return self.LV.Eta()

    def rapidity(self):
        return self.LV.Rapidity()

    def phi(self):
        return ROOT.TVector2.Phi_mpi_pi(self.LV.Phi())

    def AbsdeltaPhi(self):
        return abs(deltaPhi(self.l1.Phi(),self.l2.Phi()))

    def AbsdeltaR(self):
        return abs(deltaR(self.l1.Eta(),self.l1.Phi(),self.l2.Eta(),self.l2.Phi()))

    def deltaPhi(self):
        return deltaPhi(self.l1.Phi(),self.l2.Phi())

    def deltaR(self):
        return deltaR(self.l1.Eta(),self.l1.Phi(),self.l2.Eta(),self.l2.Phi())

    def __getattr__(self, name):
        return getattr(self.LV,name)

