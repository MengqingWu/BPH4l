import math
import ROOT
from deltar import deltaR, deltaPhi

class Pair(object):
    # based on ntuple event:
    def __init__(self,leg1,leg2,pdg = 0):
        self.l1=ROOT.TLorentzVector()
        self.l1.SetPtEtaPhiM(leg1.pt(),leg1.eta(),leg1.phi(),leg1.mass())
        self.l2=ROOT.TLorentzVector()
        self.l2.SetPtEtaPhiM(leg2.p4().pt(),leg2.p4().eta(),leg2.p4().phi(),leg2.p4().mass())

        self.leg1 = leg1
        self.leg2 = leg2
        self.pdg = pdg
        self.LV = self.l1+self.l2
        et1 = math.sqrt(self.l1.M()*self.l1.M()+self.l1.Pt()*self.l1.Pt())
        et2 = math.sqrt(self.l1.M()*self.l1.M()+self.l2.Pt()*self.l2.Pt())

        self.MT  =math.sqrt(self.l1.M()*self.l1.M()+self.l1.M()*self.l1.M()+2*(et1*et2-self.l1.Px()*self.l2.Px()-self.l1.Py()*self.l2.Py()))


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
    
    def mt2(self):
        return self.MT*self.MT

    def mt(self):
        return self.MT

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

