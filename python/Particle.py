import ROOT

class Particle(object):


    def __init__(self, tchain, lep):
        self.LV = ROOT.TLorentzVector()
        self.LV.SetPtEtaPhiM(tchain.lep_pt[lep],tchain.lep_eta[lep],tchain.lep_phi[lep],tchain.lep_mass[lep])
        self.pdg = tchain.lep_pdgId[lep]

    def p4(self):
        return self.LV
    
    def pdgId(self):
        return self.pdg
    
    def mass(self):
        return self.LV.M()
    
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

    def phi(self):
        return self.LV.Phi()

    def rapidity(self):
        return self.LV.Rapidity()

    def __str__(self):
        tmp = '{className} : {pdgId:>3}, pt={pt:5.1f}, eta={eta:5.2f}, phi={phi:5.2f}, mass={mass:5.2f}'
        return tmp.format( className = self.__class__.__name__,
                           pdgId = self.pdgId(),
                           pt = self.pt(),
                           eta = self.eta(),
                           phi = self.phi(),
                           mass = self.mass() )
    