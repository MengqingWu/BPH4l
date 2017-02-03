import math
import ROOT
#from Pair import Pair
from Particle import Particle

def GetPairsIn4Mu( tree, mu4_i, massErr_sf = 1.105):
    
        # _1a=( tree.mu4_1a_fitMass[mu4_i],  math.sqrt(tree.mu4_1a_fitMassErr2[mu4_i])*massErr_sf )
        # _1b=( tree.mu4_1b_fitMass[mu4_i],  math.sqrt(tree.mu4_1b_fitMassErr2[mu4_i])*massErr_sf )
        # _2a=( tree.mu4_2a_fitMass[mu4_i],  math.sqrt(tree.mu4_2a_fitMassErr2[mu4_i])*massErr_sf )
        # _2b=( tree.mu4_2b_fitMass[mu4_i],  math.sqrt(tree.mu4_2b_fitMassErr2[mu4_i])*massErr_sf )
    
        _1a = Particle(pt=tree.mu4_1a_pt[mu4_i], eta=tree.mu4_1a_eta[mu4_i], phi=tree.mu4_1a_phi[mu4_i], mass=tree.mu4_1a_mass[mu4_i])
        _1b = Particle(pt=tree.mu4_1b_pt[mu4_i], eta=tree.mu4_1b_eta[mu4_i], phi=tree.mu4_1b_phi[mu4_i], mass=tree.mu4_1b_mass[mu4_i])
        _2a = Particle(pt=tree.mu4_2a_pt[mu4_i], eta=tree.mu4_2a_eta[mu4_i], phi=tree.mu4_2a_phi[mu4_i], mass=tree.mu4_2a_mass[mu4_i])
        _2b = Particle(pt=tree.mu4_2b_pt[mu4_i], eta=tree.mu4_2b_eta[mu4_i], phi=tree.mu4_2b_phi[mu4_i], mass=tree.mu4_2b_mass[mu4_i])
        
        if tree.mu4_1a_fitMassErr2[mu4_i]<0 or tree.mu4_1b_fitMassErr2[mu4_i]<0 or tree.mu4_2a_fitMassErr2[mu4_i]<0 or tree.mu4_2b_fitMassErr2[mu4_i]<0:
                ''' should be really RARE! if you want to use fitMass, this is a bad evt to discard '''
                pass
        else:
                setattr(_1a,'fitMass',tree.mu4_1a_fitMass[mu4_i]); setattr(_1a,'fitMassErr',math.sqrt(tree.mu4_1a_fitMassErr2[mu4_i])*massErr_sf)
                setattr(_1b,'fitMass',tree.mu4_1b_fitMass[mu4_i]); setattr(_1b,'fitMassErr',math.sqrt(tree.mu4_1b_fitMassErr2[mu4_i])*massErr_sf)
                setattr(_2a,'fitMass',tree.mu4_2a_fitMass[mu4_i]); setattr(_2a,'fitMassErr',math.sqrt(tree.mu4_2a_fitMassErr2[mu4_i])*massErr_sf)
                setattr(_2b,'fitMass',tree.mu4_2b_fitMass[mu4_i]); setattr(_2b,'fitMassErr',math.sqrt(tree.mu4_2b_fitMassErr2[mu4_i])*massErr_sf)

        setattr(_1a,'index',[tree.mu4_1a_l1_index[mu4_i],tree.mu4_1a_l2_index[mu4_i]]); 
        setattr(_1b,'index',[tree.mu4_1b_l1_index[mu4_i],tree.mu4_1b_l2_index[mu4_i]]);         
        setattr(_2a,'index',[tree.mu4_2a_l1_index[mu4_i],tree.mu4_2a_l2_index[mu4_i]]); 
        setattr(_2b,'index',[tree.mu4_2b_l1_index[mu4_i],tree.mu4_2b_l2_index[mu4_i]]);         
        
        out={'1a':_1a, '1b':_1b, '2a':_2a, '2b':_2b}
        return out

    
class Quad(object):
    def __init__(self, tchain, imu4):
            
        self.LV = ROOT.TLorentzVector()
        self.LV.SetPtEtaPhiM(tchain.mu4_quad_pt[imu4],tchain.mu4_quad_eta[imu4],tchain.mu4_quad_phi[imu4],tchain.mu4_quad_mass[imu4])

        self.vtxProb = tchain.mu4_quad_vtxProb[imu4]
        self.vtxChi2 = tchain.mu4_quad_vtxChi2[imu4]
        self.tIndex = imu4
    
    def p4(self):
        return self.LV
        
    def mass(self):
        return self.LV.M()
    
    def pt(self):
        return self.LV.Pt()

    # def px(self):
    #     return self.LV.Px()

    # def py(self):
    #     return self.LV.Py()

    # def pz(self):
    #     return self.LV.Pz()

    def eta(self):
        return self.LV.Eta()

    def phi(self):
        return self.LV.Phi()

    # def rapidity(self):
    #     return self.LV.Rapidity()
    
    # def __str__(self):
    #     tmp = '{className} : {pdgId:>3}, pt={pt:5.1f}, eta={eta:5.2f}, phi={phi:5.2f}, mass={mass:5.2f}'
    #     return tmp.format( className = self.__class__.__name__,
    #                        pdgId = self.pdgId(),
    #                        pt = self.pt(),
    #                        eta = self.eta(),
    #                        phi = self.phi(),
    #                        mass = self.mass() )
    
