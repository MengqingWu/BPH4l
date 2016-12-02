#!/usr/bin/env python
## ---------------------------------------------
## Author: Mengqing Wu @2016Dec02
## Cuts functions used for MuonAna.py
## ---------------------------------------------
import math
from python.Particle import Particle


jpsi_mass = 3.096916 # pdg mass in GeV
upsilon_mass = 9.46030 # pdg mass in GeV

massErr_sf = 1.105 # see AN2013_099_v17: derived from a fit to the mass pull distribution

def passJpsiSubYCut(tree, mu4_i):
        # pair = [(m, m_sigma),(m, m_sigma)], ordered in mass (high to low)
        m_JpsiSubYCut = (lambda pair: abs(pair[0].fitMass-jpsi_mass)<3*pair[0].fitMassErr and\
                         abs(pair[1].fitMass-upsilon_mass)>3*pair[1].fitMassErr)

        dimuPass=[]

        if tree.mu4_1a_fitMassErr2[mu4_i]<0 or tree.mu4_1b_fitMassErr2[mu4_i]<0 or tree.mu4_2a_fitMassErr2[mu4_i]<0 or tree.mu4_2b_fitMassErr2[mu4_i]<0:
            return dimuPass

        # _1a=( tree.mu4_1a_fitMass[mu4_i],  math.sqrt(tree.mu4_1a_fitMassErr2[mu4_i])*massErr_sf )
        # _1b=( tree.mu4_1b_fitMass[mu4_i],  math.sqrt(tree.mu4_1b_fitMassErr2[mu4_i])*massErr_sf )
        # _2a=( tree.mu4_2a_fitMass[mu4_i],  math.sqrt(tree.mu4_2a_fitMassErr2[mu4_i])*massErr_sf )
        # _2b=( tree.mu4_2b_fitMass[mu4_i],  math.sqrt(tree.mu4_2b_fitMassErr2[mu4_i])*massErr_sf )

        _1a = Particle(pt=tree.mu4_1a_pt[mu4_i], eta=tree.mu4_1a_eta[mu4_i], phi=tree.mu4_1a_phi[mu4_i], mass=tree.mu4_1a_mass[mu4_i])
        _1b = Particle(pt=tree.mu4_1b_pt[mu4_i], eta=tree.mu4_1b_eta[mu4_i], phi=tree.mu4_1b_phi[mu4_i], mass=tree.mu4_1b_mass[mu4_i])
        _2a = Particle(pt=tree.mu4_2a_pt[mu4_i], eta=tree.mu4_2a_eta[mu4_i], phi=tree.mu4_2a_phi[mu4_i], mass=tree.mu4_2a_mass[mu4_i])
        _2b = Particle(pt=tree.mu4_2b_pt[mu4_i], eta=tree.mu4_2b_eta[mu4_i], phi=tree.mu4_2b_phi[mu4_i], mass=tree.mu4_2b_mass[mu4_i])

        setattr(_1a,'fitMass',tree.mu4_1a_fitMass[mu4_i]); setattr(_1a,'fitMassErr',math.sqrt(tree.mu4_1a_fitMassErr2[mu4_i])*massErr_sf)
        setattr(_1b,'fitMass',tree.mu4_1b_fitMass[mu4_i]); setattr(_1b,'fitMassErr',math.sqrt(tree.mu4_1b_fitMassErr2[mu4_i])*massErr_sf)
        setattr(_2a,'fitMass',tree.mu4_2a_fitMass[mu4_i]); setattr(_2a,'fitMassErr',math.sqrt(tree.mu4_2a_fitMassErr2[mu4_i])*massErr_sf)
        setattr(_2b,'fitMass',tree.mu4_2b_fitMass[mu4_i]); setattr(_2b,'fitMassErr',math.sqrt(tree.mu4_2b_fitMassErr2[mu4_i])*massErr_sf)
        
        pair1=[_1a, _1b]
        pair2=[_2a, _2b]
        pair1.sort(key = lambda x: x.fitMass, reverse = True)
        pair2.sort(key = lambda x: x.fitMass, reverse = True)
        #print "[output] pair1 mass(%.4f, %.4f)" % (pair1[0][0], pair1[1][0])
        #print "[output] pair2 mass(%.4f, %.4f)" % (pair2[0][0], pair2[1][0])
        
        if m_JpsiSubYCut(pair1): dimuPass.append(pair1)
        if m_JpsiSubYCut(pair2): dimuPass.append(pair2)
        return dimuPass

    
