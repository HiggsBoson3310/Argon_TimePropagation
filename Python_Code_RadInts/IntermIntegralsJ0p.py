
import NewHydrogenic as nh
import NewRomberg as nr
import numpy as np
import Radial_Integrals as RI
import os

CMperAU = 2.1947463136320e5
ne, Elo, dEcm = np.loadtxt("EcontEV_det.dat")
EcontAU = np.arange(Elo, Elo+ne*dEcm, dEcm)/CMperAU

J0p_p = os.path.join( os.getcwd(), '..', 'interm_bounds_J0p.dat')
J0p = np.transpose(np.loadtxt(J0p_p))
RI.radint_interm(EcontAU, 1, J0p[0]/CMperAU, 0, 1)