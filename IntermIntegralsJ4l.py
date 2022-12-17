
import NewHydrogenic as nh
import NewRomberg as nr
import numpy as np
import Radial_Integrals as RI
import os

CMperAU = 2.1947463136320e5
ne, Elo, dEcm = np.loadtxt("EcontEV_det.dat")
EcontAU = np.arange(Elo, Elo+ne*dEcm, dEcm)/CMperAU

J4l_p = os.path.join( os.getcwd(), '..', 'interm_bounds_J4l.dat')
J4l = np.transpose(np.loadtxt(J4l_p))
RI.radint_interm(EcontAU, 3, J4l[0,:]/CMperAU, 4, 3)
RI.radint_interm(EcontAU, 5, J4l[0,:]/CMperAU, 4, 3)