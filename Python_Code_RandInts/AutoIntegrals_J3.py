import NewHydrogenic as nh
import NewRomberg as nr
import numpy as np
import Radial_Integrals as RI
import os
#### radint_init_c(Econt, Jc, lc, Eo, Jo, noffset=0)

CMperAU = RI.CMperAU
ne, Elo, dEcm = np.loadtxt("EcontEV_AI.dat")
EcontAU = np.arange(Elo, Elo+ne*dEcm, dEcm)/CMperAU

J3_p = os.path.join( os.getcwd(),"..",'init_bounds_J3.dat')
J3 = np.transpose(np.loadtxt(J3_p))

J2f_p = os.path.join( os.getcwd(),"..","test_of_zs_J2f.dat")
J2f = np.transpose(np.loadtxt(J2f_p))
RI.radint_init_c(J2f[0]/CMperAU, 2, 3, J3[0]/CMperAU, 3, noffset=0)

J2p_p = os.path.join( os.getcwd(),"..","test_of_zs_J2p.dat")
J2p = np.transpose(np.loadtxt(J2p_p))
RI.radint_init_c(J2p[0]/CMperAU, 2, 1, J3[0]/CMperAU, 3, noffset=0)
