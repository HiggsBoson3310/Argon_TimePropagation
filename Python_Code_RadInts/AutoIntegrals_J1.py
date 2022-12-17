import NewHydrogenic as nh
import NewRomberg as nr
import numpy as np
import Radial_Integrals as RI
import os
#### radint_init_c(Econt, Jc, lc, Eo, Jo, noffset=0)

CMperAU = RI.CMperAU
ne, Elo, dEcm = np.loadtxt("EcontEV_AI.dat")
EcontAU = np.arange(Elo, Elo+ne*dEcm, dEcm)/CMperAU

J1_p = os.path.join( os.getcwd(),"..",'init_bounds_J1.dat')
J1 = np.transpose(np.loadtxt(J1_p))

J2p_p = os.path.join( os.getcwd(),"..","test_of_zs_J2p.dat")
J2p = np.transpose(np.loadtxt(J2p_p))
RI.radint_init_c(J2p[0]/CMperAU, 2, 1, J1[0]/CMperAU, 1, noffset=0)

J2f_p = os.path.join( os.getcwd(),"..","test_of_zs_J2f.dat")
J2f = np.transpose(np.loadtxt(J2f_p))
RI.radint_init_c(J2f[0]/CMperAU, 2, 3, J1[0]/CMperAU, 1, noffset=0)

J0_p = os.path.join( os.getcwd(),"..","test_of_zs_J0.dat")
J0 = np.transpose(np.loadtxt(J0_p))
RI.radint_init_c(J0[0]/CMperAU, 0, 1, J1[0]/CMperAU, 1, noffset=0)
