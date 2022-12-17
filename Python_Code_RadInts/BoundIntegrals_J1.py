import NewHydrogenic as nh
import NewRomberg as nr
import numpy as np
import Radial_Integrals as RI
import os


CMperAU = RI.CMperAU

J1_p = os.path.join( os.getcwd(),"..",'init_bounds_J1.dat')
J1 = np.transpose(np.loadtxt(J1_p))
J2f_p = os.path.join( os.getcwd(), '..', 'interm_bounds_J2f.dat')
J2f = np.transpose(np.loadtxt(J2f_p))
J2p_p = os.path.join( os.getcwd(), '..', 'interm_bounds_J2p.dat')
J2p = np.transpose(np.loadtxt(J2p_p))
J0p_p = os.path.join( os.getcwd(), '..', 'interm_bounds_J0p.dat')
J0p = np.transpose(np.loadtxt(J0p_p))

RI.wfunc(2, 3, J2f[0]/CMperAU, 1, J1[0]/CMperAU)
RI.wfunc(2, 1, J2p[0]/CMperAU, 1, J1[0]/CMperAU)
RI.wfunc(0, 1, J0p[0]/CMperAU, 1, J1[0]/CMperAU)

