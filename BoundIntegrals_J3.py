import NewHydrogenic as nh
import NewRomberg as nr
import numpy as np
import Radial_Integrals as RI
import os


CMperAU = RI.CMperAU

J3_p = os.path.join( os.getcwd(),"..",'init_bounds_J3.dat')
J3 = np.transpose(np.loadtxt(J3_p))
J2f_p = os.path.join( os.getcwd(), '..', 'interm_bounds_J2f.dat')
J2f = np.transpose(np.loadtxt(J2f_p))
J2p_p = os.path.join( os.getcwd(), '..', 'interm_bounds_J2p.dat')
J2p = np.transpose(np.loadtxt(J2p_p))
J4l_p = os.path.join( os.getcwd(), '..', 'interm_bounds_J4l.dat')
J4l = np.transpose(np.loadtxt(J4l_p))

for eee in J3[0]:
    print(RI.fnu(eee,RI.IP12*CMperAU), RI.fnu(eee,RI.IP32*CMperAU), 3, 5)

RI.wfunc(2, 3, J2f[0]/CMperAU, 3, J3[0]/CMperAU)
RI.wfunc(2, 1, J2p[0]/CMperAU, 3, J3[0]/CMperAU)
RI.wfunc(4, 3, J4l[0]/CMperAU, 3, J3[0]/CMperAU)

