import numpy as np
import sympy.physics.wigner as wg
import matplotlib.pyplot as plt
import sys
import fft_processing as fp

# As reference the order of the channels is the following for the argon case:
# in Jcs we write as [(jc)Jcs,l]J:
# J=1: [(1/2)1,2]1 [(3/2)2,2]1 [(3/2)1,2]1 [(1/2)1,0]1 [(3/2)1,0]1
# J=3: [(3/2)1,2]3 [(3/2)2,2]3 [(1/2)1,2]3 [(3/2)1,4]J [(3/2)2,4]3 [(1/2)1,4]3
# J=5: [(3/2)1,4]5 [(3/2)2,4]5 [(1/2)1,4]5 [(3/2)1,6]5 [(3/2)2,6]5 [(1/2)1,6]5

l = [[2,2,2,0,0], 
     [2, 2, 2, 4, 4, 4], 
     [4,4,4,6,6,6]]

jc = [[1/2., 3/2., 3/2.,1/2., 3/2.],
      [3/2., 3/2., 1/2., 3/2., 3/2., 1/2.],
      [3/2., 3/2., 1/2., 3/2., 3/2., 1/2.]]

jcs = [[1, 2, 1, 1, 1], 
       [1, 2, 1, 1, 2, 1],
       [1, 2, 1, 1, 2, 1]]

I_12 = [[0,3],[2,5],[2,5]]

# number of continuum states 
#ncon = 200
J = [1,3,5]
nchan = [5,6,6]

def ak12(k,cs_con,ncon):
    ak = np.zeros(ncon,dtype=complex)
    # Start with the sum over the Jc=1/2 channels
    # There is only one value of Jcs=1
    Jcs = 1
    Jindices = [[0,3],[2,5],[2,5]]
    for j in range(len(J)):
        for jp in range(len(J)):
            for i in Jindices[j]:
                for ip in Jindices[jp]:
                    #print("For term with jcs %i, %i %i, wiht partial waves %i %i using range %i %i"
                    # %(Jcs,J[j],J[jp],l[j][i],l[jp][ip],ncon*(sum(nchan[:j])+i),ncon*(sum(nchan[:jp])+ip)))
                    term = np.sqrt((2*l[j][i]+1)*(2*l[jp][ip]+1)*(2*J[j]+1)*(2*J[jp]+1))
                    term = term * (-1)**(l[j][i]-l[jp][ip]+Jcs+k)
                    term = term * float(wg.wigner_3j(l[j][i],l[jp][ip],k,0,0,0) * wg.wigner_3j(J[jp],J[j],k,0,0,0))
                    term = term * float(wg.wigner_6j(l[jp][ip],l[j][i],k,J[j],J[jp],Jcs))
                    term = term *cs_con[ncon*(sum(nchan[:j])+i):ncon*(sum(nchan[:j])+i+1)]*np.conjugate(cs_con[ncon*(sum(nchan[:jp])+ip):ncon*(sum(nchan[:jp])+ip+1)])
                    ak += term
    ak = ak*(2*k+1)/(4*np.pi)
    return ak

def p12(lpw,m,cs_con,ncon):
    p12 = np.zeros(ncon,dtype=complex)
    Jc=1/2
    Jcs = 1
    Jindices = [[0,3],[2,5],[2,5]]
    for j in range(len(J)):
        for jp in range(len(J)):
            for i in Jindices[j]:
                for ip in Jindices[jp]:
                    if(lpw==l[j][i] and lpw==l[jp][ip]):
                        term = float(wg.clebsch_gordan(Jcs,lpw,J[j],-m,m,0))*float(wg.clebsch_gordan(Jcs,lpw,J[jp],-m,m,0))
                        term = term * cs_con[ncon*(sum(nchan[:j])+i):ncon*(sum(nchan[:j])+i+1)]
                        term = term * np.conjugate(cs_con[ncon*(sum(nchan[:jp])+ip):ncon*(sum(nchan[:jp])+ip+1)])
                        p12 += term
    return p12

def pJ12(Jpw,cs_con,ncon):
    p12 = np.zeros(ncon,dtype=complex)
    Jc=1/2
    Jcs = 1
    Jindices = [[0,3],[2,5],[2,5]]
    for jcs in range(len(Jcs)):
        j=J.index(Jpw)
        jp = J.index(Jpw)
        for i in Jindices[j]:
            for ip in Jindices[jp]:
                term = cs_con[ncon*(sum(nchan[:j])+i):ncon*(sum(nchan[:j])+i+1)]
                term = term * np.conjugate(cs_con[ncon*(sum(nchan[:jp])+ip):ncon*(sum(nchan[:jp])+ip+1)])
                p12 += term
    return p12

I_32 = [[1,2,4],[0,1,3,4],[0,1,3,4]]

def ak32(k,cs_con,ncon):
    ak = np.zeros(ncon,dtype=complex)
    # Start with the sum over the Jc=1/2 channels
    # There is only one value of Jcs=1
    Jcs = [1,2]
    Jindices = [[[2,4],[0,3],[0,3]],[[1],[1,4],[1,4]]]
    for jcs in range(len(Jcs)):
        for j in range(len(J)):
            for jp in range(len(J)):
                for i in Jindices[jcs][j]:
                    for ip in Jindices[jcs][jp]:
                        #print("For term with jcs %i, %i %i, wiht partial waves %i %i using range %i %i"%
                        # (Jcs[jcs],J[j],J[jp],l[j][i],l[jp][ip],ncon*(sum(nchan[:j])+i),ncon*(sum(nchan[:jp])+ip)))
                        term = np.sqrt((2*l[j][i]+1)*(2*l[jp][ip]+1)*(2*J[j]+1)*(2*J[jp]+1))
                        term = term * (-1)**(l[j][i]-l[jp][ip]+Jcs[jcs]+k)
                        term = term * float(wg.wigner_3j(l[j][i],l[jp][ip],k,0,0,0) * wg.wigner_3j(J[jp],J[j],k,0,0,0))
                        term = term * float(wg.wigner_6j(l[jp][ip],l[j][i],k,J[j],J[jp],Jcs[jcs]))
                        term = term *cs_con[ncon*(sum(nchan[:j])+i):ncon*(sum(nchan[:j])+i+1)]*np.conjugate(cs_con[ncon*(sum(nchan[:jp])+ip):ncon*(sum(nchan[:jp])+ip+1)])
                        ak += term
    ak = ak*(2*k+1)/(4*np.pi)
    return ak


def p32(lpw,m,cs_con,ncon):
    p32 = np.zeros(ncon,dtype=complex)
    Jc=3/2
    Jcs = [1,2]
    Jindices = [[[2,4],[0,3],[0,3]],[[1],[1,4],[1,4]]]
    for jcs in range(len(Jcs)):
        for j in range(len(J)):
            for jp in range(len(J)):
                for i in Jindices[jcs][j]:
                    for ip in Jindices[jcs][jp]:
                        if(lpw==l[j][i] and lpw==l[jp][ip]):
                            term = float(wg.clebsch_gordan(Jcs[jcs],lpw,J[j],-m,m,0))*float(wg.clebsch_gordan(Jcs[jcs],lpw,J[jp],-m,m,0))
                            term = term * cs_con[ncon*(sum(nchan[:j])+i):ncon*(sum(nchan[:j])+i+1)]
                            term = term * np.conjugate(cs_con[ncon*(sum(nchan[:jp])+ip):ncon*(sum(nchan[:jp])+ip+1)])
                            p32 += term
    return p32

def pJ32(Jpw,cs_con,ncon):
    p32 = np.zeros(ncon,dtype=complex)
    Jc=3/2
    Jcs = [1,2]
    Jindices = [[[2,4],[0,3],[0,3]],[[1],[1,4],[1,4]]]
    for jcs in range(len(Jcs)):
        j=J.index(Jpw)
        jp = J.index(Jpw)
        for i in Jindices[jcs][j]:
            for ip in Jindices[jcs][jp]:
                term = cs_con[ncon*(sum(nchan[:j])+i):ncon*(sum(nchan[:j])+i+1)]
                term = term * np.conjugate(cs_con[ncon*(sum(nchan[:jp])+ip):ncon*(sum(nchan[:jp])+ip+1)])
                p32 += term
    return p32

def generate_a_specs(cdl,tos,esampled,D,ncon):
    nlim = 8+24+30+17+26+5*36
    ccon = cdl[:,nlim:]
    #print(np.shape(ccon))
    a012 = np.zeros((ncon,len(tos)),dtype=complex)
    a032 = np.zeros((ncon,len(tos)),dtype=complex)


    a212 = np.zeros((ncon,len(tos)),dtype=complex)
    a232 = np.zeros((ncon,len(tos)),dtype=complex)

    a412 = np.zeros((ncon,len(tos)),dtype=complex)
    a432 = np.zeros((ncon,len(tos)),dtype=complex)

    a612 = np.zeros((ncon,len(tos)),dtype=complex)
    a632 = np.zeros((ncon,len(tos)),dtype=complex)

    a812 = np.zeros((ncon,len(tos)),dtype=complex)
    a832 = np.zeros((ncon,len(tos)),dtype=complex)

    a1012 = np.zeros((ncon,len(tos)),dtype=complex)
    a1032 = np.zeros((ncon,len(tos)),dtype=complex)

    # compute for each
    for p in range(len(tos)):
        a012[:,p] = ak12(0,ccon[p],ncon)
        a032[:,p] = ak32(0,ccon[p],ncon)

        a212[:,p] = ak12(2,ccon[p],ncon)
        a232[:,p] = ak32(2,ccon[p],ncon)

        a412[:,p] = ak12(4,ccon[p],ncon)
        a432[:,p] = ak32(4,ccon[p],ncon)

        a612[:,p] = ak12(6,ccon[p],ncon)
        a632[:,p] = ak32(6,ccon[p],ncon)

        a812[:,p] = ak12(8,ccon[p],ncon)
        a832[:,p] = ak32(8,ccon[p],ncon)

        a1012[:,p] = ak12(10,ccon[p],ncon)
        a1032[:,p] = ak32(10,ccon[p],ncon)

    esam = esampled[nlim:nlim+ncon]
    egrid = np.arange(min(esampled[nlim]-fp.IP12,esampled[nlim]-fp.IP32),
                    max(esampled[nlim+ncon-1]-fp.IP12,esampled[nlim+ncon-1]-fp.IP32),
                    esampled[nlim+1]-esampled[nlim])

    eexp, tdels, a0 = fp.one_spectrum_w(esam, tos, np.real(a012), np.real(a032),D)
    eexp, tdels, a2 = fp.one_spectrum_w(esam, tos, np.real(a212), np.real(a232),D)
    eexp, tdels, a4 = fp.one_spectrum_w(esam, tos, np.real(a412), np.real(a432),D)
    eexp, tdels, a6 = fp.one_spectrum_w(esam, tos, np.real(a612), np.real(a632),D)
    eexp, tdels, a8 = fp.one_spectrum_w(esam, tos, np.real(a812), np.real(a832),D)
    eexp, tdels, a10 = fp.one_spectrum_w(esam, tos, np.real(a1012), np.real(a1032),D)
    return eexp, egrid, tdels, a0,a2,a4,a6,a8,a10

def generate_a_specs_sep(cdl,tos,esampled,D,ncon):
    nlim = 8+24+30+17+26+5*36
    ccon = cdl[:,nlim:]
    #print(np.shape(ccon))
    a012 = np.zeros((ncon,len(tos)),dtype=complex)
    a032 = np.zeros((ncon,len(tos)),dtype=complex)


    a212 = np.zeros((ncon,len(tos)),dtype=complex)
    a232 = np.zeros((ncon,len(tos)),dtype=complex)

    a412 = np.zeros((ncon,len(tos)),dtype=complex)
    a432 = np.zeros((ncon,len(tos)),dtype=complex)

    a612 = np.zeros((ncon,len(tos)),dtype=complex)
    a632 = np.zeros((ncon,len(tos)),dtype=complex)

    a812 = np.zeros((ncon,len(tos)),dtype=complex)
    a832 = np.zeros((ncon,len(tos)),dtype=complex)

    a1012 = np.zeros((ncon,len(tos)),dtype=complex)
    a1032 = np.zeros((ncon,len(tos)),dtype=complex)

    # compute for each
    for p in range(len(tos)):
        a012[:,p] = ak12(0,ccon[p],ncon)
        a032[:,p] = ak32(0,ccon[p],ncon)

        a212[:,p] = ak12(2,ccon[p],ncon)
        a232[:,p] = ak32(2,ccon[p],ncon)

        a412[:,p] = ak12(4,ccon[p],ncon)
        a432[:,p] = ak32(4,ccon[p],ncon)

        a612[:,p] = ak12(6,ccon[p],ncon)
        a632[:,p] = ak32(6,ccon[p],ncon)

        a812[:,p] = ak12(8,ccon[p],ncon)
        a832[:,p] = ak32(8,ccon[p],ncon)

        a1012[:,p] = ak12(10,ccon[p],ncon)
        a1032[:,p] = ak32(10,ccon[p],ncon)

    esam = esampled[nlim:nlim+ncon]
    egrid = np.arange(min(esampled[nlim]-fp.IP12,esampled[nlim]-fp.IP32),
                    max(esampled[nlim+199]-fp.IP12,esampled[nlim+199]-fp.IP32),
                    esampled[nlim+1]-esampled[nlim])

    eexp, tdels, a0_12, a0_32 = fp.one_spectrum_w_sep(esam, tos, np.real(a012), np.real(a032),D)
    eexp, tdels, a2_12, a2_32 = fp.one_spectrum_w_sep(esam, tos, np.real(a212), np.real(a232),D)
    eexp, tdels, a4_12, a4_32 = fp.one_spectrum_w_sep(esam, tos, np.real(a412), np.real(a432),D)
    eexp, tdels, a6_12, a6_32 = fp.one_spectrum_w_sep(esam, tos, np.real(a612), np.real(a632),D)
    eexp, tdels, a8_12, a8_32 = fp.one_spectrum_w_sep(esam, tos, np.real(a812), np.real(a832),D)
    eexp, tdels, a10_12, a10_32 = fp.one_spectrum_w_sep(esam, tos, np.real(a1012), np.real(a1032),D)
    return eexp, egrid, tdels, a0_12,a0_32,a2_12,a2_32,a4_12,a4_32,a6_12,a6_32,a8_12,a8_32,a10_12,a10_32

def betas(itime,e,a0,a2,a4,a6,a8,a10,egrid):
    b2 = a2[:,itime]/a0[:,itime]
    b4 = a4[:,itime]/a0[:,itime]
    b6 = a6[:,itime]/a0[:,itime]
    b8 = a8[:,itime]/a0[:,itime]
    b10 = a10[:,itime]/a0[:,itime]
    if(e>egrid[0] and e<egrid[-1]):
        nms = np.argsort(np.abs(egrid-e))[:2]
        ms = np.array([(b2[nms[0]]-b2[nms[1]])/(egrid[nms[0]]-egrid[nms[1]]),(b4[nms[0]]-b4[nms[1]])/(egrid[nms[0]]-egrid[nms[1]]),
              (b6[nms[0]]-b6[nms[1]])/(egrid[nms[0]]-egrid[nms[1]]),(b8[nms[0]]-b8[nms[1]])/(egrid[nms[0]]-egrid[nms[1]]),
              (b10[nms[0]]-b10[nms[1]])/(egrid[nms[0]]-egrid[nms[1]])])
        bs = np.array([b2[nms[0]]-ms[0]*egrid[nms[0]],b4[nms[0]]-ms[1]*egrid[nms[0]],b6[nms[0]]-ms[2]*egrid[nms[0]],b8[nms[0]]-ms[3]*egrid[nms[0]],
              b10[nms[0]]-ms[4]*egrid[nms[0]]])
        return ms*e+bs        
    else:
        return np.zeros(5)

def generate_plm(lw,cdl,esampled,D,ncon):
    nlim = 8+24+30+17+26+5*36
    ccon = cdl[nlim:]
    p0_32 = p32(lw,0,ccon,ncon)
    p0_12 = p12(lw,0,ccon,ncon)

    p1_32 = p32(lw,1,ccon,ncon)
    p1_12 = p12(lw,1,ccon,ncon)

    p2_32 = p32(lw,2,ccon,ncon)
    p2_12 = p12(lw,2,ccon,ncon)

    esam = esampled[nlim:nlim+ncon]
    #print("Splitting: ",(fp.IP12-fp.IP32)*fp.CMperAU)
    egrid = np.arange(min(esampled[nlim]-fp.IP12,esampled[nlim]-fp.IP32),
                    max(esampled[nlim+ncon-1]-fp.IP12,esampled[nlim+ncon-1]-fp.IP32),
                    esampled[nlim+1]-esampled[nlim])

    eexp, p012,p032 = fp.single_conv_w_sep(esam, np.real(p0_12), np.real(p0_32),D)
    eexp, p112,p132 = fp.single_conv_w_sep(esam, np.real(p1_12), np.real(p1_32),D)
    eexp, p212,p232 = fp.single_conv_w_sep(esam, np.real(p2_12), np.real(p2_32),D)

    return eexp, egrid, p012,p032, p112,p132, p212, p232    

def generate_plm(lw,cdl,esampled,D,ncon):
    nlim = 8+24+30+17+26+5*36
    ccon = cdl[nlim:]
    p0_32 = p32(lw,0,ccon,ncon)
    p0_12 = p12(lw,0,ccon,ncon)

    p1_32 = p32(lw,1,ccon,ncon)
    p1_12 = p12(lw,1,ccon,ncon)

    p2_32 = p32(lw,2,ccon,ncon)
    p2_12 = p12(lw,2,ccon,ncon)

    esam = esampled[nlim:nlim+ncon]
    #print("Splitting: ",(fp.IP12-fp.IP32)*fp.CMperAU)
    egrid = np.arange(min(esampled[nlim]-fp.IP12,esampled[nlim]-fp.IP32),
                    max(esampled[nlim+ncon-1]-fp.IP12,esampled[nlim+ncon-1]-fp.IP32),
                    esampled[nlim+1]-esampled[nlim])

    eexp, p012,p032 = fp.single_conv_w_sep(esam, np.real(p0_12), np.real(p0_32),D)
    eexp, p112,p132 = fp.single_conv_w_sep(esam, np.real(p1_12), np.real(p1_32),D)
    eexp, p212,p232 = fp.single_conv_w_sep(esam, np.real(p2_12), np.real(p2_32),D)

    return eexp, egrid, p012,p032, p112,p132, p212, p232

def acs(itime,e,a0,a2,a4,a6,a8,a10,egrid):
    aa0 = a0[:,itime]
    b2 = a2[:,itime]
    b4 = a4[:,itime]
    b6 = a6[:,itime]
    b8 = a8[:,itime]
    b10 = a10[:,itime]
    if(e>egrid[0] and e<egrid[-1]):
        nms = np.argsort(np.abs(egrid-e))[:2]
        ms = np.array([(aa0[nms[0]]-aa0[nms[1]])/(egrid[nms[0]]-egrid[nms[1]]),(b2[nms[0]]-b2[nms[1]])/(egrid[nms[0]]-egrid[nms[1]]),
                       (b4[nms[0]]-b4[nms[1]])/(egrid[nms[0]]-egrid[nms[1]]),(b6[nms[0]]-b6[nms[1]])/(egrid[nms[0]]-egrid[nms[1]]),
                       (b8[nms[0]]-b8[nms[1]])/(egrid[nms[0]]-egrid[nms[1]]),(b10[nms[0]]-b10[nms[1]])/(egrid[nms[0]]-egrid[nms[1]])])
        bs = np.array([aa0[nms[0]]-ms[0]*egrid[nms[0]],b2[nms[0]]-ms[1]*egrid[nms[0]],
                       b4[nms[0]]-ms[2]*egrid[nms[0]],b6[nms[0]]-ms[3]*egrid[nms[0]],
                       b8[nms[0]]-ms[4]*egrid[nms[0]],b10[nms[0]]-ms[5]*egrid[nms[0]]])
        return ms*e+bs        
    else:
        return np.zeros(6)
