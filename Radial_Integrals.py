from matplotlib.pyplot import cm
import NewHydrogenic as nh
import NewRomberg as nr
import numpy as np
import os

#Constants
RydAr = 109735.81
eVperAU = 27.211
CMperAU = 2.1947463136320e5
IP32 = 127109.88/CMperAU
IP12= IP32+(1431.41)/CMperAU
print(IP32,IP12)
#Energy of the boundstates
E5s = 14.0864375/eVperAU  
E3d = 14.1458475/eVperAU
E5sp= 14.2565375/eVperAU
E3dp= 14.3124075/eVperAU
Eb = [E5s, E3d, E5sp, E3dp]
EconteV = np.linspace(0.05,1.5,300)+IP12*eVperAU
EcontAU =  EconteV/eVperAU
EcontCM = EcontAU*CMperAU
sp = ["s","p","d","f"]

# Write the continuum energies sampled wavlengths
with open("EcontEV_det.dat", mode='w') as fileEv:
    fileEv.write("%i    %.8e    %.8e \n"%(len(EconteV),EcontCM[0], EcontCM[1]-EcontCM[0]))

# Effective quantum numbers
def fnu(Ecm, Icm):
    return np.sqrt(RydAr/(Icm-Ecm))

# String formatter
def formatme(float_tuple):
    str = ""
    for f in float_tuple:
        str = str+"%.8e "%f
    return str 

# Bound to continuum radial integrals for the five channel J=1 states.
def array_cont(l, J):
    dirname = "RI b(J=1) c(J%i%s)"%(J,sp[l])
    os.makedirs(dirname, exist_ok=True)
    bnd_chan = lambda x: [[x-IP12,2],[x-IP32,2],[x-IP32,2],[x-IP12,0],[x-IP32,0]]
    if(J==2 and l==3):
        cont_chan = lambda x: [[x-IP32,3],[x-IP32,3],[x-IP12,3]]
    elif(J==2 and l==1):
        cont_chan = lambda x: [[x-IP32,1],[x-IP32,1],[x-IP12,1]]
    elif(J==0 and l==1):
        cont_chan = lambda x: [[x-IP32,1],[x-IP12,1]]
    else:
        print("Wrong J value for this process")
        return None

    for ii,ee in enumerate(Eb):
        i_chan = bnd_chan(ee)
        file_r = open(dirname+"/cont_rad_ints_regb%i.dat"%(ii), mode='w')
        file_i = open(dirname+"/cont_rad_ints_iregb%i.dat"%(ii), mode='w')
        for ec in EcontAU:
            c_chan = cont_chan(ec)
            liner = ""
            linei = ""
            for ecc, lcc in c_chan:
                liner += formatme([nr.lcont_RadIntRombRegWhit(ecc,lcc,ebb,lbb,False) for ebb, lbb in i_chan])+"    "
                linei += formatme([nr.lcont_RadIntRombIregWhit(ecc,lcc,ebb,lbb,False) for ebb, lbb in i_chan])+"    "
            liner += "\n"
            linei += "\n"
            file_r.write(liner)
            file_i.write(linei)
        file_r.close()
        file_i.close()
        print("Done with initial bound state %i"%(ii))
    return 0

# Regular function integrals for the one photon code.
def regular_int(l):
    dir_name = "Rads %s"%sp[l]
    os.makedirs(dir_name, exist_ok=True)
    file1 = open(dir_name+"/Radial_e%s_bs_5s.dat"%(sp[l]), mode="w")
    file2 = open(dir_name+"/Radial_e%s_bd_5s.dat"%(sp[l]), mode="w")

    file3 = open(dir_name+"/Radial_e%s_bs_3d.dat"%(sp[l]), mode="w")
    file4 = open(dir_name+"/Radial_e%s_bd_3d.dat"%(sp[l]), mode='w')

    file5 = open(dir_name+"/Radial_e%s_bs_5sp.dat"%(sp[l]), mode="w")
    file6 = open(dir_name+"/Radial_e%s_bd_5sp.dat"%(sp[l]), mode="w")

    file7 = open(dir_name+"/Radial_e%s_bs_3dp.dat"%(sp[l]), mode="w")
    file8 = open(dir_name+"/Radial_e%s_bd_3dp.dat"%(sp[l]), mode='w')

    file = [[file1,file2],[file3,file4],[file5,file6], [file7, file8]]

    for e, files in zip(Eb, file):
        print("Computing for %.8f"%(e))
        for ee in EcontAU:
            e12 = ee-IP12
            e32 = ee-IP32
            e12_b = e-IP12
            e32_b = e-IP32
            r1 = nr.lcont_RadIntRombRegWhit(e12,l,e12_b,0, False)
            r2 = nr.lcont_RadIntRombRegWhit(e32,l,e12_b,0, False)
            r3 = nr.lcont_RadIntRombRegWhit(e12,l,e32_b,0, False)
            r4 = nr.lcont_RadIntRombRegWhit(e32,l,e32_b,0, False)
            files[0].write("%.8f    %.8e    %.8e    %.8e    %.8e \n"%(ee*CMperAU, r1, r2, r3, r4))

            r1 = nr.lcont_RadIntRombRegWhit(e12,l,e12_b,2, False)
            r2 = nr.lcont_RadIntRombRegWhit(e32,l,e12_b,2, False)
            r3 = nr.lcont_RadIntRombRegWhit(e12,l,e32_b,2, False)
            r4 = nr.lcont_RadIntRombRegWhit(e32,l,e32_b,2, False)
            files[1].write("%.8f    %.8e    %.8e    %.8e    %.8e \n"%(ee*CMperAU, r1, r2, r3, r4))
        
        files[1].close()
        files[0].close()

# Irregular function integrals for the one photon code.
def iregular_int(l):
    dir_name = "Rads %s"%sp[l]
    os.makedirs(dir_name, exist_ok=True)
    file1 = open(dir_name+"/Radial_Ir_e%s_bs_5s.dat"%(sp[l]), mode="w")
    file2 = open(dir_name+"/Radial_Ir_e%s_bd_5s.dat"%(sp[l]), mode="w")

    file3 = open(dir_name+"/Radial_Ir_e%s_bs_3d.dat"%(sp[l]), mode="w")
    file4 = open(dir_name+"/Radial_Ir_e%s_bd_3d.dat"%(sp[l]), mode='w')

    file5 = open(dir_name+"/Radial_Ir_e%s_bs_5sp.dat"%(sp[l]), mode="w")
    file6 = open(dir_name+"/Radial_Ir_e%s_bd_5sp.dat"%(sp[l]), mode="w")

    file7 = open(dir_name+"/Radial_Ir_e%s_bs_3dp.dat"%(sp[l]), mode="w")
    file8 = open(dir_name+"/Radial_Ir_e%s_bd_3dp.dat"%(sp[l]), mode='w')

    file = [[file1,file2],[file3,file4],[file5,file6], [file7, file8]]

    for e, files in zip(Eb, file):
        print("Computing for %.8f"%(e))
        for ee in EcontAU:
            e12 = ee-IP12
            e32 = ee-IP32
            e12_b = e-IP12
            e32_b = e-IP32
            
            r1 = nr.lcont_RadIntRombIregWhit(e12,l,e12_b,0, False)
            r2 = nr.lcont_RadIntRombIregWhit(e32,l,e12_b,0, False)
            r3 = nr.lcont_RadIntRombIregWhit(e12,l,e32_b,0, False)
            r4 = nr.lcont_RadIntRombIregWhit(e32,l,e32_b,0, False)
            files[0].write("%.8f    %.8e    %.8e    %.8e    %.8e \n"%(ee*CMperAU, r1, r2, r3, r4))

            r1 = nr.lcont_RadIntRombIregWhit(e12,l,e12_b,2, False)
            r2 = nr.lcont_RadIntRombIregWhit(e32,l,e12_b,2, False)
            r3 = nr.lcont_RadIntRombIregWhit(e12,l,e32_b,2, False)
            r4 = nr.lcont_RadIntRombIregWhit(e32,l,e32_b,2, False)
            files[1].write("%.8f    %.8e    %.8e    %.8e    %.8e \n"%(ee*CMperAU, r1, r2, r3, r4))
        
        files[1].close()
        files[0].close()

# Integrals between bound states, inial o to intermediate int
def wfunc(Jint, linti, Eint, Jo, Eo):
    dir_name = "RI b(J=%i)  b(J=%i lint %s)"%(Jo, Jint, sp[linti])
    os.makedirs(dir_name, exist_ok=True)
    if(abs(Jint-Jo)!=1): 
        raise ValueError("Wrong combination of initial and intermediate J values")

    else:  
        for label_i, ee in enumerate(Eo):
            ee12 = ee-IP12
            ee32 = ee-IP32
            file = open(dir_name+"/interm_rad_ints_%i.dat"%(label_i), mode='w')
            for e in Eint:
                e12 = e-IP12
                e32 = e-IP32

                if(Jint==0):
                    int_chan = [[e32,1], [e12,1]]
                elif(Jint==2):
                    int_chan = [[e32, linti],[e32,linti],[e12,linti]]
                elif(Jint==4):
                    int_chan = [[e32,3],[e32,3],[e12,3],[e32,5],[e32,5],[e12,5]]
                else:
                    raise ValueError("Incorrect Intermediate J")
                
                if(Jo==1):
                    ini_chan = [[ee12,2], [ee32,2], [ee32,2], [ee12,0], [ee32,0]]
                elif(Jo==3):
                    ini_chan = [[ee32,2], [ee32,2], [ee12,2], [ee32,4], [ee32,4], [ee12,4]]
                else:
                    raise ValueError("Incorrect Initial J")

                line = ""
                for i, [ei,lint] in enumerate(int_chan):
                    results = np.zeros(len(ini_chan))
                    for kk, [eo, lo] in enumerate(ini_chan):
                        results[kk] = nr.lcont_RadIntRombWhitWhit(eo,lo, ei,lint,False) if abs(lo-lint)==1 else 0
                    line = line + formatme(tuple(results))
                    if(i+1==len(int_chan)): line = line+'\n'
                file.write(line)
            file.close()

# Integrals from intermediate states to continuum
def radint_interm(Econt, Jc, Eint, Jint, lint, noffset=0):
    dir_name = "RI b(J=%i,l=%s) c(J=%i)"%(Jint,sp[lint], Jc)
    os.makedirs(dir_name, exist_ok=True)
    if(abs(Jc-Jint)!=1):
        raise ValueError("Incorrect combination of continuum and intermediate J values")
    if(Jc==1):
        cont_chan = lambda x: [[x-IP12,2],[x-IP32,2],[x-IP32,2],[x-IP12,0],[x-IP32,0]]
    elif(Jc==3):
        cont_chan = lambda z: [[z-IP32,2],[z-IP32,2],[z-IP12,2],[z-IP32,4],[z-IP32,4],[z-IP12,4]]
    elif(Jc==5):
        cont_chan = lambda z: [[z-IP32,4],[z-IP32,4],[z-IP12,4],[z-IP32,6],[z-IP32,6],[z-IP12,6]]
    else:
        raise ValueError("Incorrect value of J for the continuum")
        return None
    
    if(Jint==2):
        int_chan = lambda z: [[z-IP32,lint],[z-IP32,lint],[z-IP12,lint]]
    elif(Jint==0):
        int_chan = lambda z: [[z-IP32, 1],[z-IP12, 1]]
    elif(Jint==4):
        int_chan = lambda z: [[z-IP32, 3],[z-IP32, 3],[z-IP12, 3],[z-IP32,5],[z-IP32,5],[z-IP12,5]]
    else:
        raise ValueError("Incorrect value of J for intermediate state")
    
    for eei in Eint:
        i_chan = int_chan(eei)
        iii = np.where(Eint==eei)
        iii = iii[0]+noffset
        print("Currently at %i"%(iii))
        filer = open(dir_name+"/cont_rad_ints_reg_b%i.dat"%(iii), mode='w')
        filei = open(dir_name+"/cont_rad_ints_ireg_b%i.dat"%(iii), mode='w')
        for eec in Econt:
            c_chan = cont_chan(eec)
            liner = ""
            linei = ""
            for i, [ec,lc] in enumerate(c_chan):
                resultsr = [nr.lcont_RadIntRombRegWhit(ec,lc,ei,li,False) for [ei, li] in i_chan]
                resultsi = [nr.lcont_RadIntRombIregWhit(ec,lc,ei,li,False) for [ei, li] in i_chan]
                liner = liner + formatme(tuple(resultsr))
                linei = linei + formatme(tuple(resultsi))
                if(i+1==len(c_chan)): 
                    linei = linei+'\n'
                    liner = liner+'\n'
            filer.write(linei)
            filei.write(liner)
        filer.close()
        filei.close()
        print("done with %i"%(iii))
    return 0

# Integrals from initial bound to continuum
def radint_init_c(Econt, Jc, lc, Eo, Jo, noffset=0):
    dir_name = "RI b(J=%i) c(J=%i,,l=%s)"%(Jo,Jc,sp[lc])
    os.makedirs(dir_name, exist_ok=True)
    if(abs(Jc-Jo)!=1):
        raise ValueError("Incorrect combination of continuum and intermediate J values")
    if(Jc==0):
        cont_chan = lambda z: [[z-IP32, 1],[z-IP12, 1]] 
    elif(Jc==2):
        cont_chan = lambda z: [[z-IP32,lc],[z-IP32,lc],[z-IP12,lc]]
    elif(Jc==4):
        cont_chan = lambda z: [[z-IP32, 3],[z-IP32, 3],[z-IP12, 3],[z-IP32,5],[z-IP32,5],[z-IP12,5]]
    else:
        raise ValueError("Incorrect value of J for the continuum")
        return None
    
    if(Jo==1):
        o_chan = lambda x: [[x-IP12,2],[x-IP32,2],[x-IP32,2],[x-IP12,0],[x-IP32,0]]
    elif(Jo==3):
        o_chan = lambda z:[[z-IP32,2],[z-IP32,2],[z-IP12,2],[z-IP32,4],[z-IP32,4],[z-IP12,4]]
    else:
        raise ValueError("Incorrect value of deep initial J")
    
    for eei in Eo:
        o_chan = o_chan(eei)
        iii = np.where(Eo==eei)
        iii = iii[0]+noffset
        print("Currently at %i"%(iii))
        filer = open(dir_name+"/cont_rad_ints_reg_b%i.dat"%(iii), mode='w')
        filei = open(dir_name+"/cont_rad_ints_ireg_b%i.dat"%(iii), mode='w')
        for eec in Econt:
            c_chan = cont_chan(eec)
            liner = ""
            linei = ""
            for i, [ec,lc] in enumerate(c_chan):
                resultsr = [nr.lcont_RadIntRombRegWhit(ec,lc,ei,li,False) for [ei, li] in o_chan]
                resultsi = [nr.lcont_RadIntRombIregWhit(ec,lc,ei,li,False) for [ei, li] in o_chan]
                liner = liner + formatme(tuple(resultsr))
                linei = linei + formatme(tuple(resultsi))
                if(i+1==len(c_chan)): 
                    linei = linei+'\n'
                    liner = liner+'\n'
            filer.write(linei)
            filei.write(liner)
        filer.close()
        filei.close()
        print("done with %i"%(iii))
    return 0

# Integrals with information printed in between numbers to identify appropriate number of lines. 
def radint_interm_labelled(Econt, Jc, Eint, Jint, lint, labels):
    dir_name = "RI b(J=%i,l=%s) c(J=%i)"%(Jint,sp[lint], Jc)
    os.makedirs(dir_name, exist_ok=True)
    
    if(Jc==1):
        cont_chan = lambda x: [[x-IP12,2],[x-IP32, 2],[x-IP32,2], [x-IP12,0],[x-IP32,0]]
    elif(Jc==3):
        cont_chan = lambda z: [[z-IP32,2],[z-IP32,2],[z-IP12,2],[z-IP32,4],[z-IP32,4],[z-IP12,4]]
    else:
        print("Correct the value of J for the continuum")
        return None
    
    if(Jint==2):
        int_chan = lambda z: [[z-IP32,lint],[z-IP32,lint],[z-IP12,lint]]
    elif(Jint==0 and lint==1):
        int_chan = lambda z: [[z-IP32, lint], [z-IP12, lint]]
    else:
        print("Correct the valu of J intermediate")
        return None
    
    for eei in Eint:
        i_chan = int_chan(eei)
        iii = np.where(Eint==eei)
        iii = labels[iii]
        print("Currently at "+str(iii))
        filer = open(dir_name+"/cont_rad_ints_reg_b%s.dat"%(iii), mode='w')
        filei = open(dir_name+"/cont_rad_ints_ireg_b%s.dat"%(iii), mode='w')
        for eec in Econt:
            c_chan = cont_chan(eec)
            liner = ""
            linei = ""
            for i, [ec,lc] in enumerate(c_chan):
                print("%i Energia del continuo %i es %.8f"%(i,lc,ec))
                resultsr = [nr.lcont_RadIntRombRegWhit(ec,lc,ei,li,False) for [ei, li] in i_chan]
                resultsi = [nr.lcont_RadIntRombIregWhit(ec,lc,ei,li,False) for [ei, li] in i_chan]
                liner = liner + formatme(tuple(resultsr))
                linei = linei + formatme(tuple(resultsi))
                if(i+1==len(c_chan)): 
                    print("fin de la linea")
                    linei = linei+'\n'
                    liner = liner+'\n'
            filer.write(linei)
            filei.write(liner)
        filer.close()
        filei.close()
        print("done with %s"%(iii))
    return 0

# Def plot by closed channel pair

def plot_pair(E1, l1, E2, l2):
    ch1 = [[E1-IP12, l1],[E1-IP32, l1]]
    ch2 = [[E2-IP12, l2],[E2-IP32, l2]]
    print("Channel 1")
    print(ch1)
    print("Channel 2")
    print(ch2)
    for i1,c1 in enumerate(ch1):
        for i2,c2 in enumerate(ch2):
            name = "Trouble Shoot/chan%i%i_l%i%i.dat"%(i1,i2,l1,l2)
            nr.lcont_RadIntRombWhitWhit(c1[0],c1[1],c2[0],c2[1], True, name=name)
    
    return None





# Data to plot continuum to bound behaviour of the integrals.
def plotem():
    e12 = EcontAU[-1]-IP12
    e32 = EcontAU[-1]-IP32
    e12_b = E5s-IP12
    e32_b = E5s-IP32
    print("s")
    r1 = nr.lcont_RadIntRombIregWhit(e12,3,e12_b,0, True)
    r2 = nr.lcont_RadIntRombIregWhit(e12,1,e12_b,0, True)
    print(r1,r2)
    r3 = nr.lcont_RadIntRombRegWhit(e12,3,e12_b,0, True)
    r4 = nr.lcont_RadIntRombRegWhit(e12,1,e12_b,0, True)
    print(r3,r4)
    print("d")
    r1 = nr.lcont_RadIntRombIregWhit(e12,3,e12_b,2, True)
    r2 = nr.lcont_RadIntRombIregWhit(e12,1,e12_b,2, True)
    print(r1,r2)
    r3 = nr.lcont_RadIntRombRegWhit(e12,3,e12_b,2, True)
    r4 = nr.lcont_RadIntRombRegWhit(e12,1,e12_b,2, True)
    print(r3,r4)

    return 0


