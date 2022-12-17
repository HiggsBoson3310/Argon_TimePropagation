import numpy as np

def formatme(float_tuple):
    str = ""
    for f in float_tuple:
        str = str+"%.8e "%f
    return str 


no=6
nc=3
ne = 200
for i in range(24):
    print(i)
    name = "RI b(J=2,l=p) c(J=3)/cont_rad_ints_ireg_b%i.dat"%(i)
    fileun = open(name)
    elems = []
    for line in fileun.readlines():
        ll = line.replace(" \n", "").split(" ")
        if ("" in ll):
            ll.remove("")
        #print(ll)
        [elems.append(k) for k in np.array(ll,dtype=float)]
        #print(elems)
    fileun.close()
    elems = np.array(elems).reshape(200,(nc*no))
    file = open(name, mode="w")
    for ee in elems:
        file.write(formatme(tuple(ee))+"\n")
    file.close()

    name = "RI b(J=2,l=p) c(J=3)/cont_rad_ints_reg_b%i.dat"%(i)
    fileun = open(name)
    elems = []
    for line in fileun.readlines():
        ll = line.replace(" \n", "").split(" ")
        if ("" in ll):
            ll.remove("")
        #print(ll)
        [elems.append(k) for k in np.array(ll,dtype=float)]
        #print(elems)
    fileun.close()
    elems = np.array(elems).reshape(200,(nc*no))
    file = open(name, mode="w")
    for ee in elems:
        file.write(formatme(tuple(ee))+"\n")

    file.close()

