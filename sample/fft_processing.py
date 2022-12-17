from multiprocessing.sharedctypes import Value
import numpy as np
import scipy.fft as fft
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib as mpl
#from matplotlib.collections import LineCollection
#from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
import matplotlib.pylab as pl
import matplotlib.scale as mscale
import matplotlib.transforms as mtransforms
import matplotlib.ticker as ticker
from numba import jit
from numba.typed import List


def trapz_d(y, x):
    return 0.5*((x[1:]-x[:-1])*(y[1:]+y[:-1])).sum()

trapz = jit(trapz_d)


class SquareRootScale(mscale.ScaleBase):
    """
    ScaleBase class for generating square root scale.
    """
 
    name = 'squareroot'
 
    def __init__(self, axis, **kwargs):
        # note in older versions of matplotlib (<3.1), this worked fine.
        # mscale.ScaleBase.__init__(self)

        # In newer versions (>=3.1), you also need to pass in `axis` as an arg
        mscale.ScaleBase.__init__(self, axis)
 
    def set_default_locators_and_formatters(self, axis):
        axis.set_major_locator(ticker.AutoLocator())
        axis.set_major_formatter(ticker.ScalarFormatter())
        axis.set_minor_locator(ticker.NullLocator())
        axis.set_minor_formatter(ticker.NullFormatter())
 
    def limit_range_for_scale(self, vmin, vmax, minpos):
        return  max(0., vmin), vmax
 
    class SquareRootTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
 
        def transform_non_affine(self, a): 
            return np.array(a)**0.5
 
        def inverted(self):
            return SquareRootScale.InvertedSquareRootTransform()
 
    class InvertedSquareRootTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
 
        def transform(self, a):
            return np.array(a)**2
 
        def inverted(self):
            return SquareRootScale.SquareRootTransform()
 
    def get_transform(self):
        return self.SquareRootTransform()
 
mscale.register_scale(SquareRootScale)



auoftime = 2.4188843265857e-17
CMperAU =  2.194746313632e5
auofI = 3.50944758e16
IP32 = 127109.9/CMperAU
IP12= IP32+1432./CMperAU
Ebs = np.array([113614.9731,114096.3442,114986.918,115429.9305])/CMperAU
Ebs_exp = np.array([113643.260, 114147.732, 114975.019, 115366.866 ])/CMperAU

#tdels = np.linspace(0.0, 7.5e4, 1800) 
#tlist, test = np.transpose(np.loadtxt('testing_the_beats.dat'))

#spectrogram_1ph = np.transpose(np.loadtxt("timde_delays_photoionization.dat"))

#Order of the bound states: 5s, 3d, 5s' and 3d'
#tst = fft.fftshift(fft.fftfreq(15000,tdels[1]-tdels[0])*27.211*2*np.pi)


#Energy resolution function. For the autoionization region we used D=0.005
# e is the reported energy and eo is the real energy. One must integrate over eo
# everything assumes atomic units.

def format_me(x):
    return ("%.8f   "*len(x)+"\n")%tuple(x)

def Eresol(e,eo,D):
    if(e<0):
        print("lel")
    return (1/(np.sqrt(np.pi)*D*np.sqrt(e)))*np.exp(-(((e-eo)/(D*np.sqrt(e)))**2))

#Convolution for the photoionization yield
# photoprob should be an array with a row for each sampled energy and a colum for each time delay sampled
# it will be transformed into an arry with a row for each energy in the grid and a column for each time delay.
def convolphoto(photoprob, Egrid, EgridSample, tgrid, D,first=""):
    photoprobconvolved = np.zeros((len(Egrid),len(tgrid)))
    for j in range(len(Egrid)):
        for i in range(len(tgrid)):
            photoprobconvolved[j,i] = trapz(Eresol(Egrid[j],EgridSample,D)*photoprob[:,i],EgridSample)
    np.save("%s.npy"%first,np.column_stack((EgridSample,photoprob[:,0])))
    np.save("%s_cc.npy"%first,np.column_stack((Egrid,photoprobconvolved[:,0])))
    return photoprobconvolved


# Signal convolution for a single time delay or just a driver 
def single_convol(photoprob, Egrid, EgridSample, D):
    photoprobconvolved = np.zeros((len(Egrid)))

    for j in range(len(Egrid)):
        photoprobconvolved[j] = trapz(Eresol(Egrid[j],EgridSample,D)*photoprob,EgridSample)
    
    return photoprobconvolved
def single_conv_w(Esampled, probsJ12, probsJ32,D):
    Es12 = Esampled - IP12
    Es32 = Esampled - IP32
    Eexp = np.arange(min(Es12[0],Es32[0]),max(Es12[-1],Es32[-1]),Esampled[1]-Esampled[0])
    conv_12 = single_convol(probsJ12, Eexp, Es12, D)
    conv_32 = single_convol(probsJ32, Eexp, Es32, D)
    return [Eexp, conv_12+conv_32]

def single_conv_w_sep(Esampled, probsJ12, probsJ32,D):
    Es12 = Esampled - IP12
    Es32 = Esampled - IP32
    Eexp = np.arange(min(Es12[0],Es32[0]),max(Es12[-1],Es32[-1]),Esampled[1]-Esampled[0])
    conv_12 = single_convol(probsJ12, Eexp, Es12, D)
    conv_32 = single_convol(probsJ32, Eexp, Es32, D)
    return [Eexp, conv_12,conv_32]

def fft_eV(tlist, Ard):
    auoftime = 2.4188843265857e-17
    T = (tlist[1]-tlist[0])*1e-12/auoftime
    N = len(tlist)
    freq40 = fft.fftshift(fft.fftfreq(N,T)*27.211*2*np.pi)
    dat_fft = fft.fftshift(fft.fft(Ard)/N)
    return (freq40[np.where(freq40>=0)], dat_fft[np.where(freq40>=0)])

#Atomic units assumed for sampled energies.
def one_spectrum(Esampled, tdels, probsJ12, probsJ32):
    D = 0.005
    Es12 = Esampled - IP12
    Es32 = Esampled - IP32
    Eexp = np.arange(min(Es12[0],Es32[0]),max(Es12[-1],Es32[-1]),Esampled[1]-Esampled[0])
    #print("Crunching numbers and doing some convolutions...")
    conv_12 = convolphoto(probsJ12, Eexp, Es12, tdels, D,first="uno")
    conv_32 = convolphoto(probsJ32, Eexp, Es32, tdels, D,first="duo")
    return [Eexp, tdels, conv_12+conv_32]

def one_spectrum_w(Esampled, tdels, probsJ12, probsJ32,D):
    Es12 = Esampled - IP12
    Es32 = Esampled - IP32
    Eexp = np.arange(min(Es12[0],Es32[0]),max(Es12[-1],Es32[-1]),Esampled[1]-Esampled[0])
    #print("Crunching numbers and doing some convolutions...")
    conv_12 = convolphoto(probsJ12, Eexp, Es12, tdels, D,first="uno")
    conv_32 = convolphoto(probsJ32, Eexp, Es32, tdels, D,first="duo")
    return [Eexp, tdels, conv_12+conv_32]

def one_spectrum_w_sep(Esampled, tdels, probsJ12, probsJ32,D):
    Es12 = Esampled - IP12
    Es32 = Esampled - IP32
    Eexp = np.arange(min(Es12[0],Es32[0]),max(Es12[-1],Es32[-1]),Esampled[1]-Esampled[0])
    #print("Crunching numbers and doing some convolutions...")
    conv_12 = convolphoto(probsJ12, Eexp, Es12, tdels, D,first="uno")
    conv_32 = convolphoto(probsJ32, Eexp, Es32, tdels, D,first="duo")
    return [Eexp, tdels, conv_12,conv_32]

def fft_spectrum_shifted(tps, spectrum):
    freq, first = fft_eV(tps, spectrum[0]-sum(spectrum[0])/len(spectrum[0]))
    fft_spectrum = np.zeros((len(spectrum), len(freq)), dtype=complex)
    fft_spectrum[0] = first
    #print(first[0])
    for line in range(1,len(spectrum)):
        first, fft_spectrum[line] = fft_eV(tps, spectrum[line]-sum(spectrum[line])/len(spectrum[line]))
    return freq, fft_spectrum

def fft_spectrum(tps, spectrum):
    freq, first = fft_eV(tps, spectrum[0])
    fft_spectrum = np.zeros((len(spectrum), len(freq)), dtype=complex)
    fft_spectrum[0] = first
    for line in range(1,len(spectrum)):
        first, fft_spectrum[line] = fft_eV(tps, spectrum[line])
    return freq, fft_spectrum

def fft_spectrum_meaned(tps, spectrum):
    meaned = np.zeros(len(spectrum),dtype=complex)
    meaned[0] = sum(spectrum[0])/len(spectrum[0])
    freq, first = fft_eV(tps, spectrum[0]-meaned[0])
    fft_spectrum = np.zeros((len(spectrum), len(freq)), dtype=complex)
    fft_spectrum[0] = first
    print(first[0])
    for line in range(1,len(spectrum)):
        meaned[line] = sum(spectrum[line])/len(spectrum[line])
        first, fft_spectrum[line] = fft_eV(tps, spectrum[line]-meaned[line])
        
    return freq, fft_spectrum, meaned

# Plotting routines

def one_photon_analysis(spectrogram, tdels, tlist,test):
    Ebs = np.array([0.1136149731E+6,0.1140963442E+6,0.1149869180E+6,0.1154299305E+6])/CMperAU
    freq40, arg = fft_eV(tdels*auoftime/1e-12, spectrogram[0])
    spectro_fft = np.zeros((len(spectrogram), len(freq40)),dtype=complex)
    for i in range(len(spectrogram)):
        freq40, spectro_fft[i] = fft_eV(tdels*auoftime/1e-12, spectrogram[i])

    plt.plot(tlist, test)
    plt.savefig('testing.png')
    plt.close()

    flist, tstfft = fft_eV(tlist*auoftime/1e-12, test)

    plt.plot(flist, abs(tstfft)**2, 'go')
    plt.yscale('log')
    plt.axvline(x=27.211/CMperAU*abs(Ebs[0]-Ebs[1]), c='C0')
    plt.axvline(x=27.211/CMperAU*abs(Ebs[0]-Ebs[2]), c='C1')
    plt.axvline(x=27.211/CMperAU*abs(Ebs[0]-Ebs[3]), c='C2')
    plt.axvline(x=27.211/CMperAU*abs(Ebs[1]-Ebs[2]), c='C3')
    plt.axvline(x=27.211/CMperAU*abs(Ebs[1]-Ebs[3]), c='C4')
    plt.axvline(x=27.211/CMperAU*abs(Ebs[3]-Ebs[2]), c='C5')
    plt.xlim(0,0.5)
    plt.savefig('testing_fft.png')
    plt.close()

    plt.plot(range(len(spectrogram[:,0])),spectrogram[:,0])
    plt.xlim(25,50)
    plt.savefig("vertical_line.png")
    plt.close()

    plt.imshow(spectrogram[::-1,::], aspect=10)
    plt.savefig("sepectrogram.png")
    plt.close()

    beat = abs(Ebs[0]-Ebs[1])*27.211/CMperAU
    df = freq40[1]-freq40[0]
    nf = 30
    plt.imshow(np.log10(np.abs(spectro_fft)[::-1,:nf:1]),aspect=2*nf/len(spectro_fft[:,0]))
    plt.xlim(-2,nf)
    plt.savefig('spectrogram_fft.png')
    plt.close()

    plt.plot(range(len(spectrogram[:,0])),spectrogram[:,0])
    #plt.xlim(25,50)
    plt.savefig("vertical_line.png")
    plt.close()


    plt.plot(freq40,abs(spectro_fft[80,:])**2)
    plt.axvline(beat)
    plt.yscale('log')
    plt.xlim(0,beat*1.05)
    plt.savefig("horizontal_line.png")
    plt.close()

    plt.yscale('linear')
    plt.xlim(0,tdels[-1])
    plt.plot(tdels,spectrogram[80,:])
    plt.axvline(beat)
    plt.savefig("horizontal_line_realtime.png")
    plt.close()
    return 0

def make_the_plots_first_time(name12, name32, outname):
    ncon, Elo, dE = np.loadtxt("Python_Code_RadInts/EcontEV_det.dat")
    Esam= np.arange(Elo, Elo+ncon*dE, dE)/CMperAU
    Esam12 = Esam - IP12
    Esam32 = Esam - IP32
    yaxis= np.arange(min(Esam12[0],Esam32[0]),max(Esam12[-1],Esam32[-1]),Esam[1]-Esam[0])

    #Load the one photon data
    data32_op = np.loadtxt(name32+".dat")
    data12_op = np.loadtxt(name12+".dat")
    xaxis, probs32 = data32_op[:,0], np.transpose(data32_op[:,1:])
    tdels, probs12 = data12_op[:,0], np.transpose(data12_op[:,1:])
    yaxis, xaxis, spec= one_spectrum(Esam, tdels, probs12, probs32)
    np.savetxt(outname+".dat", spec)
    np.savetxt(outname+"xaxis.dat", xaxis)
    np.savetxt(outname+"yaxis.dat", yaxis)

    res = (xaxis[-1]*auoftime/1e-12-xaxis[0]*auoftime/1e-12)/(yaxis[-1]*27.211-yaxis[0]*27.211)

    freqs, spec_fft = fft_spectrum_shifted(xaxis*auoftime/1e-12, spec)
    nlim = int(0.1/(freqs[1]-freqs[0]))

    resfft = (freqs[nlim] - freqs[0])/(yaxis[-1]*27.211-yaxis[0]*27.211)

    
    im = plt.imshow(spec[::-1,:], aspect=res, extent=[xaxis[0]*auoftime/1e-12, xaxis[-1]*auoftime/1e-12, yaxis[0]*27.211, yaxis[-1]*27.211])
    plt.title(outname+" ionization")
    plt.ylabel("Electron kin energy (eV)")
    plt.xlabel("Time delay (ps)")
    plt.savefig(outname+"Spectrum_nocbar.png")
    plt.colorbar(im)
    plt.savefig(outname+"Spectrum.png")
    plt.close()

    
    im = plt.imshow(np.abs(spec_fft[::-1,:nlim]), aspect=resfft, extent=[freqs[0], freqs[nlim], yaxis[0]*27.211, yaxis[-1]*27.211], interpolation=None)
    plt.title(outname+" ionization (fft)")
    plt.ylabel("Electron kin energy (eV)")
    plt.xlabel("freqs (eV)")
    plt.savefig(outname+"Spectrum_fftshifted_nocbar.png")
    plt.colorbar(im)
    plt.savefig(outname+"Spectrum_fftshifted.png")
    plt.close()
    return 0

def make_the_plots_first_timeprop(name12, name32, outname, gam):
    smto = 4*gam*np.sqrt(np.log(2))
    ncon, Elo, dE = np.loadtxt("Python_Code_RadInts/EcontEV_det.dat")
    Esam= np.arange(Elo, Elo+ncon*dE, dE)/CMperAU
    Esam12 = Esam - IP12
    Esam32 = Esam - IP32
    yaxis= np.arange(min(Esam12[0],Esam32[0]),max(Esam12[-1],Esam32[-1]),Esam[1]-Esam[0])

    #Load the one photon data
    data32_op = np.loadtxt(name32+".dat")
    data12_op = np.loadtxt(name12+".dat")
    xaxis = np.linspace(smto, 4.0e-12/auoftime,1000)[0:840]
    tdels = np.linspace(smto, 4.0e-12/auoftime,1000)[0:840]
    print("First")
    print(tdels[0]*auoftime/1e-12,tdels[-1]*auoftime/1e-12)
    probs32 = np.transpose(data32_op[:,:])
    probs12 = np.transpose(data12_op[:,:])
    yaxis, xaxis, spec= one_spectrum(Esam, tdels, probs12, probs32)
    print("Second")
    print(xaxis[0]*auoftime/1e-12,xaxis[-1]*auoftime/1e-12)
    np.savetxt(outname+".dat", spec)
    np.savetxt(outname+"xaxis.dat", xaxis)
    np.savetxt(outname+"yaxis.dat", yaxis)

    res = (xaxis[-1]*auoftime/1e-12-xaxis[0]*auoftime/1e-12)/(yaxis[-1]*27.211-yaxis[0]*27.211)

    freqs, spec_fft = fft_spectrum_shifted(xaxis*auoftime/1e-12, spec)
    nlim = int(0.1/(freqs[1]-freqs[0]))
    print(nlim, len(freqs))
    resfft = (freqs[nlim] - freqs[0])/(yaxis[-1]*27.211-yaxis[0]*27.211)
    print(resfft)

    plt.figure()
    ax = plt.gca()
    print(xaxis[0]*auoftime/1e-12, xaxis[-1]*auoftime/1e-12)
    im = plt.imshow(spec[::-1,:], aspect=res, extent=[xaxis[0]*auoftime/1e-12, xaxis[-1]*auoftime/1e-12, yaxis[0]*27.211, yaxis[-1]*27.211], interpolation=None)
    plt.title(outname+" ionization")
    plt.ylabel("Electron kin energy (eV)")
    plt.xlabel("Time delay (ps)")
    plt.savefig(outname+"Spectrum_nocbar.png")
    plt.colorbar(im)
    plt.savefig(outname+"Spectrum.png")
    plt.close()

    plt.figure()
    im = plt.imshow(np.abs(spec_fft[::-1,:nlim]), aspect=resfft, extent=[freqs[0], freqs[nlim], yaxis[0]*27.211, yaxis[-1]*27.211], interpolation='none')
    plt.title(outname+" ionization (fft)")
    plt.ylabel("Electron kin energy (eV)")
    plt.xlabel("freqs (eV)")
    plt.savefig(outname+"Spectrum_fftshifted_nocbar.png")
    plt.colorbar(im)
    plt.savefig(outname+"Spectrum_fftshifted.png")
    plt.close()
    return 0

def make_the_plots(name, exp_flag=False):
    ncon, Elo, dE = np.loadtxt("Python_Code_RadInts/EcontEV_det.dat")
    Esam= np.arange(Elo, Elo+ncon*dE, dE)/CMperAU
    Esam12 = Esam - IP12
    Esam32 = Esam - IP32
    yaxis= np.arange(min(Esam12[0],Esam32[0]),max(Esam12[-1],Esam32[-1]),Esam[1]-Esam[0])
    custom_lines = [Line2D([0], [0], color="C0", lw=4),
                Line2D([0], [0], color="C1", lw=4),Line2D([0], [0], color="C2", lw=4),
                Line2D([0], [0], color="C3", lw=4),Line2D([0], [0], color="C4", lw=4),
                Line2D([0], [0], color="C5", lw=4)]

    #Load the one photon data
    spec = np.loadtxt(name+".dat")
    xaxis = np.loadtxt(name+"xaxis.dat")
    yaxis = np.loadtxt(name+'yaxis.dat')

    res = (xaxis[-1]*auoftime/1e-12-xaxis[0]*auoftime/1e-12)/(yaxis[-1]*27.211-yaxis[0]*27.211)

    freqs, spec_fft = fft_spectrum_shifted(xaxis*auoftime/1e-12, spec)
    nlim = int(0.25/(freqs[1]-freqs[0]))

    resfft = (freqs[nlim] - freqs[0])/(yaxis[-1]*27.211-yaxis[0]*27.211)


    #Spectrogram construction
    plt.figure()
    ax = plt.gca()
    im = plt.imshow(spec[::-1,:], aspect=res, extent=[xaxis[0]*auoftime/1e-12, xaxis[-1]*auoftime/1e-12, yaxis[0]*27.211, yaxis[-1]*27.211], interpolation='none')
    plt.title(name+" ionization")
    plt.ylabel("Electron kin energy (eV)")
    plt.xlabel("Time delay (ps)")
    plt.savefig(name+"Spectrum_nocbar.png")
    plt.colorbar(im)
    #plt.tight_layout()
    plt.savefig(name+"Spectrum.png")
    plt.close()

    plt.figure()
    ax = plt.gca()
    im = plt.imshow(np.abs(spec_fft[::-1,:nlim]), aspect=resfft, extent=[freqs[0], freqs[nlim], yaxis[0]*27.211, yaxis[-1]*27.211], interpolation=None)
    plt.title(name+" ionization (fft)")
    plt.ylabel("Electron kin energy (eV)")
    plt.xlabel("freqs (eV)")
    if(exp_flag):
        print(Ebs_exp)
        plt.axvline(27.211*abs(Ebs_exp[0]-Ebs_exp[1]),color='C0',alpha=1.0,linestyle=":")
        plt.axvline(27.211*abs(Ebs_exp[2]-Ebs_exp[3]),color='C1',alpha=1.0,linestyle=":")
        plt.axvline(27.211*abs(Ebs_exp[0]-Ebs_exp[2]),color='C2',alpha=1.0,linestyle=":")
        plt.axvline(27.211*abs(Ebs_exp[1]-Ebs_exp[3]),color='C3',alpha=1.0,linestyle=":")
        plt.axvline(27.211*abs(Ebs_exp[0]-Ebs_exp[3]),color='C4',alpha=1.0,linestyle=":")
        plt.axvline(27.211*abs(Ebs_exp[1]-Ebs_exp[2]),color='C5',alpha=1.0,linestyle=":")
    else:
        print(Ebs)
        plt.axvline(27.211*abs(Ebs[0]-Ebs[1]),color='C0',alpha=0.2)
        plt.axvline(27.211*abs(Ebs[2]-Ebs[3]),color='C1',alpha=0.2)
        plt.axvline(27.211*abs(Ebs[0]-Ebs[2]),color='C2',alpha=0.2)
        plt.axvline(27.211*abs(Ebs[1]-Ebs[3]),color='C3',alpha=0.2)
        plt.axvline(27.211*abs(Ebs[0]-Ebs[3]),color='C4',alpha=0.2)
        plt.axvline(27.211*abs(Ebs[1]-Ebs[2]),color='C5',alpha=0.2)
    ax.legend(custom_lines, ['3d-5s',"3d'-5s'","5s-5s'","3d-3d'","5s-3d'","3d-5s'"])
    plt.savefig(name+"Spectrum_fftshifted_nocbar.png")
    #
    # divider = make_axes_locatable(ax)
    plt.colorbar(im)
    #plt.tight_layout()
    plt.savefig(name+"Spectrum_fftshifted.png")
    plt.close()
    return 0

def OP_vs_TP_prof(dirw):
    #Load the spectrongrams and average at each energy
    OP = np.laodtxt(dirw+"/One Photon.dat")
    OPxaxis = np.loadtxt(dirw+"/One Photonxaxis.dat")
    OPyaxis = np.loadtxt(dirw+"/One Photonyaxis.dat")
    #The two photon data
    TP = np.laodtxt(dirw+"/Two Photon.dat")
    TPxaxis = np.loadtxt(dirw+"/Two Photonxaxis.dat")
    TPyaxis = np.loadtxt(dirw+"/Two Photonyaxis.dat")

    #Make a simple plot just to try
    fig, ax = plt.subplots(1,1)
    ax.plot(OPyaxis, OP[:,0])
    ax.plot(TPyaxis, TP[:,0])
    plt.savefig(dirw+"/Profile_Comp.png", dpi=150)
    
    return 0

def lineouts(dirw, bf):
    #Load the one photon data
    specop = np.loadtxt(dirw+"One Photon.dat")
    xaxis = np.loadtxt(dirw+"One Photonxaxis.dat")
    yaxis = np.loadtxt(dirw+'One Photonyaxis.dat')
    freqs, specop_fft = fft_spectrum_shifted(xaxis*auoftime/1e-12, specop)
    #Getting the lineouts.                                                                                                                                                                        
    fbin =  np.argmin(np.abs(freqs-bf))
    yim = specop_fft[::,fbin]
    fig, axs = plt.subplots(1,1)
    ax2 = axs.twinx()
    axs.set_xlabel("Kinetic Energy (eV)")
    axs.set_ylabel("Fourier amplitude")
    axs.set_title("Lineout at beat at %.3e eV"%bf)
    ax2.set_ylabel("Phase /$\\pi$")
    ax2.set_ylim(0,2)
    axs.plot(yaxis*27.211, np.abs(yim))
    ax2.plot(yaxis*27.211, np.mod(np.unwrap(np.angle(yim))/np.pi,2),"r--", alpha=0.6)
    plt.savefig(dirw+'OnePhoton_lineout.png')
    plt.close()
    #Load the one photon data
    spectp = np.loadtxt(dirw+"Two Photon.dat")
    xaxis = np.loadtxt(dirw+"Two Photonxaxis.dat")
    yaxis = np.loadtxt(dirw+'Two Photonyaxis.dat')
    freqs, spectp_fft = fft_spectrum_shifted(xaxis*auoftime/1e-12, spectp)
    #Getting the lineouts.
    fbin =  np.argmin(np.abs(freqs-bf))
    yim = spectp_fft[::,fbin]
    fig, axs = plt.subplots(1,1)
    ax2 = axs.twinx()
    axs.set_xlabel("Kinetic Energy (eV)")
    axs.set_ylabel("Fourier amplitude")
    axs.set_title("Lineout at beat at %.3e eV"%bf)
    ax2.set_ylabel("Phase /$\\pi$")
    ax2.set_ylim(0,2)
    axs.plot(yaxis*27.211, np.abs(yim))
    ax2.plot(yaxis*27.211, np.mod(np.unwrap(np.angle(yim))/np.pi,2),"r--", alpha=0.6)
    plt.savefig(dirw+'TwoPhoton_lineout.png')
    plt.close()
    return 0

def lineouts_tp(dirw,bf):
    #Load the one photon data
    spectp = np.loadtxt(dirw+"Two Photon.dat")
    xaxis = np.loadtxt(dirw+"Two Photonxaxis.dat")
    yaxis = np.loadtxt(dirw+'Two Photonyaxis.dat')
    freqs, spectp_fft = fft_spectrum_shifted(xaxis*auoftime/1e-12, spectp)
    #Getting the lineouts.
    fbin =  np.argmin(np.abs(freqs-bf))
    yim = spectp_fft[::,fbin]
    fig, axs = plt.subplots(1,1)
    ax2 = axs.twinx()
    axs.set_xlabel("Kinetic Energy (eV)")
    axs.set_ylabel("Fourier amplitude")
    axs.set_title("Lineout at beat")
    ax2.set_ylabel("Phase /$\\pi$")
    axs.plot(yaxis*27.211, np.abs(yim))
    ax2.plot(yaxis*27.211, np.angle(yim)/np.pi,"r--", alpha=0.6)
    plt.savefig(dirw+'TwoPhoton_lineout.png')
    plt.close()
    return 0

def lineouts_tp_fullname(dirw,bf):
    #Load the one photon data
    spectp = np.loadtxt(dirw+".dat")
    xaxis = np.loadtxt(dirw+"xaxis.dat")
    yaxis = np.loadtxt(dirw+'yaxis.dat')
    freqs, spectp_fft = fft_spectrum_shifted(xaxis*auoftime/1e-12, spectp)
    #Getting the lineouts.
    fbin =  np.argmin(np.abs(freqs-bf))
    yim = spectp_fft[::,fbin]
    fig, axs = plt.subplots(1,1)
    ax2 = axs.twinx()
    axs.set_xlabel("Kinetic Energy (eV)")
    axs.set_ylabel("Fourier amplitude")
    axs.set_title("Lineout at beat")
    ax2.set_ylabel("Phase /$\\pi$")
    axs.plot(yaxis*27.211, np.abs(yim))
    ax2.plot(yaxis*27.211, np.unwrap(np.angle(yim))/np.pi,"r--", alpha=0.6)
    plt.savefig(dirw+'TwoPhoton_lineout.png')
    plt.close()
    return 0

def probability_plots(dirw):
    OnePh = np.loadtxt(dirw+"One Photon.dat")
    OnePhx = np.loadtxt(dirw+"One Photonxaxis.dat")
    OnePhy = np.loadtxt(dirw+"One Photonyaxis.dat")
    OnePhprob = np.zeros_like(OnePhx, dtype=float)

    TwoPh = np.loadtxt(dirw+"Two Photon.dat")
    TwoPhx = np.loadtxt(dirw+"Two Photonxaxis.dat")
    TwoPhy = np.loadtxt(dirw+"Two Photonyaxis.dat")
    TwoPhprob = np.zeros_like(TwoPhx, dtype=float)

    for t in range(1,len(OnePhx)):
        OnePhprob[t] = trapz(y=OnePh[:,t], x=OnePhy)
    

    for t in range(len(TwoPhx)):
        TwoPhprob[t] = trapz(y=TwoPh[:,t], x=TwoPhy)
    
    plt.plot(OnePhx[1:], OnePhprob[1:], label='One Photon')
    plt.plot(TwoPhx, TwoPhprob, label='Two Photon')
    plt.title(dirw+"Photionization probability")
    plt.xlabel("Time Delay")
    plt.ylabel("Probability of photionization")
    plt.legend(loc=0)
    plt.savefig(dirw+"Probability_plots.png")
    plt.close()


    return 0

def make_SI_line_dat(nl,dirw):
    print("J=0 states: ")
    file0 = open("IndividualInterms/"+dirw+"/J0Int2p12.dat",'w')
    file1 = open("IndividualInterms/"+dirw+"/J0Int2p32.dat",'w')
    sp1 = [0,0,0]
    sp2 = [0,0,0]
    for i in range(1,14):
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J0p%i-%i_tdel_probsJ12.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        file0.write(format_me(xav))
        if(max(*xav) > sp1[0]):
            sp1 = [max(*xav), i, np.argmax(xav)]
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J0p%i-%i_tdel_probsJ32.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        file1.write(format_me(xav))
        if(max(*xav) > sp2[0]):
            sp2 = [max(*xav), i, np.argmax(xav)]
    file0.close()
    file1.close()
    print("Ionization to p12 is dominated by %i intermediate state with peak density %.8f for energy %i"%(sp1[1],sp1[0], sp1[2]))
    print("Ionization to p32 is dominated by %i intermediate state with peak density %.8f for energy %i"%(sp2[1],sp2[0], sp2[2]))
    print("J=2 f states: ")
    file0 = open("IndividualInterms/"+dirw+"/J2fInt2p12.dat",'w')
    file1 = open("IndividualInterms/"+dirw+"/J2fInt2p32.dat",'w')
    sp1 = [0,0]
    sp2 = [0,0]
    for i in range(1,25):
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J2f%i-%i_tdel_probsJ12.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        file0.write(format_me(xav))
        if(max(*xav) > sp1[0]):
            sp1 = [max(*xav), i, np.argmax(xav)]
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J2f%i-%i_tdel_probsJ32.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        file1.write(format_me(xav))
        if(max(*xav) > sp2[0]):
            sp2 = [max(*xav), i, np.argmax(xav)]
    file0.close()
    file1.close()
    print("Ionization to p12 is dominated by %i intermediate state with peak density %.8f at continuum energy %i"%(sp1[1],sp1[0],sp1[2]))
    print("Ionization to p32 is dominated by %i intermediate state with peak density %.8f at continuum energy %i"%(sp2[1],sp2[0],sp2[2]))
    
    print("J=2 p states: ")
    file0 = open("IndividualInterms/"+dirw+"/J2pInt2p12.dat",'w')
    file1 = open("IndividualInterms/"+dirw+"/J2pInt2p32.dat",'w')
    sp1 = [0,0]
    sp2 = [0,0]
    for i in range(1,25):
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J2p%i-%i_tdel_probsJ12.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        file0.write(format_me(xav))
        if(max(*xav) > sp1[0]):
            sp1 = [max(*xav), i, np.argmax(xav)]
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J2p%i-%i_tdel_probsJ32.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        file1.write(format_me(xav))
        if(max(*xav) > sp2[0]):
            sp2 = [max(*xav), i, np.argmax(xav)]
    file0.close()
    file1.close()
    print("Ionization to p12 is dominated by %i intermediate state with peak density %.8f at energy %i"%(sp1[1],sp1[0],sp1[2]))
    print("Ionization to p32 is dominated by %i intermediate state with peak density %.8f at energy %i"%(sp2[1],sp2[0],sp2[2]))
    return "Files written"

def make_SI_lineouts(nl,dirw):
    fig, ax = plt.subplots(1,3,sharey=False,gridspec_kw={'width_ratios':[1,1,0.05]})
    colors = pl.cm.jet(np.linspace(0,1,13))
    ncon, Elo, dEcont = np.loadtxt("Python_Code_RadInts/EcontEV_det.dat")
    yaxis = np.arange(Elo,Elo+ncon*dEcont, dEcont)/CMperAU
    plt.suptitle("%i run interm J=0 state contribution"%nl)
    sp1 = 0
    sp2 = 0
    for i in range(1,14):
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J0p%i-%i_tdel_probsJ12.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        ax[0].plot(27.211*(yaxis-IP12),xav, c=colors[i-1])
        ax[0].set_yscale("squareroot")
        sp1 = max(sp1,*xav)
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J0p%i-%i_tdel_probsJ32.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        ax[1].plot(27.211*(yaxis-IP32),xav,c=colors[i-1])
        ax[1].set_yscale("squareroot")
        sp2 = max(sp2,*xav)
    ax[0].set_title("$3p_{1/2}$ core channels")
    ax[1].set_title("$3p_{3/2}$ core channels")
    ax[0].set_ylim(0,1.1*max(sp1,sp2))
    ax[1].set_ylim(0,1.1*max(sp1,sp2))
    ax[1].set_yticks([])
    cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=1, vmax=13,clip=True)

    cb1 = mpl.colorbar.ColorbarBase(ax[2], cmap=cmap,
                                norm=norm,
                                orientation='vertical',boundaries=range(1,15))
    cb1.set_label('Interm State')
    plt.savefig("IndividualInterms/%iSIJ0p_profiles.png"%(nl),dpi=150)
    
    plt.close()
    
    colors = pl.cm.jet(np.linspace(0,1,24))
    fig1, ax1 = plt.subplots(1,3,sharey=False,gridspec_kw={'width_ratios':[1,1,0.05]})
    plt.suptitle("%i run interm J=2p state contribution"%nl)
    sp1 = 0
    sp2 = 0
    for i in range(1,25):
        #Plotting the 2p states
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J2p%i-%i_tdel_probsJ12.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        ax1[0].plot(27.211*(yaxis-IP12),xav, c=colors[i-1])
        ax1[0].set_yscale("squareroot")
        sp1 = max(sp1,*xav)
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J2p%i-%i_tdel_probsJ32.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        ax1[1].plot(27.211*(yaxis-IP32),xav,c=colors[i-1])
        ax1[1].set_yscale("squareroot")
        sp2 = max(sp2,*xav)
    ax1[0].set_title("$3p_{1/2}$ core channels")
    ax1[1].set_title("$3p_{3/2}$ core channels")
    ax1[0].set_ylim(0,1.1*max(sp1,sp2))
    ax1[1].set_ylim(0,1.1*max(sp1,sp2))
    ax1[1].set_yticks([])
    cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=1, vmax=24,clip=True)

    cb1 = mpl.colorbar.ColorbarBase(ax1[2], cmap=cmap,
                                norm=norm,
                                orientation='vertical',boundaries=range(1,26))
    cb1.set_label('Interm State')
    plt.savefig("IndividualInterms/%iSIJ2p_profiles.png"%(nl),dpi=150)
    plt.close()

    fig1, ax1 = plt.subplots(1,3,sharey=False,gridspec_kw={'width_ratios':[1,1,0.05]})
    plt.suptitle("%i run interm J=2f state contribution"%nl)
    sp1 = 0
    sp2 = 0
    for i in range(1,25):
        #Plotting the 2p states
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J2f%i-%i_tdel_probsJ12.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        ax1[0].plot(27.211*(yaxis-IP12),xav, c=colors[i-1])
        ax1[0].set_yscale("squareroot")
        sp1 = max(sp1,*xav)
        x = np.loadtxt("IndividualInterms/"+dirw+"/'%iSI_J2f%i-%i_tdel_probsJ32.dat"%(nl,i,i))
        ts = np.array([xx[0] for xx in x])
        xav = np.zeros_like(x[0,1:])
        for k in range(1,len(x[0])):
            temp = np.array([xx[k] for xx in x])
            xav[k-1] = trapz(y=temp, x=ts)/(ts[-1]-ts[0])
        ax1[1].plot(27.211*(yaxis-IP32),xav,c=colors[i-1])
        ax1[1].set_yscale("squareroot")
        sp2 = max(sp2,*xav)
    ax1[0].set_title("$3p_{1/2}$ core channels")
    ax1[1].set_title("$3p_{3/2}$ core channels")
    ax1[0].set_ylim(0,1.1*max(sp1,sp2))
    ax1[1].set_ylim(0,1.1*max(sp1,sp2))
    ax1[1].set_yticks([])
    cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=1, vmax=24,clip=True)

    cb1 = mpl.colorbar.ColorbarBase(ax1[2], cmap=cmap,
                                norm=norm,
                                orientation='vertical',boundaries=range(1,26))
    cb1.set_label('Interm State')
    plt.savefig("IndividualInterms/%iSIJ2f_profiles.png"%(nl), dpi=150)

def initiate_plots(data12_op, data32_op, bfs, axs):
    ax_spec,ax_fft, ax_bs = axs[0],axs[1],axs[2:]
    ncon, Elo, dE = np.loadtxt("Python_Code_RadInts/EcontEV_det.dat")
    Esam= np.arange(Elo, Elo+ncon*dE, dE)/CMperAU
    Esam12 = Esam - IP12
    Esam32 = Esam - IP32
    yaxis= np.arange(min(Esam12[0],Esam32[0]),max(Esam12[-1],Esam32[-1]),Esam[1]-Esam[0])

    
    xaxis, probs32 = data32_op[:,0], np.transpose(data32_op[:,1:])
    tdels, probs12 = data12_op[:,0], np.transpose(data12_op[:,1:])
    yaxis, xaxis, spec= one_spectrum(Esam, tdels, probs12, probs32)

    line_ffts = np.array([[None, None, None]]*len(bfs),dtype=object)

    freqs, spec_fft, mean = fft_spectrum_meaned(xaxis*auoftime/1e-12, spec)
    fm = max([1.01*bf[0] for bf in bfs])
    nlim = int(fm/(freqs[1]-freqs[0]))
    if(nlim>=len(freqs)):
        raise ValueError("Increase the sampling frequency.")

    res = 0.5*(xaxis[-1]*auoftime/1e-12-xaxis[0]*auoftime/1e-12)/(yaxis[-1]*27.211-yaxis[0]*27.211)
    im = ax_spec.imshow(spec[::-1,:], aspect=res, extent=[xaxis[0]*auoftime/1e-12, xaxis[-1]*auoftime/1e-12, yaxis[0]*27.211, yaxis[-1]*27.211])
    
    resfft = 0.5*(freqs[nlim]-freqs[0])/(yaxis[-1]*27.211-yaxis[0]*27.211)
    im2 = ax_fft.imshow(np.abs(spec_fft[::-1,:nlim])**1/8.,aspect=resfft, extent=[freqs[0],freqs[nlim],yaxis[0]*27.211, yaxis[-1]*27.211])
    for i in range(len(bfs)):
        line_ffts[i,2] = ax_fft.axvline(bfs[i,0],c='C%i'%i)
        if(i%2==0):
            ax_fft.text(x=bfs[i,0], y=yaxis[-1]*27.211,s="$%s$"%bfs[i,1],c='k',rotation=90,fontsize=8,ha='center',va='bottom')
        else:
            ax_fft.text(x=bfs[i,0], y=yaxis[0]*27.211,s="$%s$"%bfs[i,1],c='w',rotation=90,fontsize=8,ha='center',va='bottom')
    custom_lines = [Line2D([0], [0], color='C%i'%i, lw=1) for i in range(len(bfs))]
    #ax_fft.legend(custom_lines, bfs[:,1],bbox_to_anchor=(1.0,0.5)) 

    nbs = np.array([int(bf[0]/(freqs[1]-freqs[0])) for bf in bfs])
    
    nrm  = max([max(np.abs(spec_fft[:,nbs[i]])) for i in range(len(bfs))])
    ax_bs[-1].set_title('Fourier Amplitudes')
    ax_bs[-2].set_title('Fourier Phase')
    line_mean, = ax_bs[-1].plot(yaxis*27.211,np.abs(mean/max(np.abs(mean))),c='k')
    for i in range(len(bfs)):
        ax_bs[-1].set_ylim(-1e-3,1.1)
        #if(i==0): ax2 = ax_bs[-1].twinx()
        phs = np.unwrap(np.angle(spec_fft[:,nbs[i]]))/np.pi
        sig = np.abs(spec_fft[:,nbs[i]])
        if(max(sig)/nrm>0.2):
            line_ffts[i,0],  = ax_bs[-1].plot(yaxis*27.211, sig/nrm,label="%s"%(bfs[i,1]),c="C%i"%i,alpha=1.0)
            line_ffts[i,1],  = ax_bs[-2].plot(yaxis*27.211, np.mod(phs-phs[0],2),'--',c="C%i"%i, alpha=1.0)
            line_ffts[i,2].set(alpha=1.0)
        else:
            line_ffts[i,0],  = ax_bs[-1].plot(yaxis*27.211, sig/nrm,label="%s"%(bfs[i,1]),c="C%i"%i,alpha=0.0)
            line_ffts[i,1],  = ax_bs[-2].plot(yaxis*27.211, np.mod(phs-phs[0],2),'--',c="C%i"%i, alpha=0.0)
            line_ffts[i,2].set(alpha=0.0)
        #ax2.set_ylim(-1,1)
        if(i==len(bfs)-1): 
            ax_bs[-1].set_xlabel("Kinetic Energy (eV)",fontsize=15)
            ax_bs[-2].set_ylabel("Phase / $\\pi$")

    return im,im2,line_ffts, line_mean

def update_plots(data12_op, data32_op, ax_spec, ax_fft, line_ffts, bfs, line_mean):
    ncon, Elo, dE = np.loadtxt("Python_Code_RadInts/EcontEV_det.dat")
    Esam= np.arange(Elo, Elo+ncon*dE, dE)/CMperAU
    
    xaxis, probs32 = data32_op[:,0], np.transpose(data32_op[:,1:])
    tdels, probs12 = data12_op[:,0], np.transpose(data12_op[:,1:])
    yaxis, xaxis, spec= one_spectrum(Esam, tdels, probs12, probs32)

    freqs, spec_fft,mean = fft_spectrum_meaned(xaxis*auoftime/1e-12, spec)
    fm = max([1.01*bf[0] for bf in bfs])
    nlim = int(fm/(freqs[1]-freqs[0]))
    

    
    ax_spec.set_data(spec[::-1,:])
    ax_fft.set_data(np.abs(spec_fft[::-1,:nlim])**1/8.)
    
    nbs = np.array([int(bf[0]/(freqs[1]-freqs[0])) for bf in bfs])
    nrm  = max([max(np.abs(spec_fft[:,nbs[i]])) for i in range(len(bfs))])
    line_mean.set_ydata(mean/np.max(mean))
    for i in range(len(bfs)):
        sig = np.abs(spec_fft[:,nbs[i]])
        if(max(sig)/nrm>0.2):
            line_ffts[i,0].set_ydata(np.abs(spec_fft[:,nbs[i]])/nrm)
            line_ffts[i,0].set(alpha=1.0)
            phs = np.unwrap(np.angle(spec_fft[:,nbs[i]]))/np.pi
            line_ffts[i,1].set_ydata(np.mod(phs-phs[0],2))
            line_ffts[i,1].set(alpha=1.0)
            line_ffts[i,2].set(alpha=1.0)

        else:
            line_ffts[i,0].set_ydata(np.abs(spec_fft[:,nbs[i]])/nrm)
            line_ffts[i,0].set(alpha=0.0)
            phs = np.unwrap(np.angle(spec_fft[:,nbs[i]]))/np.pi
            line_ffts[i,1].set_ydata(np.mod(phs-phs[0],2))
            line_ffts[i,1].set(alpha=0.0)
            line_ffts[i,2].set(alpha=0.0)
    
    return 0

def initiate_resol(data12_op,data32_op,axs):
    ax12, ax32 = axs[0], axs[1]
    ncon, Elo, dE = np.loadtxt("Python_Code_RadInts/EcontEV_det.dat")
    Esam= np.arange(Elo, Elo+ncon*dE, dE)/CMperAU
    Js = ["J=1","J=3","J=5"]
    l12 = [None]*3
    l32 = [None]*3
    Esam12 = Esam - IP12
    Esam32 = Esam - IP32
    for i in range(3):
        l12[i], = ax12.plot(Esam12*27.211, data12_op[i],label=Js[i],c="C%i"%i)
        l32[i], = ax32.plot(Esam32*27.211, data32_op[i],label=Js[i],c="C%i"%i)
    ax12.legend(loc=0)
    ax12.set_ylim(-1e-2,0.3)
    ax32.legend(loc=0)
    ax32.set_ylim(-1e-2,0.3)
    return l12, l32

def initiate_cResol(data12_op,data32_op,axs):
    j12chan = ["nd J1",'ns J1', 'nd J3', 'ng J3', 'ng J5', 'ni J5']
    j32chan = ["nd J1 (jcs=1)",'nd J1 (jcs=2)', 'ns J1 (jcs=1)',
               'nd J3 (jcs=1)','nd J3 (jcs=2)', 'ng J3 (jcs=1)','ng J3 (jcs=2)',
               'ng J5 (jcs=1)','ng J5 (jcs=2)', 'ni J5 (jcs=1)','ni J5 (jcs=1)']
    ax12, ax32 = axs[0], axs[1]
    ncon, Elo, dE = np.loadtxt("Python_Code_RadInts/EcontEV_det.dat")
    Esam= np.arange(Elo, Elo+ncon*dE, dE)/CMperAU
    Esam12 = Esam - IP12
    Esam32 = Esam - IP32
    l12 = [None]*6
    l32 = [None]*11
    for i in range(6):
        if(i<2):
            l12[i], = ax12.plot(Esam12*27.211,data12_op[i],linestyle='--',label=j12chan[i])
        elif(i<4):
            l12[i], = ax12.plot(Esam12*27.211,data12_op[i],linestyle='--',label=j12chan[i])
        else:
            l12[i], = ax12.plot(Esam12*27.211,data12_op[i],linestyle='--',label=j12chan[i])
    ax12.legend(loc=0)
    ax12.set_ylim(-1e-2,0.5)
    for i in range(11):
        if(i<3):
            l32[i], = ax32.plot(Esam32*27.211,data32_op[i],linestyle='--',label=j32chan[i])
        elif(i<7):
            l32[i], = ax32.plot(Esam32*27.211,data32_op[i],linestyle='--',label=j32chan[i])
        else:
            l32[i], = ax32.plot(Esam32*27.211,data32_op[i],linestyle='--',label=j32chan[i])
    ax32.legend(loc=0)
    ax32.set_ylim(-1e-2, 0.5)
    return l12, l32



def update_resol(data12_op, data32_op, l12,l32):
    for i in range(len(l12)):
        l12[i].set_ydata(data12_op[i])
    for i in range(len(l32)):
        l32[i].set_ydata(data32_op[i])

