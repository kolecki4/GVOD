#!/usr/bin/env python
# coding: utf-8

from Settings import *
from scipy.signal import peak_widths
from scipy.optimize import curve_fit
from scipy.special import wofz
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# IF THERE ARE BLENDED LINES, SET THIS TO TRUE. IF THE LINES ARE FAIRLY ISOLATED, SET THIS TO FALSE
lotsofblendedlines = True

# Defines
def gaussian(x, a1, mu1, sig1,background):
    return background-a1*np.exp(-np.power(x - mu1, 2.) / (2 * np.power(sig1, 2.)))

def Voigt(x, x0, y0, a, sigma, gamma,background):
    return background - a * np.real(wofz((x - x0 + 1j*gamma)/sigma/np.sqrt(2))) / sigma /np.sqrt(2*np.pi)

# Finds where the spectral lines of interest are in the observed spectrum
indices = np.asarray([],dtype = int)
for i in linelist['Wavelength']:
    if min(wav) <= i and  i <= max(wav):
        indices = np.append(indices,(np.abs(wav - i)).argmin())
        
# A list of the line wavelengths in the observed spectrum
importantwavs = np.asarray([], dtype = float)
for i in indices:
    importantwavs = np.append(importantwavs, wav[i])

adjustedindices = ([flx[i-halfwindow:i+halfwindow].argmin()+(i-halfwindow) for i in indices])

# I used the FWHM calculations to define the window within which the eqws were calculated
FWHM = peak_widths(1-flx, rel_height = FWHeight, peaks = adjustedindices)[0]
FWHMS = np.array([int(round(i,0)) for  i in FWHM])

# Initializes the eqws array
eqws = np.asarray([], dtype = float)

# This loop takes the data surrounding each peak, calculates the background level,
# and calculates the eqw based on the formula from the wikipedia page on equivalent width
# where the integral is replaced with a Riemann sum
for i in range(len(adjustedindices)):
    #Prepare for plotting
    #plt.figure(dpi=200)
    
    # Gets relevant wavelength and flux data
    yesflx = flx[adjustedindices[i]-(2*FWHMS[i]):adjustedindices[i]+(2*FWHMS[i])] 
    yeswav = wav[adjustedindices[i]-(2*FWHMS[i]):adjustedindices[i]+(2*FWHMS[i])]
    
    #Plot the observed data
    if len(yesflx) > 0:
        maxflx = max(yesflx)
        minflx = min(yesflx)
    else:
        maxflx = max(flx)
        minflx = min(flx)
        
    if method == 'obs':
        if not lotsofblendedlines:
            # takes the continuum to be the average of the points toward the edge of the area of interest
            background = np.average( np.append( yesflx[int(len(yesflx)/2- 1.5*FWHMS[i]-5):int(len(yesflx)/2- 1.5*FWHMS[i])-0], yesflx[int(len(yesflx)/2 +1.5*FWHMS[i]+0):int(len(yesflx)/2 +1.5*FWHMS[i]+5)])) 
        else:
            # takes the continuum to be the average of all flux data > 0.9
            background = np.average(yesflx[yesflx >= 0.9*maxflx])
        
        EW = 0
        for j in range(len(yesflx)-1):
            EW += (1 - (yesflx[j]/background))*(yeswav[j+1] - yeswav[j]) 
        eqws = np.append(eqws,EW)
        
        plt.plot(yeswav,yesflx, color ='#666666', zorder = 2)
        plt.vlines(wav[adjustedindices[i]],0,1.1, zorder = 1, color = '#7B8fB9')
        plt.vlines(wav[adjustedindices[i]]-(eqws[i]/2),0,background, zorder = 1, color = '#041E42')
        plt.vlines(wav[adjustedindices[i]]+(eqws[i]/2),0,background, zorder = 1, color = '#041E42')
        if len(yeswav) > 0: plt.hlines(background,yeswav[0], yeswav[-1], color = '#041E42')
        plt.xticks(fontsize = '8')
        #plt.ylim(max([flx[adjustedindices[i]] - 0.5*(1-flx[adjustedindices[i]]),0]),1.025*maxflx)
        try:
            plt.ylim(0,1.1*maxflx)
        except:
            pass
        idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
        plt.title(linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B (EW = ' +'{:.3f}'.format(1000*EW) + ' m\u212B)')
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel('Normalized Flux')
        plt.show()
    
    if method == 'gau':
        EW = 0
        OEW = 0
        try:
            popt2, pcov=curve_fit(gaussian,yeswav,yesflx,p0=[maxflx - minflx,wav[adjustedindices][i],.2,maxflx])
            background=popt2[3]
            for j in range(len(yeswav)-1):
                EW += (1 - (gaussian(yeswav,*popt2)[j]/background))*(yeswav[j+1] - yeswav[j]) 
                #print('Equivalent width from Gaussian model: ' + '{:.2f}'.format(1000*EW) + 'm\u212B')
            
            if not lotsofblendedlines:
            # takes the continuum to be the average of the points toward the edge of the area of interest
                Obackground = np.average( np.append( yesflx[int(len(yesflx)/2- 1.5*FWHMS[i]-5):int(len(yesflx)/2- 1.5*FWHMS[i])-0], yesflx[int(len(yesflx)/2 +1.5*FWHMS[i]+0):int(len(yesflx)/2 +1.5*FWHMS[i]+5)])) 
            else:
            # takes the continuum to be the average of all flux data > 0.9
                Obackground = np.average(yesflx[yesflx >= 0.9*maxflx])
            
            for j in range(len(yesflx)-1):
                OEW += (1 - (yesflx[j]/Obackground))*(yeswav[j+1] - yeswav[j]) 
                
                
            plt.plot(yeswav, gaussian(yeswav,*popt2),ls='--', c='#4444FF',lw=3, label='Model')
        except:
            print('Could not model with Gaussian.')
            
        if np.abs(OEW-EW) < .5*EW:
            plt.plot(yeswav,yesflx, color ='#666666', zorder = 2)    
            plt.vlines(wav[adjustedindices[i]],0,1.1, zorder = 1, color = '#7B8fB9')
            plt.vlines(wav[adjustedindices[i]]-(EW/2),0,background, zorder = 1, color = '#041E42')
            plt.vlines(wav[adjustedindices[i]]+(EW/2),0,background, zorder = 1, color = '#041E42')
            if len(yeswav) > 0: plt.hlines(background,yeswav[0], yeswav[-1], color = '#041E42')
            plt.xticks(fontsize = '8')
            #plt.ylim(max([flx[adjustedindices[i]] - 0.5*(1-flx[adjustedindices[i]]),0]),1.025*maxflx)
            try:
                plt.ylim(0,1.1*maxflx)
            except:
                pass
            idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
            plt.title(str(linelist['Element'][idk]) + ' ' + str(linelist['Ionization'][idk]) + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B (EW = ' +'{:.3f}'.format(1000*EW) + ' m\u212B)')
            plt.xlabel('Wavelength (Angstroms)')
            plt.ylabel('Normalized Flux')
            plt.show()
        else:
            EW = 0
        
        eqws = np.append(eqws,EW)
        
    if method == 'voi':
        try:
            popt1, pcov=curve_fit(Voigt,yeswav,yesflx,p0=[wav[adjustedindices][i],maxflx,maxflx - minflx,1,.1,maxflx]) 
            EW = 0
            background=popt1[5]
            for j in range(len(yeswav)-1):
                EW += (1 - (Voigt(yeswav,*popt1)[j]/background))*(yeswav[j+1] - yeswav[j])
            eqws = np.append(eqws,EW)
            plt.plot(yeswav, Voigt(yeswav,*popt1),ls='--', c='#FF4444',lw=3, label='Model')
            
            if not lotsofblendedlines:
            # takes the continuum to be the average of the points toward the edge of the area of interest
                Obackground = np.average( np.append( yesflx[int(len(yesflx)/2- 1.5*FWHMS[i]-5):int(len(yesflx)/2- 1.5*FWHMS[i])-0], yesflx[int(len(yesflx)/2 +1.5*FWHMS[i]+0):int(len(yesflx)/2 +1.5*FWHMS[i]+5)])) 
            else:
            # takes the continuum to be the average of all flux data > 0.9
                Obackground = np.average(yesflx[yesflx >= 0.9*maxflx])
            OEW = 0
            for j in range(len(yesflx)-1):
                OEW += (1 - (yesflx[j]/Obackground))*(yeswav[j+1] - yeswav[j]) 
            
        except:
            print('Could not model with Voigt distribution.')
            EW=-1
            OEW=-1
            
        if 0.5 < EW/OEW < 2:
            eqws = np.append(eqws,EW)
            plt.plot(yeswav,yesflx, color ='#666666', zorder = 2)    
            plt.vlines(wav[adjustedindices[i]],0,1.1, zorder = 1, color = '#7B8fB9')
            plt.vlines(wav[adjustedindices[i]]-(EW/2),0,background, zorder = 1, color = '#041E42')
            plt.vlines(wav[adjustedindices[i]]+(EW/2),0,background, zorder = 1, color = '#041E42')
            if len(yeswav) > 0: plt.hlines(background,yeswav[0], yeswav[-1], color = '#041E42')
            plt.xticks(fontsize = '8')
            #plt.ylim(max([flx[adjustedindices[i]] - 0.5*(1-flx[adjustedindices[i]]),0]),1.025*maxflx)
            try:
                plt.ylim(0,1.1*maxflx)
            except:
                pass
            idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
            plt.title(linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B (EW = ' +'{:.3f}'.format(1000*EW) + ' m\u212B)')
            plt.xlabel('Wavelength (Angstroms)')
            plt.ylabel('Normalized Flux')
            plt.show()
        else:
            eqws = np.append(eqws,-1)
    
    
    if method == 'ask':
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (12,4))
        idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
        fig.suptitle(linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
        
        counter = 1
        guesses = [wav[adjustedindices][i],maxflx,maxflx - minflx,1,.1,maxflx]        
        while counter == 1:
            #plt.figure()
            try:
                counter = 0
                VEW = 0
                
                popt1, pcov=curve_fit(Voigt,yeswav,yesflx,p0=guesses) 
                Vbackground=popt1[5]
                for j in range(len(yeswav)-1):
                    VEW += (1 - (Voigt(yeswav,*popt1)[j]/Vbackground))*(yeswav[j+1] - yeswav[j])
                #eqws = np.append(eqws,VEW)
                
                
                ax1.plot(yeswav,yesflx, color ='#666666', zorder = 2)
                
                #plt.title(linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
                ax1.set_xlabel('Wavelength (Angstroms)')
                ax1.set_ylabel('Flux')
                ax1.plot(yeswav, Voigt(yeswav,*popt1),ls='--', c='#FF4444',lw=3, label='Model')
                ax1.vlines(wav[adjustedindices[i]],0,1.1*Vbackground, zorder = 1, color = '#7B8fB9')
                ax1.vlines(wav[adjustedindices[i]]-(VEW/2),0,Vbackground, zorder = 1, color = '#FF4444')
                ax1.vlines(wav[adjustedindices[i]]+(VEW/2),0,Vbackground, zorder = 1, color = '#FF4444')
                ax1.set_title('Voigt Model EW = ' + '{:.3f}'.format(1000*VEW) + 'm\u212b')
                if len(yeswav) > 0: ax1.hlines(Vbackground,yeswav[0], yeswav[-1], color = '#FF4444')
                ax1.set_ylim(0,1.1*maxflx)
                    
            except:
                counter = 1
                ax1.plot(yeswav,yesflx, color ='#666666', zorder = 2)
                idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
                #plt.title(linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
                ax1.set_xlabel('Wavelength (Angstroms)')
                ax1.set_ylabel('Flux')
                ax1.set_title(linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
                #plt.show(block = False)
                
                # hlp, ax = plt.subplots()
                # ax.plot(yeswav,yesflx, color ='#666666', zorder = 2)
                # hlp.show()
                instring = ''
                if not autoDiscard:
                    print('Could not model with Voigt distribution. Enter new parameters')
                    while instring != 'p' and len(instring) < 11:
                        instring = input('Separated by spaces, Input the following.\n(x0, y0, a, sigma, gamma, background).\nOr enter "p" to skip this line: ')
                    if instring == "p":
                        VEW = -1
                        #eqws = np.append(eqws,EW)
                        break
                    guesses = [float(i) for i in instring.split()]
                else:
                    VEW = -1
                    counter=0
            
        counter = 1
        guesses = [maxflx - minflx,wav[adjustedindices][i],.2,maxflx]
        while counter == 1:
            #plt.figure()
            try:
                counter = 0
                GEW = 0
                
                popt2, pcov=curve_fit(gaussian,yeswav,yesflx,p0=guesses)
                Gbackground=popt2[3]
                for j in range(len(yeswav)-1):
                    GEW += (1 - (gaussian(yeswav,*popt2)[j]/Gbackground))*(yeswav[j+1] - yeswav[j]) 
                #eqws = np.append(eqws,GEW)
                ax2.plot(yeswav, gaussian(yeswav,*popt2),ls='--', c='#4444FF',lw=3, label='Model')
                
                
                ax2.plot(yeswav,yesflx, color ='#666666', zorder = 2)
                idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
                #eqws = np.append(eqws,VEW)ustedindices[i]]).argmin()
                #plt.title(linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
                ax2.set_xlabel('Wavelength (Angstroms)')
                ax2.set_ylabel('Flux')
                ax2.vlines(wav[adjustedindices[i]],0,1.1*Gbackground, zorder = 1, color = '#7B8fB9')
                ax2.vlines(wav[adjustedindices[i]]-(GEW/2),0,Gbackground, zorder = 1, color = '#4444FF')
                ax2.vlines(wav[adjustedindices[i]]+(GEW/2),0,Gbackground, zorder = 1, color = '#4444FF')
                ax2.set_title('Gaussian Model EW = ' + '{:.3f}'.format(1000*GEW) + 'm\u212b')
                if len(yeswav) > 0: ax2.hlines(Gbackground,yeswav[0], yeswav[-1], color = '#4444FF')
                ax2.set_ylim(0,1.1*maxflx)
                #plt.show(block = False)
                    
            except:
                counter = 1
                ax2.plot(yeswav,yesflx, color ='#666666', zorder = 2)
                idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
                ax2.set_title(linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
                ax2.set_xlabel('Wavelength (Angstroms)')
                ax2.set_ylabel('Flux')
                #plt.show(block = False)
                instring = ''
                if not autoDiscard:
                    print('Could not model with Gaussian distribution. Enter new parameters')
                    while instring != 'p' and len(instring) < 7:
                        instring = input('Separated by spaces, Input the following.\n(a, sigma, mu, background).\nOr enter "p" to skip this line: ')
                    if instring == "p":
                        GEW = -1
                        #eqws = np.append(eqws,EW)
                        break
                    guesses = [float(i) for i in instring.split()]
                else:
                    GEW=-1
                    counter = 0
        
        
        if not lotsofblendedlines:
            # takes the continuum to be the average of the points toward the edge of the area of interest
            Obackground = np.average( np.append( yesflx[int(len(yesflx)/2- 1.5*FWHMS[i]-5):int(len(yesflx)/2- 1.5*FWHMS[i])-0], yesflx[int(len(yesflx)/2 +1.5*FWHMS[i]+0):int(len(yesflx)/2 +1.5*FWHMS[i]+5)])) 
        else:
            # takes the continuum to be the average of all flux data > 0.9
            Obackground = np.average(yesflx[yesflx >= 0.9*maxflx])
        
        OEW = 0
        for j in range(len(yesflx)-1):
            OEW += (1 - (yesflx[j]/Obackground))*(yeswav[j+1] - yeswav[j]) 
        
        #plt.figure()    
        ax3.plot(yeswav,yesflx, color ='#666666', zorder = 2)
        ax3.vlines(wav[adjustedindices[i]],0,1.1, zorder = 1, color = '#7B8fB9')
        ax3.vlines(wav[adjustedindices[i]]-(OEW/2),0,Obackground, zorder = 1, color = '#041E42')
        ax3.vlines(wav[adjustedindices[i]]+(OEW/2),0,Obackground, zorder = 1, color = '#041E42')
        if len(yeswav) > 0: ax3.hlines(Obackground,yeswav[0], yeswav[-1], color = '#041E42')
        #ax3.set_xticks(fontsize = '8')
        #plt.ylim(max([flx[adjustedindices[i]] - 0.5*(1-flx[adjustedindices[i]]),0]),1.025*maxflx)
        try:
            plt.ylim(0,1.1*maxflx)
        except:
            pass
        idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
        ax3.set_title('Observed Data EW = ' + '{:.3f}'.format(1000*OEW) + 'm\u212b')
        ax3.set_xlabel('Wavelength (Angstroms)')
        ax3.set_ylabel('Flux')
        plt.show(block = False)
        
        
        EWs = np.sort([OEW,GEW,VEW])
        if np.std(EWs) <  .2*EWs[1]:
            eqws = np.append(eqws,np.average(EWs))
            print('\n' + linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) +  '\u212B')
            print('EW = ' + '{:.3f}'.format(1000*np.average(EWs)))
            plt.close('all')
        elif not autoDiscard:
            instring = ''
            while instring != 'G' and instring != 'V' and instring != 'O' and instring != 'D':
                
                print('\n' + linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) +  '\u212B')
                print('EW_G = ' +'{:.3f}'.format(1000*GEW) + ' m\u212B)')
                print('EW_V = ' +'{:.3f}'.format(1000*VEW) + ' m\u212B)')
                print('EW_O = ' +'{:.3f}'.format(1000*OEW) + ' m\u212B)')
                instring = input('Use EW from (G)aussian, (V)oigt, or (O)bserved measurement? Or, (D)iscard line?').capitalize()
            if instring == 'G':
                eqws = np.append(eqws,GEW)
            if instring == 'V':
                eqws = np.append(eqws,VEW)
            if instring == 'O':
                eqws = np.append(eqws,OEW)
            if instring == 'D':
                eqws = np.append(eqws,-1)
            plt.close('all')
        else:
            eqws = np.append(eqws,-1)
            print('Discarded ' + linelist['Element'][idk] + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) +  '\u212B')
            plt.close('all')
            
# Outputting the line list to a file
if MOOGformat:
    
    # Replace name with number because MOOG demands it
    for i in range(len(linelist['Element'])):
        for j in range(len(elnames)):
            if linelist['Element'][i] == elnames[j]:
                linelist['Element'][i] = str(elnums[j])
                
    for i in range(len(linelist['Ionization'])):
        if str(linelist['Ionization'][i]) == 'I' or str(linelist['Ionization'][i]) == '1':
            linelist['Ionization'][i] = '.0'
            
        elif str(linelist['Ionization'][i]) == 'II' or str(linelist['Ionization'][i]) == '2':
            linelist['Ionization'][i] = '.1'

        elif str(linelist['Ionization'][i]) == 'III' or str(linelist['Ionization'][i]) == '3':
            linelist['Ionization'][i] = '.2'
            
        elif str(linelist['Ionization'][i]) == 'IV' or str(linelist['Ionization'][i]) == '4':
            linelist['Ionization'][i] = '.3'            
    
    linelist = linelist[linelist['Wavelength'] <= max(wav)]
    linelist = linelist[linelist['Wavelength'] >= min(wav)]
    linelist = linelist[1000*eqws>=mineqw]
    eqws = eqws[1000*eqws>=mineqw]
    
    with open(outfile, 'w') as f:
        print('NLINES: %s, OBJECT: %s' % (str(len(eqws)), obj), file = f)
        for i in range(len(eqws)):
            print('{:10}'.format(linelist['Wavelength'][i]) + '{:>10}'.format(linelist['Element'][i]+linelist['Ionization'][i]) + '{:>10}'.format(str(linelist['ExcitationPot'][i])) + '{:>10}'.format(str(linelist['Loggf'][i])) + '{:>30}'.format('{:10f}'.format(eqws[i]*1000)), file = f)

    
else:
    
    linelist = linelist[linelist['Wavelength'] <= max(wav)]
    linelist = linelist[linelist['Wavelength'] >= min(wav)]
    linelist = linelist[1000*eqws>=mineqw]
    eqws = eqws[1000*eqws>=mineqw]
    
    with open(outfile, 'w') as f:
        print(obj + ' WAV,ELEM,ION,EP,LOGGF,EQW', file = f)
        for i in range(len(eqws)):
            print(str(linelist['Wavelength'][i]) + ',' + str(linelist['Element'][i]) + ',' + str(linelist['Ionization'][i]) + ',' + str(linelist['ExcitationPot'][i]) + ',' + str(linelist['Loggf'][i]) + ',' + str(eqws[i]*1000), file = f)

