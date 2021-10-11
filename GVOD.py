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

print(halfwindow)
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
    
    # Gets relevant wavelength and flux data
    yesflx = flx[adjustedindices[i]-(2*FWHMS[i]):adjustedindices[i]+(2*FWHMS[i])] 
    yeswav = wav[adjustedindices[i]-(2*FWHMS[i]):adjustedindices[i]+(2*FWHMS[i])]
    
    
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
        plt.title(str(str(linelist['Element'][idk])) + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B (EW = ' +'{:.3f}'.format(1000*EW) + ' m\u212B)')
        plt.xlabel('Wavelength (Angstroms)')
        plt.ylabel('Normalized Flux')
        plt.show()
    
    if method == 'gau':
        EW = 0
        OEW = 0
        try:
            popt2, pcov=curve_fit(gaussian,yeswav,yesflx,p0=[maxflx - minflx,wav[adjustedindices[i]],FWHMS[i]*(wav[adjustedindices[i]]-wav[adjustedindices[i]-1]),maxflx])
            background=popt2[3]
            for j in range(len(yeswav)-1):
                EW += (1 - (gaussian(yeswav,*popt2)[j]/background))*(yeswav[j+1] - yeswav[j]) 
            
            if not lotsofblendedlines:
            # takes the continuum to be the average of the points toward the edge of the area of interest
                Obackground = np.average( np.append( yesflx[int(len(yesflx)/2- 1.5*FWHMS[i]-5):int(len(yesflx)/2- 1.5*FWHMS[i])-0], yesflx[int(len(yesflx)/2 +1.5*FWHMS[i]+0):int(len(yesflx)/2 +1.5*FWHMS[i]+5)])) 
            else:
            # takes the continuum to be the average of all flux data > 0.9
                Obackground = np.average(yesflx[yesflx >= 0.9*maxflx])
            
            for j in range(len(yesflx)-1):
                OEW += (1 - (yesflx[j]/Obackground))*(yeswav[j+1] - yeswav[j]) 
                
                
            #plt.plot(yeswav, gaussian(yeswav,*popt2),ls='--', c='#4444FF',lw=3, label='Model')
        except:
            print('Could not model with Gaussian.')
            EW=0
            
        if 1:#np.abs(OEW-EW) < .5*EW:
            #plt.plot(yeswav,yesflx, color ='#666666', zorder = 2)    
            #plt.vlines(wav[adjustedindices[i]],0,1.1, zorder = 1, color = '#7B8fB9')
            #plt.vlines(wav[adjustedindices[i]]-(EW/2),0,background, zorder = 1, color = '#041E42')
            #plt.vlines(wav[adjustedindices[i]]+(EW/2),0,background, zorder = 1, color = '#041E42')
            #if len(yeswav) > 0: plt.hlines(background,yeswav[0], yeswav[-1], color = '#041E42')
            #plt.xticks(fontsize = '8')
            #plt.ylim(max([flx[adjustedindices[i]] - 0.5*(1-flx[adjustedindices[i]]),0]),1.025*maxflx)
            try:
                plt.ylim(0,1.1*maxflx)
            except:
                pass
            idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
            #plt.title(str(str(linelist['Element'][idk])) + ' ' + str(linelist['Ionization'][idk]) + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B (EW = ' +'{:.3f}'.format(1000*EW) + ' m\u212B)')
            #plt.xlabel('Wavelength (Angstroms)')
            #plt.ylabel('Normalized Flux')
            #plt.show()
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
            try:
                plt.ylim(0,1.1*maxflx)
            except:
                pass
            idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
            plt.title(str(linelist['Element'][idk]) + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B (EW = ' +'{:.3f}'.format(1000*EW) + ' m\u212B)')
            plt.xlabel('Wavelength (Angstroms)')
            plt.ylabel('Normalized Flux')
            plt.show()
        else:
            eqws = np.append(eqws,-1)
    
    
    if method == 'ask':
        instring = 'W'
        while instring in ['W','J','C','I']:
            if not lotsofblendedlines:
                # takes the continuum to be the average of the points toward the edge of the area of interest
                Obackground = np.average( np.append( yesflx[int(len(yesflx)/2- 1.5*FWHMS[i]-5):int(len(yesflx)/2- 1.5*FWHMS[i])-0], yesflx[int(len(yesflx)/2 +1.5*FWHMS[i]+0):int(len(yesflx)/2 +1.5*FWHMS[i]+5)])) 
            else:
                # takes the continuum to be the average of all flux data > 0.9
                Obackground = np.average(yesflx[yesflx >= 0.9*maxflx])
                
            if instring == 'C' or instring == 'I':
                Obackground = Cont

            
            fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (16,4))
            idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
            try:        
                fig.suptitle(str(linelist['Element'][idk]) + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
            except:
                fig.suptitle(str(linelist['Element'][idk]) + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
            counter = 1
            guesses = [wav[adjustedindices][i],maxflx,maxflx - minflx,1,.1,Obackground]  
            # Voigt EW      
            while counter == 1:
                try:
                    counter = 0
                    VEW = 0
                    
                    popt1, pcov=curve_fit(Voigt,yeswav,yesflx,p0=guesses) 
                    Vbackground=popt1[5]
                    for j in range(len(yeswav)-1):
                        VEW += (1 - (Voigt(yeswav,*popt1)[j]/Vbackground))*(yeswav[j+1] - yeswav[j])
                        if VEW > maxeqw/1000:
                            VEW = -1
                            break
                    
                    ax2.plot(yeswav,yesflx, color ='#666666', zorder = 2)
                    ax2.set_xlabel('Wavelength (Angstroms)')
                    ax2.set_ylabel('Flux')
                    ax2.plot(yeswav, Voigt(yeswav,*popt1),ls='--', c='#FF4444',lw=3, label='Model')
                    ax2.vlines(wav[adjustedindices[i]],0,1.1*Vbackground, zorder = 1, color = '#7B8fB9')
                    ax2.vlines(wav[adjustedindices[i]]-(VEW/2),0,Vbackground, zorder = 1, color = '#FF4444')
                    ax2.vlines(wav[adjustedindices[i]]+(VEW/2),0,Vbackground, zorder = 1, color = '#FF4444')
                    ax2.set_title('Voigt Model EW = ' + '{:.3f}'.format(1000*VEW) + 'm\u212b')
                    if len(yeswav) > 0: ax2.hlines(Vbackground,yeswav[0], yeswav[-1], color = '#FF4444')
                    if minflx < 0.9*maxflx: ax2.set_ylim(0,1.1*maxflx)
                    else: ax2.set_ylim(0.99*minflx,1.01*maxflx)
                        
                except:
                    counter = 1
                    ax2.plot(yeswav,yesflx, color ='#666666', zorder = 2)
                    idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
                    ax2.set_xlabel('Wavelength (Angstroms)')
                    ax2.set_ylabel('Flux')
                    try: 
                        ax2.set_title(str(linelist['Element'][idk]) + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
                    except:
                        ax2.set_title(str(linelist['Element'][idk]) + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
                    instring = ''
                    VEW = -1
                    counter=0

            # Gaussian EW    
            counter = 1
            guesses = [maxflx - minflx,wav[adjustedindices][i],FWHMS[i]*(wav[adjustedindices[i]]-wav[adjustedindices[i]-1]),Obackground]
            while counter == 1:
                try:
                    counter = 0
                    GEW = 0
                    
                    popt2, pcov=curve_fit(gaussian,yeswav,yesflx,p0=guesses)
                    Gbackground=popt2[3]
                    for j in range(len(yeswav)-1):
                        GEW += (1 - (gaussian(yeswav,*popt2)[j]/Gbackground))*(yeswav[j+1] - yeswav[j]) 
                        if GEW > maxeqw/1000:
                            GEW = -1
                            break
                    ax1.plot(yeswav, gaussian(yeswav,*popt2),ls='--', c='#4444FF',lw=3, label='Model')
                    ax1.plot(yeswav,yesflx, color ='#666666', zorder = 2)
                    idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
                    ax1.set_xlabel('Wavelength (Angstroms)')
                    ax1.set_ylabel('Flux')
                    ax1.vlines(wav[adjustedindices[i]],0,1.1*Gbackground, zorder = 1, color = '#7B8fB9')
                    ax1.vlines(wav[adjustedindices[i]]-(GEW/2),0,Gbackground, zorder = 1, color = '#4444FF')
                    ax1.vlines(wav[adjustedindices[i]]+(GEW/2),0,Gbackground, zorder = 1, color = '#4444FF')
                    ax1.set_title('Gaussian Model EW = ' + '{:.3f}'.format(1000*GEW) + 'm\u212b')
                    if len(yeswav) > 0: ax1.hlines(Gbackground,yeswav[0], yeswav[-1], color = '#4444FF')
                    if minflx < 0.9*maxflx: ax1.set_ylim(0,1.1*maxflx)
                    else: ax1.set_ylim(0.99*minflx,1.01*maxflx)
                        
                except:
                    counter = 1
                    ax1.plot(yeswav,yesflx, color ='#666666', zorder = 2)
                    idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
                    try: 
                        ax1.set_title(str(linelist['Element'][idk]) + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
                    except:
                        ax1.set_title(str(linelist['Element'][idk]) + ' at ' + str(linelist['Wavelength'][idk]) + ' \u212B')
                    ax1.set_xlabel('Wavelength (Angstroms)')
                    ax1.set_ylabel('Flux')
                    instring = ''
                    GEW=-1
                    counter = 0
            
            # Observed EW
            OEW = 0
            for j in range(len(yesflx)-1):
                OEW += (1 - (yesflx[j]/Obackground))*(yeswav[j+1] - yeswav[j])
                if OEW > maxeqw/1000:
                    OEW = -1
                    break
            
            ax3.plot(yeswav,yesflx, color ='#666666', zorder = 2)
            ax3.vlines(wav[adjustedindices[i]],0,1.1*Obackground, zorder = 1, color = '#7B8fB9')
            ax3.vlines(wav[adjustedindices[i]]-(OEW/2),0,Obackground, zorder = 1, color = '#041E42')
            ax3.vlines(wav[adjustedindices[i]]+(OEW/2),0,Obackground, zorder = 1, color = '#041E42')
            if len(yeswav) > 0: ax3.hlines(Obackground,yeswav[0], yeswav[-1], color = '#041E42')
            if minflx < 0.9*maxflx: plt.ylim(0,1.1*maxflx)
            else: plt.ylim(0.99*minflx,1.01*maxflx)

            idk = np.abs(linelist['Wavelength'] - wav[adjustedindices[i]]).argmin()
            ax3.set_title('Observed Data EW = ' + '{:.3f}'.format(1000*OEW) + 'm\u212b')
            ax3.set_xlabel('Wavelength (Angstroms)')
            ax3.set_ylabel('Flux')
            
            
            ax4x = wav[(wav[adjustedindices[i]]-4 < wav) & (wav < wav[adjustedindices[i]]+4)]
            ax4y = flx[(wav[adjustedindices[i]]-4 < wav) & (wav < wav[adjustedindices[i]]+4)]
            
            ax4.plot(ax4x,ax4y, c = '#666666')
            if len(yeswav) > 0: ax4.plot(np.append(np.append(yeswav[0],yeswav),yeswav[-1]),np.append(np.append(Obackground,yesflx), Obackground),color ='red', zorder = 2)
            try:
                ax4.set_ylim(0,1.05*max(flx[(wav[adjustedindices[i]]-4 < wav) & (wav < wav[adjustedindices[i]]+4)]))
            except:
                pass
            plt.show(block = False)
            
            
            EWs = np.sort([OEW,GEW,VEW])
            if autoAccept and np.std(EWs) <  .2*EWs[1]:
                instring = ''
                eqws = np.append(eqws,np.average(EWs))
                try:
                    print('\n' + str(linelist['Element'][idk]) + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) +  '\u212B')
                except:
                    print('\n' + str(linelist['Element'][idk]) + ' at ' + str(linelist['Wavelength'][idk]) +  '\u212B')
                print('EW = ' + '{:.3f}'.format(1000*np.average(EWs)))
                plt.close('all')
            elif not autoDiscard:
                instring = ''
                while instring not in ['G','V','O','D','W','M','K','L',';','\'','J','P','C','I']:     
                    try: 
                        print('\n' + str(linelist['Element'][idk]) + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) +  '\u212B')
                    except:
                        print('\n' + str(linelist['Element'][idk]) + ' at ' + str(linelist['Wavelength'][idk]) +  '\u212B')
                    try:
                        print('EW_G = ' +'{:.3f}'.format(1000*GEW) + ' m\u212B)' + ', sum(resids) = ' + str(np.sum(np.abs(gaussian(yeswav,*popt2)-yesflx))))
                    except:     
                        print('EW_G = ' +'{:.3f}'.format(1000*GEW) + ' m\u212B)')
                    try:                   
                        print('EW_V = ' +'{:.3f}'.format(1000*VEW) + ' m\u212B)' + ', sum(resids) = ' + str(np.sum(np.abs(Voigt(yeswav,*popt1)-yesflx))))
                    except:
                        print('EW_V = ' +'{:.3f}'.format(1000*VEW) + ' m\u212B)')                    
                    print('EW_O = ' +'{:.3f}'.format(1000*OEW) + ' m\u212B)')
                    instring = input('Use EW from (G,L)aussian, (V,;)oigt, or (O,\')bserved measurement?\nAdditionally, you may change the (W,J)indow size of the calculations,\n(M,P)anually define an EW, change the (C,I)ontinuum level, or (D,K)iscard line.\n').capitalize()
                if instring == 'G' or instring == 'L':
                    eqws = np.append(eqws,GEW)
                elif instring == 'V' or instring == ';':
                    eqws = np.append(eqws,VEW)
                elif instring == 'O' or instring == '\'':
                    eqws = np.append(eqws,OEW)
                elif instring == 'D' or instring == 'K':
                    eqws = np.append(eqws,-1)
                elif instring == 'M' or instring == 'P':
                    validMEW = False
                    while not validMEW:                    
                        try:
                            MEW = float(input('New EW in m\u212B: '))/1000
                            validMEW = True
                        except:
                            pass
                    eqws = np.append(eqws,MEW)
                    
                elif instring == 'C' or instring == 'I':
                    validCont = False
                    while not validCont:                    
                        try:
                            Cont = float(input('Old continuum level: ' + str(Obackground) + '\nNew continuum: '))
                            validCont = True
                        except:
                            pass  
                    
                elif instring == 'W' or instring == 'J':
                    validFWH = False
                    while not validFWH:                    
                        try:
                            FWHeight = float(input('Old window size: ' + str(FWHeight) + '\nNew window size (FWHeight, range 0:1 exclusive): '))
                            validFWH = True
                        except:
                            pass
                    # Gets relevant wavelength and flux data
                    FWHM = peak_widths(1-flx, rel_height = FWHeight, peaks = adjustedindices)[0]
                    FWHMS = np.array([int(round(i,0)) for  i in FWHM])
                    yesflx = flx[adjustedindices[i]-(2*FWHMS[i]):adjustedindices[i]+(2*FWHMS[i])] 
                    yeswav = wav[adjustedindices[i]-(2*FWHMS[i]):adjustedindices[i]+(2*FWHMS[i])]
                    if len(yesflx) > 0:
                       maxflx = max(yesflx)
                       minflx = min(yesflx)
                    else:
                       maxflx = max(flx)
                       minflx = min(flx)

                plt.close('all')
            else:
                instring = ''
                eqws = np.append(eqws,-1)
                try:
                    print('Discarded ' + str(linelist['Element'][idk]) + ' ' + linelist['Ionization'][idk] + ' at ' + str(linelist['Wavelength'][idk]) +  '\u212B')
                except:
                    print('Discarded ' + str(linelist['Element'][idk]) + ' at ' + str(linelist['Wavelength'][idk]) +  '\u212B')
                plt.close('all')
            
# Outputting the line list to a file
if MOOGformat:
    
    # Replace name with number because MOOG demands it
    try:
        for i in range(len(linelist['Ionization'])):
            if str(linelist['Ionization'][i]) == 'I' or str(linelist['Ionization'][i]) == '1':
                linelist['Ionization'][i] = '.0'
                
            elif str(linelist['Ionization'][i]) == 'II' or str(linelist['Ionization'][i]) == '2':
                linelist['Ionization'][i] = '.1'
    
            elif str(linelist['Ionization'][i]) == 'III' or str(linelist['Ionization'][i]) == '3':
                linelist['Ionization'][i] = '.2'
                
            elif str(linelist['Ionization'][i]) == 'IV' or str(linelist['Ionization'][i]) == '4':
                linelist['Ionization'][i] = '.3' 
      
        for i in range(len(linelist['Element'])):
            for j in range(len(elnames)):
                if linelist['Element'][i] == elnames[j]:
                    linelist['Element'][i] = str(elnums[j])
    except:
        pass                     
    
    linelist = linelist[linelist['Wavelength'] <= max(wav)]
    linelist = linelist[linelist['Wavelength'] >= min(wav)]
    linelist = linelist[1000*eqws>=mineqw]
    eqws = eqws[1000*eqws>=mineqw]
    eqws = eqws[1000*eqws<=maxeqw]
    
    with open(outfile, 'w') as f:
        print('NLINES: %s, OBJECT: %s' % (str(len(eqws)), obj), file = f)
        try:        
            for i in range(len(eqws)):
                print('{:10}'.format(linelist['Wavelength'][i]) + '{:>10}'.format(linelist['Element'][i]+linelist['Ionization'][i]) + '{:>10}'.format(str(linelist['ExcitationPot'][i])) + '{:>10}'.format(str(linelist['Loggf'][i])) + '{:>30}'.format('{:10f}'.format(eqws[i]*1000)), file = f)
        except:
            print("Input linelist from linemake or other similar MOOG-like format")
            for i in range(len(eqws)):
                print('{:10}'.format(linelist['Wavelength'][i]) + '{:>10}'.format(linelist['Element'][i]) + '{:>10}'.format(str(linelist['ExcitationPot'][i])) + '{:>10}'.format(str(linelist['Loggf'][i])) + '{:>30}'.format('{:10f}'.format(eqws[i]*1000)), file = f)
    
else:
    
    linelist = linelist[linelist['Wavelength'] <= max(wav)]
    linelist = linelist[linelist['Wavelength'] >= min(wav)]
    linelist = linelist[1000*eqws>=mineqw]
    eqws = eqws[1000*eqws>=mineqw]
    eqws = eqws[1000*eqws<=maxeqw]
    
    with open(outfile, 'w') as f:
        print(obj + ' WAV,ELEM,ION,EP,LOGGF,EQW', file = f)
        for i in range(len(eqws)):
            print(str(linelist['Wavelength'][i]) + ',' + str(linelist['Element'][i]) + ',' + str(linelist['Ionization'][i]) + ',' + str(linelist['ExcitationPot'][i]) + ',' + str(linelist['Loggf'][i]) + ',' + str(eqws[i]*1000), file = f)
