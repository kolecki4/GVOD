import numpy as np
import astropy.io.fits as fits


#=============================================================================
# This program does most of the heavy-lifting of measuring absorption lines for you.
# It calculates the equivalent width of all the given lines in the spectrum supplied
# in units of milliAngstroms with graphs for verification purposes. 

# CHANGE THE PARAMETERS BELOW ACCORDINGLY
#=============================================================================

# FULL PATH OF SPECTRUM IN FITS FORMAT
#filename = "/media/jared/Acer/Users/jared/Desktop/ResearchFolder/spectra/BD+183423bspectrum.dxt.nor"
filename='1D_Spec_HR7672.txt'

# OPEN THE 1D SPECTRUM FILE HERE (MAKE SURE YOU IMPORT ANY AND ALL LIBRARIES NECESSARY)

#fitsfile= fits.open(filename')
wav =np.genfromtxt(filename, usecols=0)
flx =np.genfromtxt(filename, usecols=1)

# FULL PATH OF LINE LIST TEXT FILE (COLUMNS Element,Ionization,Wavelength,ExcitationPotential,Loggf)
linelistname = "Master Line List.txt"

# CHOOSE FROM THE FOLLOWING METHODS OF CALCULATING EQUIVALENT WIDTHS:
# 'obs': Use the observed data
# 'gau': Fit a Gaussian profile to each line
# 'voi': Fit a Voigt profile to each line
# 'ask': Ask for one of the above for lines where the above differ by more than 20%
method = 'ask'

# FOR ask MODE ONLY, AUTO-DISCARD ANY LINE WHERE THE ABOVE DIFFER BY MORE THAN 20%
autoDiscard = False

# RAISE THIS NUMBER IF YOUR CALCULATED EQUIVALENT WIDTHS ARE TOO LOW (RANGE 0:1 EXCLUSIVE)
# GENERALLY USE 0.5
FWHeight = 0.5

# INCREASE halfwindow TO SEARCH A WIDER RANGE FOR THE MIN FLUX VALUES
# (USEFUL IF THE SPECTRUM IS NOT RV OR REDSHIFT ADJUSTED, FOR EXAMPLE)
dlambda = np.average([wav[i+1] - wav[i] for i in range(len(wav)-1)])
halfwindow = int(.35/dlambda)

# MINIMUM EQW IN MILLIANGSTROMS TO BE OUTPUT TO THE LINE LIST FILE
# (MEASUREMENTS OF LINES WEAKER THAN THIS ARE CONSIDERED TO BE TOO NOISY FOR USE)
mineqw = 5

# DO YOU WANT YOUR LINE LIST TO BE FORMATTED FOR THE MOOG SPECTRAL ANALYSIS PROGRAM?
MOOGformat = True

# DIRECTORY, FILE, AND OBJECT NAME FOR SAVING THE LINE LIST
outfile = 'TestLineList.txt'
obj = 'HR 7672'

# PLAY WITH colnames SO THEY'RE IN THE CORRECT ORDER AND MAKE SURE YOU'RE GENNING FROM TXT CORRECTLY FOR YOUR FILE
# ALSO elnums MUST BE IN THE SAME ORDER AS elnames
colnames = ('Element', 'Ionization', 'Wavelength', 'ExcitationPot', 'Loggf')
linelist = np.genfromtxt(linelistname, dtype = None, encoding = None, names = colnames, skip_header = 16, usecols = [0,1,2,3,4])
elnums = np.array([2,3,6,7,8,11,12,13,14,16,20,22,23,24,26,28,39,40,56,63])
elnames = np.array(['He','Li','C','N','O','Na','Mg','Al','Si','S','Ca','Ti','V','Cr','Fe','Ni','Y','Zr','Ba','Eu'])
