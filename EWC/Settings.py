import numpy as np
import astropy.io.fits as fits


#=============================================================================
# CHANGE THE PARAMETERS BELOW ACCORDINGLY
#=============================================================================

# OPEN THE 1D SPECTRUM FILE HERE (MAKE SURE YOU IMPORT ANY AND ALL LIBRARIES NECESSARY)
filename = "1D_Spec_HR7672.txt"
#fitsfile = fits.open(filename)
wav = np.genfromtxt(filename, usecols=0)
flx = np.genfromtxt(filename, usecols=1)

# DIRECTORY, FILE, AND OBJECT NAME FOR SAVING THE LINE LIST
outfile = 'TestLineList.txt'
obj = 'HR 7672'

# FULL PATH OF LINE LIST TEXT FILE (COLUMNS Element,Ionization,Wavelength,ExcitationPotential,Loggf)
linelistname = "CNOS.txt"

# CHOOSE FROM THE FOLLOWING METHODS OF CALCULATING EQUIVALENT WIDTHS:
# 'obs': Use the observed data
# 'gau': Fit a Gaussian profile to each line
# 'voi': Fit a Voigt profile to each line
# 'ask': Ask for one of the above for lines where the above differ by more than 20%
method = 'ask'

# FOR ask MODE ONLY, AUTO-DISCARD ANY LINE WHERE THE ABOVE DIFFER BY MORE THAN 20%
autoDiscard = False
# FOR ask MODE ONLY, AUTO-ACCEPT ANY LINE WHERE THE ABOVE DIFFER BY LESS THAN 20%
autoAccept = False

# RAISE THIS NUMBER IF YOUR CALCULATED EQUIVALENT WIDTHS ARE TOO LOW (RANGE 0:1 EXCLUSIVE)
# GENERALLY USE 0.5
FWHeight = 0.4

# MINIMUM/MAXIMUM EQW IN MILLIANGSTROMS TO BE OUTPUT TO THE LINE LIST FILE
mineqw = 2
maxeqw = 150

# DO YOU WANT YOUR OUTPUT LINE LIST TO BE FORMATTED FOR THE MOOG SPECTRAL ANALYSIS PROGRAM?
MOOGformat = True

# PLAY WITH colnames SO THEY'RE IN THE CORRECT ORDER AND MAKE SURE YOU'RE GENNING FROM TXT CORRECTLY FOR YOUR FILE
# ALSO elnums MUST BE IN THE SAME ORDER AS elnames
colnames = ('Element', 'Ionization', 'Wavelength', 'ExcitationPot', 'Loggf')
# colnames = ('Wavelength', 'Element','ExcitationPot', 'Loggf')
linelist = np.genfromtxt(linelistname, dtype = None, encoding = None, names = colnames, skip_header = 0, usecols = range(len(colnames)), filling_values = '0', invalid_raise = False)
elnums = np.array([2,3,6,7,8,11,12,13,14,16,20,22,23,24,26,28,39,40,56,63])
elnames = np.array(['He','Li','C','N','O','Na','Mg','Al','Si','S','Ca','Ti','V','Cr','Fe','Ni','Y','Zr','Ba','Eu'])

# INCREASE halfwindow TO SEARCH A WIDER RANGE FOR THE MIN FLUX VALUES
dlambda = np.average([wav[i+1] - wav[i] for i in range(len(wav)-1)])
halfwindow = int(.15/dlambda)
