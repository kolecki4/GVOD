# EWC
Planning on running abfind on MOOG? You're going to need to measure equivalent widths first!

This is a program capable of measuring these absorption lines by a number of different ways: just using the observed spectrum, fitting a Gaussian model, or fitting a Voigt model. Or you can make it do all three and decide for yourself which fits best.

### Using the EWC
Firstly, you'll need to open Settings.py and modify it accordingly. In order, the parameters are as follows:

##### filename
The 1D spectrum file. A file for testing purposes has been provided.

##### wav, flx
The wavelength and flux data arrays. Depending on the format of your spectrum file, it may be necessary to add additional lines of code before the final assignment of wav and flx.

##### linelistname
The spectral line information file. The list provided was originally taken from [Fulbright (2000)](https://iopscience.iop.org/article/10.1086/301548) with data appended from [Reddy et al. (2003)](https://academic.oup.com/mnras/article/367/4/1329/1079817) and the [NIST database](https://physics.nist.gov/PhysRefData/ASD/lines_form.html).

##### method
The method by which the equivalent width is calculated
Valid arguments are:
- 'obs': Use the observed data
- 'gau': Fit a Gaussian profile to each line
- 'voi': Fit a Voigt profile to each line
- 'ask': Ask for one of the above for lines where the above differ by more than 20%

##### autoDiscard
Only used while method == 'ask', if set to true, rather than ask for user input, the program will simply discard measurements where the above differ by more than 20%

##### FWHeight
Affects the width of the EW calculation window. Default is 0.5, range is (0,1).

##### mineqw
The minimum width which lines must reach to be added to the final list. Generally should be higher for noisier data. Default is 5 milliAngstroms

##### MOOGformat
If true, the final line list will be immediately readable by MOOG, else it will be output as a csv

##### outfile, object
The outfile variable is set to the desired destination of final line list, while object is an optional label which gets written to the final line list.



After verifying that your settings are correct, run ```python EWC.py``` and let it do its thing!
