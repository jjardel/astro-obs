
# simple moving window skysubtraction
import pyfits
import numpy as np
import sys

def kappaSigma(data, kappa=2.5, maxiter=10):
	# simple kappa sigma klipping
	ii    = ~ np.isnan(data)
	count = 0
	while True:
		m = np.median(data[ii])
		s = np.std(data[ii])
		iinew = abs(data-m) <= kappa*s
 		if all(ii == iinew) or count > maxiter:
			break
		ii *= iinew
		count += 1	

	return ii


def movingWin(i, ii, width=20):
	# implements a moving window
	# i : index, i>0 and i<len(ii)
	# ii : boolean array 
	# width : window width
	# retrun boolean array, only those entries are
	#  True which lie in the window
	# At the edges the window will be smaller.
	lb = max(0,i-width/2) # lower bound
	ub = min(len(ii)-1,i+width/2) # lower bound
	jj = ii.copy()
	jj[:lb] *= False
	jj[ub:] *= False
	
	return jj

# load data
hdulist = pyfits.open("mFepesvw003927.fits")
data = hdulist[0].data

# we do some simple flatfielding here
# this needs improvement, see flat.fits to see what I mean.
# load flat data
hdulist2 = pyfits.open("Feeve_mastertrace.fits")
flatdata = hdulist2[0].data
for i in range(flatdata.shape[1]):
	ii = kappaSigma(flatdata[:,i])
	m = np.mean(flatdata[ii,i])
	flatdata[:,i] /= m
hdulist[0].data = flatdata
hdulist.writeto("flat.fits",clobber=True)
data /= flatdata
hdulist[0].data = data
hdulist.writeto("flattened.fits",clobber=True)


# create collapse image along spectal direction
cim = np.median(data, axis=1)
# reject bright or dark columns 
ii = kappaSigma(cim, kappa=2.5, maxiter=10)

# ols version, sky spectrum, by collapsing all skyfibers
# create collapsed sky spectrum
#csky = np.mean(data[ii,:],axis=0)
#skyim = np.ones(data.shape) * csky

# new version, use moving window
newdata = data.copy()
for i in range(len(ii)):
	jj =  movingWin(i, ii)
	sky = np.median(data[jj,:],axis=0)
	#sky = np.mean(data[jj,:],axis=0)
	
	newdata[i] -= sky

# save file
hdulist[0].data = newdata
hdulist.writeto("outfile.fits",clobber=True)


