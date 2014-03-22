# Copyright 2011 Dustin Lang and David W. Hogg.	 All rights reserved.
if __name__ == '__main__':
	import matplotlib
	matplotlib.use('Agg')

import os
import pyfits
import numpy as np
from astrometry.util.miscutils import *

import pylab as plt


# make normalized meshgrids
def hogg_meshgrid(W, H):
	mX, mY = np.meshgrid(np.arange(float(W)), np.arange(float(H)))
	# normalize to ~ +- 1
	mX -= W/2
	mY -= H/2
	# scale by same amount (chip is square anyway)
	mX /= float(W/2)
	mY /= float(W/2)
	return mX,mY

## HACK -- duplicated code -- produce and weight the components
# as produced for one reference image, but apply to a different image.
def synth_another(img, amps, shiftmax=1, maxorder=1):
	W, H = img.shape
	mX, mY = hogg_meshgrid(W, H)
	Isynth = np.zeros_like(img)
	i = 0
	for mname,mult in getmults(mX, mY, maxorder):
		Isynth += amps[i] * np.ones_like(img) * mult
		i += 1
		IM = img * mult
		shiftrange = range(-shiftmax, shiftmax+1)
		for dx in shiftrange: # x shift
			for dy in shiftrange: # y shift
				sh = np.zeros_like(IM)
				sx1,sx2 = get_overlapping_region(dx, dx+W-1, 0, W-1)
				sy1,sy2 = get_overlapping_region(dy, dy+H-1, 0, H-1)
				sh[sy2, sx2] = IM[sy1, sx1]
				Isynth += amps[i] * sh
				i += 1
	assert(i == len(amps))
	return Isynth
		

def getmults(mX, mY, maxorder):
	mults = [('Img', np.ones_like(mX)),
			 ('X', mX),
			 ('Y', mY)]
	if maxorder > 1:
		for order in range(2, maxorder+1):
			for ox in range(0, order+1):
				oy = order - ox
				mults.append(('X^%i Y^%i' % (ox,oy), mX**ox * mY**oy))
	return mults

# data: image to be synthesized
# invvar: inverse variances for data
# reference: image to be used to construct linear components
# shiftmax: maximum number of pixels permitted in shifts
def synthesize(data, invvar, reference, shiftmax=1, maxorder=1):
	W, H = reference.shape
	mX, mY = hogg_meshgrid(W, H)
	assert(data.shape == invvar.shape)
	assert(data.shape == reference.shape)
	# Zero outer edge of invvar.
	invvar[:shiftmax,:] = 0
	invvar[:,:shiftmax] = 0
	invvar[-shiftmax:,:] = 0
	invvar[:,-shiftmax:] = 0

	reference -= np.median(reference)

	#plt.clf()
	#plt.imshow(invvar, interpolation='nearest', origin='lower')
	#plt.gray()
	#plt.axis([-10, 510, -10, 510])
	#plt.savefig('invvarX.png')

	mults = getmults(mX, mY, maxorder)

	# component images
	Ishifts = []
	Iabout = []
	for mname,mult in mults:
		Ishifts.append(np.ones_like(reference) * mult) # piston;  X,Y gradients
		Iabout.append((mname, 0, 0, 0, '%s x 1' % (mname)))
		# shifted components
		IM = reference * mult
		shiftrange = range(-shiftmax, shiftmax+1)
		for dx in shiftrange: # x shift
			for dy in shiftrange: # y shift
				sh = np.zeros_like(IM)
				sx1,sx2 = get_overlapping_region(dx, dx+W-1, 0, W-1)
				sy1,sy2 = get_overlapping_region(dy, dy+H-1, 0, H-1)
				sh[sy2, sx2] = IM[sy1, sx1]
				Iabout.append((mname, 1, dx, dy,'%s (%+i, %+i)' % (mname, dx, dy)))
				Ishifts.append(sh)
	w = np.sqrt(invvar.ravel() / invvar.max())
	A = np.zeros((W * H, len(Ishifts)))
	for i,sh in enumerate(Ishifts):
		A[:, i] = sh.ravel() * w
	X = data.ravel() * w
	amplitudes, nil, nil, nil = np.linalg.lstsq(A, X)
	synth = np.sum( np.dstack(Ishifts) * amplitudes, axis=2)
	return synth, amplitudes, Ishifts, Iabout

def synthesize_fn(data, invvar, ref, shiftmax=1):
	return synthesize(*[pyfits.open(x)[0].data for x in [data,invvar,ref]],
					  shiftmax=shiftmax)

def percentiles(data, pcts):
	dr = data.ravel()
	I = np.argsort(dr)
	return [ dr[I[int(len(dr) * p)]] for p in pcts ]
	

if __name__ == '__main__':

	from astrometry.util.sdss_das import sdss_das_get
	from astrometry.util.sdss_cutout import cutout

	ra,dec = 53.202125, -0.365361

	# 7101-g3-0751
	if True:
		for i,(run,band,camcol,field) in enumerate([(6955, 'g', 3, 809),
													(6458, 'g', 3, 410),
													(7101, 'g', 3, 751),
													]):
			fpC = sdss_das_get('fpC', None, run, camcol, field, band=band)
			fpM = sdss_das_get('fpM', None, run, camcol, field, band=band)
			psField = sdss_das_get('psField', None, run, camcol, field, band=band)
			#tsField = sdss_das_get('tsField', None, run, camcol, field, band=band)
			#outfn = sdss_filename(filetype, run, camcol, field, band) + suffix

			print 'fpC is', fpC
			cutout(fpC, ra, dec, 200, 'data%i.fits' % i, fpM, psField, 'invvar%i.fits' % i,
				   band)


			

	import sys
	args = sys.argv[1:]

	reffn = 'data0.fits'
	ref	 = pyfits.open(reffn)[0].data

	plt.figure(1)
	plt.clf()

	ND = 2
	ploti = 1
	im = dict(interpolation='nearest', origin='lower')
	for i in range(1, ND+1):
		datafn = 'data%i.fits' % i
		invvarfn = 'invvar%i.fits' % i
		data = pyfits.open(datafn)[0].data
		invvar = pyfits.open(invvarfn)[0].data
		print 'synthesizing from:', datafn, invvarfn

		synthfn = 'synth%i.fits' % i
		if not os.path.exists(synthfn):
			synth,amps,components = synthesize(data, invvar, ref, shiftmax=4)
			pyfits.writeto(synthfn, synth, clobber=True)
		else:
			synth = pyfits.open(synthfn)[0].data

		mn,mx = percentiles(data, [0.05, 0.999])
		ima = dict(interpolation='nearest', origin='lower',
				   vmin=mn, vmax=mx)

		NC = 4
		plt.figure(1)
		plt.subplot(ND, NC, ploti)
		ploti += 1
		plt.imshow(data, **ima)
		plt.xticks([])
		plt.yticks([])
		plt.gray()

		plt.subplot(ND, NC, ploti)
		ploti += 1
		plt.imshow(synth, **ima)
		plt.xticks([])
		plt.yticks([])

		diff = data - synth
		dd = (mx - mn)/2.

		plt.subplot(ND, NC, ploti)
		ploti += 1
		plt.imshow(diff, vmin=-dd, vmax=dd, **im)
		plt.xticks([])
		plt.yticks([])

		plt.subplot(ND, NC, ploti)
		ploti += 1
		chi = (data - synth) * np.sqrt(invvar)
		plt.imshow(chi, vmin=-10, vmax=10, **im)
		plt.xticks([])
		plt.yticks([])
		#plt.colorbar()



		plt.figure(2)
		plt.clf()
		plt.imshow(data, **ima)
		plt.colorbar()
		plt.savefig('data%i.png' % i)
		plt.clf()
		plt.imshow(synth, **ima)
		plt.colorbar()
		plt.savefig('synth%i.png' % i)
		plt.clf()
		plt.imshow(diff, vmin=-dd, vmax=dd, **im)
		plt.colorbar()
		plt.savefig('diff%i.png' % i)
		plt.clf()
		plt.imshow(chi, vmin=-10, vmax=10, **im)
		plt.colorbar()
		plt.savefig('chi%i.png' % i)


	plt.figure(1)
	plt.savefig('all.png')


	if False:
		for i,(a,comp) in enumerate(zip(amps,components)):
			plt.clf()
			plt.imshow(comp, **im)
			plt.savefig('comp-%03i.png' % i)

		dr = ref.ravel()
		I = np.argsort(dr)
		rmn = dr[I[int(len(dr)*0.05)]]
		rmx = dr[I[int(len(dr)*0.999)]]
		plt.clf()
		plt.imshow(ref, interpolation='nearest', origin='lower',
				   vmin=rmn, vmax=rmx)
		plt.savefig('ref.png')
	
