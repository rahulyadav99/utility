import numpy as np
import matplotlib.pyplot as pl
import h5py
from tqdm import tqdm

def readout_1c_ch(obsf,nx,ny, invf):
	'''
	function reads only chromospheric inverted data

	inputs: 
		obsf: observed Stokes maps
		nx: xsize 
		ny: ysize
		invf: inverted maps
	output: 
		inv_ch: inverted parameters in 3D array [nx, ny, parameters]
		synpro: synthetic profiles in 4D array [nx, ny, 4stks, wav]
		chi2: chisquare map [nx, ny]
	Example:
		inv_ch, synpro,chi2 = readout_1c_ch[obsf,nx,ny, invf]
	'''
	fobs=h5py.File(obsf,'r')
	tpix, nlam, stk = fobs['stokes'].shape
	fobs.close()

	synpro = np.zeros((tpix,4,nlam), dtype=np.float64)

	chlabel = ['tau','v','Bx','By','Bz','a','beta','deltav','ff']
	nch = len(chlabel)

	maps_ch = np.zeros((tpix,nch), dtype=np.float64)

	f1 = h5py.File(invf,'r')
	key = list(f1.keys())

	print('Input inverted file -->',invf)
	print('Input observed file -->',obsf)

	for j in tqdm(range(nch)):
		maps_ch[:,j] = np.squeeze(f1['ch1'][chlabel[j]][:,:,-1,:])


	syn = np.squeeze(f1['spec1']['stokes'][:,0,-1,:,:])

	synpro[:,:,:] = syn

	inv_ch = np.reshape(maps_ch,(nx,ny,nch))
	synpro = np.reshape(synpro,(nx,ny,4,nlam))
	ch2 = f1['spec1']['chi2'][:]
	ch2 = np.squeeze(ch2[:,-1,-1])
	chi2 = np.reshape(ch2,(nx,ny))
	f1.close()
	return inv_ch,synpro,chi2

def readout_1c(obsf,nx,ny, invf):
	'''
	function reads 1-component of photospheric and chromospheric data 
	inputs: 
		obsf: observed Stokes maps
		nx: xsize 
		ny: ysize
		invf: inverted maps
	output: 
		inv_ch: chromospheric inverted parameters in 3D array [nx, ny, parameters]
		inv_ph: photospheric inverted parameters in 4D array [nx, ny, logtau, parameters]
		logtau: logtau array
		phff: photospheric filling factor
		synpro: synthetic profiles in 4D array [nx, ny, 4stks, wav]
		chi2: chisquare map [nx, ny]
	Example:
		inv_ch, inv_ph, synpro, logtau, phff,chi2 = readout_1c_ch[obsf,nx,ny, invf]
	'''
	fobs=h5py.File(obsf,'r')
	tpix, nlam, stk = fobs['stokes'].shape
	fobs.close()

	synpro = np.zeros((tpix,4,nlam), dtype=np.float64)

	phlabel = ['Bx','By','Bz','T','v','vmic']
	chlabel = ['tau','v','Bx','By','Bz','a','beta','deltav','ff']
	nch = len(chlabel)
	nph = len(phlabel)

	maps_ch = np.zeros((tpix,nch), dtype=np.float64)
	maps_ph = np.zeros((tpix,73,nph), dtype=np.float64)
	phff = np.zeros((tpix,), dtype=np.float64)

	f1 = h5py.File(invf,'r')
	key = list(f1.keys())

	print('Input inverted file -->',invf)
	print('Input observed file -->',obsf)

	for j in tqdm(range(nch)):
		maps_ch[:,j] = np.squeeze(f1['ch1'][chlabel[j]][:,:,-1,:])

	for p in tqdm(range(nph)):
		maps_ph[:,:,p] = np.squeeze(f1[key[1]][phlabel[p]][:,:,-1,:])

	phff[:] = np.squeeze(f1[key[1]]['ff'][:,:,-1,:])

	syn = np.squeeze(f1['spec1']['stokes'][:,0,-1,:,:])

	synpro[:,:,:] = syn

	inv_ch = np.reshape(maps_ch,(nx,ny,nch))
	inv_ph = np.reshape(maps_ph,(nx,ny,73,nph))
	synpro = np.reshape(synpro,(nx,ny,4,nlam))
	logtau = np.squeeze(f1[key[1]]['log_tau'][:])
	ch2 = f1['spec1']['chi2'][:]
	ch2 = np.squeeze(ch2[:,-1,-1])
	chi2 = np.reshape(ch2,(nx,ny))
	f1.close()
	return inv_ch, inv_ph,synpro, logtau, phff,chi2

def readout_2c(obsf,nx,ny, invf):
	'''
	function reads 2-comp of chromospheric and 1-comp of photospheric data
	inputs: 
		obsf: observed Stokes maps
		nx: xsize 
		ny: ysize
		invf: inverted maps
	output: 
		inv_ch: chromospheric inverted parameters in 4D array [nx, ny, parameters,2]
		inv_ph: photospheric inverted parameters in 4D array [nx, ny, logtau, parameters]
		logtau: logtau array
		phff: photospheric filling factor
		synpro: synthetic profiles in 4D array [nx, ny, 4stks, wav]
		chi2: chisquare map [nx, ny]
	Example:
		inv_ch, inv_ph,synpro, logtau, phff,chi2 = readout_1c_ch[obsf,nx,ny, invf]
	'''
	fobs=h5py.File(obsf,'r')
	tpix, nlam, stk = fobs['stokes'].shape
	fobs.close()

	synpro = np.zeros((tpix,4,nlam), dtype=np.float64)

	phlabel = ['Bx','By','Bz','T','v','vmic']
	chlabel = ['tau','v','Bx','By','Bz','a','beta','deltav','ff']
	nch = len(chlabel)
	nph = len(phlabel)

	maps_ch12 = np.zeros((tpix,nch,2), dtype=np.float64)
	maps_ph = np.zeros((tpix,73,nph), dtype=np.float64)
	phff = np.zeros((tpix,), dtype=np.float64)

	f1 = h5py.File(invf,'r')
	print('Input inverted file -->',invf)
	print('Input observed file -->',obsf)

	for j in tqdm(range(nch)):
		maps_ch12[:,j,0] = np.squeeze(f1['ch1'][chlabel[j]][:,:,-1,:])
		maps_ch12[:,j,1] = np.squeeze(f1['ch2'][chlabel[j]][:,:,-1,:])

	for p in tqdm(range(nph)):
		maps_ph[:,:,p] = np.squeeze(f1['ph1'][phlabel[p]][:,:,-1,:])	

	phff[:] = np.squeeze(f1['ph1']['ff'][:,:,-1,:])

	syn = np.squeeze(f1['spec1']['stokes'][:,0,-1,:,:])

	synpro[:,:,:] = syn

	inv_ch = np.reshape(maps_ch12,(nx,ny,nch,2))
	inv_ph = np.reshape(maps_ph,(nx,ny,73,nph))
	synpro = np.reshape(synpro,(nx,ny,4,nlam))
	logtau = np.squeeze(f1['ph1']['log_tau'][:])
	ch2 = f1['spec1']['chi2'][:]
	ch2 = np.squeeze(ch2[:,-1,-1])
	chi2 = np.reshape(ch2,(nx,ny))
	f1.close()
	return inv_ch, inv_ph,synpro, logtau, phff, chi2
