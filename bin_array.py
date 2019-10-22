import numpy as np
import matplotlib.pyplot as pl
from scipy import signal
from tqdm import tqdm
import h5py

def binavg1D(xdata,ydata,nbin):
	#input: 1D arrray, number of bins (nbin)
	#output: 1D array of nbinned data
	
	nx = len(ydata)
	binv = int((nx)/nbin)

	if (nbin/nx > 1):
		print('Error: Bin size > the array size!')
	else:
		ybin = np.zeros(binv,dtype = np.float64)
		xbin = np.zeros(binv,dtype = np.float64)
		for i in range(binv):
			ybin[i] = np.mean(ydata[i*nbin:i*nbin+nbin])
			xbin[i] = np.mean(xdata[i*nbin:i*nbin+nbin])
		#print('binned array',binned)
		return xbin, ybin


def binavg2D(data,nbin):
	#input: 2D array, number of bins (nbin, same for x & y direction) 
	#output: binned 2D array

	nx,ny = data.shape
	npx = int(nx/nbin)
	npy = int(ny/nbin)
	rebin = np.zeros((npx, npy),dtype = np.float64) 
	if (nbin/nx > 1):
		print('Error: Bin size > the array size!')
	else:
		for i in range(npx):
			for j in range(npy):
				rebin[i,j] = np.mean(data[i*nbin:i*nbin+nbin,j*nbin:j*nbin+nbin])

		#print('binned array',binned)
	return rebin

'''if __name__ == "__main__":
	#
	# routine for spectral binning
	#
	wav = np.loadtxt('wavelength_full.txt')
	ff = h5py.File('fullscan_stokes_f.h5','r')

	#f2 = ff.get('stokes_cor')
	cor = ff.get('stokes_fil')
	obs = ff.get('stokes_obs')

	lmb,nx,ny,stk = obs.shape
	data2bin = np.zeros((505,nx,ny,stk), dtype=np.float64)
	for i in tqdm(range(nx)):
		for j in range(ny):
			for k in range(4):
				y = np.squeeze(cor[:,i,j,k])
				x,yb = binavg1D(wav,y,2)
				data2bin[:,i,j,k] = yb
				
	hf = h5py.File('/scratch/rahul/grisinv/fullscan_stokes_f2bin.h5', 'w')
	hf.create_dataset('stokes_fil',data=data2bin, dtype=np.float64)
	hf.create_dataset('wavelength',data=x, dtype=np.float64)
	hf.create_dataset('stokes_obs',data=obs, dtype=np.float64)
	hf.close()
		
	gauss= signal.gaussian(100, std=7)+np.random.rand(100,)*0.4
	xx = np.linspace(-10,10,len(gauss))
	#x,y = binavg1D(xx,gauss,4)
	#pl.plot(xx,gauss)
	#pl.plot(x,y)
	#pl.show()
	
	nx = 400
	ny = 400
	data = np.zeros((nx,ny),dtype = np.float64) 
	data = data + np.random.randn(nx,ny)
	#for i in range(nx):
	#	for j in range(ny):
	#		data[i,j]=(i-200)**2+(j-200)**2
	aa = binavg2D(data,4)
	fig, ax = pl.subplots(nrows=1,ncols=2,figsize=(6,5))
	ax = ax.flatten()
	im=ax[0].imshow(data,vmax=5,vmin=-5,cmap = 'gray')
	#fig.colorbar(im, ax=ax[0])
	
	im2 = ax[1].imshow(aa,vmax=5,vmin=-5,cmap = 'gray')
	fig.colorbar(im, ax=ax)
	pl.show()'''
	
