import numpy as np
import matplotlib.pyplot as pl
import h5py
import time
import scipy.io as io
from scipy.io import readsav as readsav
import astropy.io.fits as fits
import lptools as lp
import sparsetools as sp
import glob
import imtools as mt
'''
program to see observed, fitted Stokes profiles 
and inverted maps
'''
def xyplot(ix,iy):
	ww = []
	stk = []
	for i in range(4):
		ww.append(o[i].wav)
		stk.append(o[i].dat)

	nx = int(np.ceil(ix))
	ny = int(np.ceil(iy))
	pl.close(2)
	fig, ax = pl.subplots(nrows=5, ncols=3, figsize=(10,11))
	#ax = ax.flatten()
	ytit = ['I/I$_c$', 'Q/I$_c$ [%]', 'U/I$_c$ [%]', 'V/I$_c$ [%]']
	w0 = [3933.64,8542.1,6173.34]

	for nc in range(len(ss)):

		ll = [0,2,3]
		for k in range(len(ll)):
			z=np.where(stk[ll[k]][0,0,0,:,0] !=0)

			ax[0,k].plot(ww[ll[k]][z]-w0[k], np.squeeze(stk[ll[k]][nt,nx,ny,z,0]),'k.')
			ax[0,k].plot(ss[nc][ll[k]].wav-w0[k], ss[nc][ll[k]].dat[nt,nx,ny,:,0],label='cyc_'+str(nc+1))
			ax[0,k].set_xlabel('$\lambda$ [$\AA$]')
			ax[0,k].set_ylabel(ytit[0])

			if k!=0:
				for i in range(1,4):
					ax[k,i-1].plot(ww[ll[k]][z]-w0[k], np.squeeze(stk[ll[k]][nt,nx,ny,z,i])*100,'k.')
					ax[k,i-1].plot(ss[nc][ll[k]].wav-w0[k], ss[nc][ll[k]].dat[nt,nx,ny,:,i]*100,label='cyc_'+str(nc+1))
					ax[k,i-1].set_ylabel(ytit[i])
					ax[k,i-1].set_xlabel('$\lambda$ [$\AA$]')
	ax = ax.flatten()
	ytit = [ 'T [kK]', r'V$_\mathrm{l.o.s}$ [km s$^{-1}$]',r'V$_\mathrm{turb}$ [km s$^{-1}$]', r'B$_\parallel$ [kG]', r'$|$B$_\bot|$ [kG]', r'$\chi$ [deg]' ]
	for i in range(len(mm)):
		m = mm[i]

		scl = [1.e-3, 1.e-5, 1.e-5, 1.e-3, 1.e-3, 180./3.1415926]
		var = [m.temp, m.vlos, m.vturb, m.Bln, m.Bho, m.azi]
		tau = m.ltau[0,0,0]
		ma = [10, 5, 5, 2, 3, 180]
		mi = [3.,-5, 0, -2, 0, 0]
		
		npar = len(var)
		for iax in range(npar):
			ax[9+iax].plot(tau, var[iax][0,nx,ny,:]*scl[iax],label='cyc_'+str(i+1))
			ax[9+iax].axvspan(tau[45], tau[45], linestyle = '--',ymin=-110.1, ymax=10e3, alpha=1, color='black')
			ax[9+iax].axvspan(tau[25], tau[25], linestyle ='-.',ymin=-110.1, ymax=10e3, alpha=1, color='black')
			ax[9+iax].axvspan(tau[32], tau[32],linestyle =':',ymin=-110.1, ymax=10e3, alpha=1, color='black')

			ax[9+iax].set_ylim(mi[iax], ma[iax])
			ax[9+iax].set_xlim(-6.5, 0.0)
			ax[9+iax].set_ylabel(ytit[iax])
			ax[9+iax].set_xlabel(r'$\log \tau_{500}$')
	pl.tight_layout()
	pl.legend()
	pl.show(block=False)

def chimap(synf, obsf,nt=0):
	ss = sp.profile(synf)
	ob = sp.profile(obsf)

	z=np.where(ob.dat[nt,0,0,:,0] !=0)
	syn = np.squeeze(ss.dat[nt,:,:,z,0])
	obs = np.squeeze(ob.dat[nt,:,:,z,0])
	nw,px,py = obs.shape
	chisqr = np.zeros((px,py))

	for nx in range(px):
		for ny in range(py):
			ch = (syn[:,nx,ny]-obs[:,nx,ny])**2
			chisqr[nx,ny]=np.sum(ch)

	return chisqr

def onclick(event):
	#pl.close(1)
	global ix, iy
	ix, iy = event.xdata, event.ydata
	print("x-pos:",np.ceil(ix)," y-pos:", np.ceil(iy))
	#print(ix,iy)
	xyplot(iy,ix)


if __name__ == "__main__":
	'''
	map1 & map2: 4d array [nx, ny, ns, nw]
	w1 and w2 = wavelength array
	'''
	global ss,mm, o, nt
	
	indx = '185_135_15_15sr_025_049'
	
	#--read input files (observed and synthetic profiles, and model atmosphere)
	sf = np.sort(glob.glob('/proj/solarphysics/users/x_rahya/outputs/synthe*_'+indx+'_cyc2*.nc'))
	mf = np.sort(glob.glob('/proj/solarphysics/users/x_rahya/outputs/atmos*_'+indx+'_cyc2*.nc'))
	#of = '/home/x_rahya/stic/data/observed_40_114_370_444.nc'
	of = '/proj/solarphysics/users/x_rahya/data/observed_'+indx+'.nc'
	##--time frame index---	
	nt = 0
	
	#--list synthetic profiles per cycle
	ss=[]
	for i in range(len(sf)):
		ss.append(sp.profile(sf[i]).splitRegions())

	#--list model atmosphere per cycle
	mm=[]
	for i in range(len(mf)):
		mm.append(sp.model(mf[i]))
		
	o = sp.profile(of).splitRegions()

	print(sf)
	print(mf)
	print(of)

	#call chisquare routine...
	chisq = chimap(sf[-1],of,nt=nt)

	#--read model atmosphere at different optical depth
	temp = mm[-1].temp
	tau = mm[-1].ltau[0,0,0,:]
	
	pht = temp[nt,:,:,45]
	mch = temp[nt,:,:,32]
	uch = temp[nt,:,:,25]
	tpara = [pht, mch, uch]

	temp = mm[-1].vlos
	pht = temp[nt,:,:,45]
	mch = temp[nt,:,:,32]
	uch = temp[nt,:,:,25]
	vpara = [pht, mch, uch]

	temp = mm[-1].Bln
	pht = temp[nt,:,:,54]
	mch = temp[nt,:,:,32]
	uch = temp[nt,:,:,25]
	bpara = [pht, mch, uch]

	temp = mm[-1].Bho
	pht = temp[nt,:,:,54]
	mch = temp[nt,:,:,32]
	uch = temp[nt,:,:,25]
	bhor = [pht, mch, uch]

	pp = [tpara, vpara, bpara, bhor]
	#labls = ['Photo','middle chromo', 'upper chromo']

	labls = [r'T [kK] $\log \tau_{500}$ = '+str(np.round(tau[45],2)),r'T [kK] $\log \tau_{500}$ = '+str(np.round(tau[32],2)), r'T [kK] $\log \tau_{500}$ = '+str(np.round(tau[25],2))]

	labls = [r'$\log \tau_{500}$ = '+str(np.round(tau[54],1)),r'$\log \tau_{500}$ = '+str(np.round(tau[32],1)), r'$\log \tau_{500}$ = '+str(np.round(tau[25],1))]
	paralbl = ['T [kK]',  r'V$_\mathrm{l.o.s}$ [km s$^{-1}$]', r'B$_\parallel$ [kG]',r'$|$B$_\bot|$ [kG]']

	scl = [1.e-3, 1.e-5, 1.e-3, 1.e-3]
	#var = [m.temp, m.vlos, m.vturb, m.Bln, m.Bho, m.azi]
	#tau = m.ltau[0,0,0]
	ma = [8, 5, 2, 2]
	mi = [3.,-5,-2, -2]
	cm = ['hot','bwr', 'terrain', 'terrain']
	fig, ax = pl.subplots(figsize=(7.3,7), nrows=4,ncols=3,sharex=True, sharey=True)
	#ax = ax.flatten()
	for p in range(len(pp)):
		for i in range(3):
			if (p <3 or i < 2):
				im = ax[p][i].imshow(pp[p][i]*scl[p], origin='lower', cmap=cm[p],vmin=mi[p],vmax=ma[p])
				if i ==2:fig.colorbar(im, ax=ax[p][2],label=paralbl[p])
				ax[p][i].text(2,2,labls[i],fontsize=9,color='k')
			else:
				im = ax[p][i].imshow(mt.histo_opt(chisq),origin = 'lower',vmax=1)
				fig.colorbar(im, ax=ax[p][i],label='chi2')
	pl.tight_layout()	

	cid = fig.canvas.mpl_connect('button_press_event', onclick)

	pl.show()
