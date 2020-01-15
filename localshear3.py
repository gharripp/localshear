#!/usr/bin/env python

import sys

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import _cntr

#sys.path.append('../../../misc/')
#import myplotlib as my

from matplotlib import gridspec
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d

sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
import dd as dd

sys.path.append('/afs/ipp-garching.mpg.de/home/g/gharr/python3')

import map_eq_gharr as me

norm = np.linalg.linalg.norm


shot = 35929 #28389 # 32232

#trange=np.arange(5.58,5.6,0.01)  #USE THIS IF YOU WANT TIMETRACES
trange=np.arange(5.4,5.41,0.01)

diag='IDE' #"EQH" # IDE
exp='ABOCK' #AUGD #RRF
ed=1

order=r'test34483_gharr/'


def runavg(x,y,sum=100):

	mod=y.size%sum
	ysum=0
	ysum2=0
	xsum=0
	xsum2=0
	for j in np.arange(sum):
		ysum+=y[mod+j::sum]
		ysum2+=y[mod+j::sum]**2
		xsum+=x[mod+j::sum]
		xsum2+=x[mod+j::sum]**2

	std=np.sqrt(ysum**2-ysum2)/(sum-1)
	ysum/=sum
	ysum2/=sum
	xsum/=sum
	xsum2/=sum

	#xsum=x[sum-1::sum]

	return xsum, ysum


def curl(F, dx=1, dy=1):
	dFdx, dFdy, dump = np.gradient(F, dx, dy, 1)
	x =   dFdy[:,:,2]
	y = - dFdx[:,:,2]
	z =   dFdx[:,:,1] - dFdy[:,:,0]
	return np.array([x.T, y.T, z.T]).T

def curl_paddy(F, dx=1, dy=1):
	dFdx, dump, dFdy = np.gradient(F, dx, 1, dy)
	x = - dFdy[:,:,1]
	y =   dFdy[:,:,0] - dFdx[:,:,2]
	z =   dFdx[:,:,2]

	x =   dFdy[:,:,2]
	y = - dFdx[:,:,2]
	z =   dFdx[:,:,1] - dFdy[:,:,0]
	return np.array([x.T, y.T, z.T]).T


def curl(F, dx=1, dy=1, dz=1): # in poloidal plane
	Fx = F[:,:,0]
	Fy = F[:,:,1]
	Fz = F[:,:,2]
	dFxdy, dFxdx = np.gradient(Fx, dy, dx)
	dFydy, dFydx = np.gradient(Fy, dy, dx)
	dFzdy, dFzdx = np.gradient(Fz, dy, dx)
	dFxdz = dFydz = dFzdz = 0*dFxdy
	#i = j = 5
	#print 'dFxdy', dFxdy[i,j]
	#print 'dFxdx', dFxdx[i,j]
	#print 'dFydy', dFydy[i,j]
	#print 'dFydx', dFydx[i,j]
	#print 'dFzdy', dFzdy[i,j]
	#print 'dFzdx', dFzdx[i,j]
	return np.array([dFzdy.T - dFydz.T, dFxdz.T - dFzdx.T, dFydx.T - dFxdy.T]).T
	'''dFdx, dFdy, dFdz = np.gradient(F, dx, dy, dz)
	if no_dz: dFdz *= 0
	i = j = 5
	print 'dFdx', dFdx[i, j]
	print 'dFdy', dFdy[i, j]
	print 'dFdz', dFdz[i, j]
	x = dFdy[:,:,2] - dFdz[:,:,1]
	y = dFdz[:,:,0] - dFdx[:,:,2]
	z = dFdx[:,:,1] - dFdy[:,:,0]
	return np.array([x.FT, y.T, z.T]).T'''


#eqh.open(shot, exp, diag)
#eqh.open(shot, 'MICDU', 'EQI') #works
#eqh.open(shot, 'AUGD', 'IDE')
eq = me.equ_map()
#eq.open(shot, 'EQI', 'MICDU') #works
#eq.open(shot, 'IDE', 'AUGD')
eq.open(shot, exp, diag, ed)

eq.read_pfm()
eq.read_scalars()
eq.read_ssq()
ti=0
shotfile = True #False #True


R = eq.Rmesh #res['Ri']
z = eq.Zmesh

Bl = np.array(eq.Bmesh(trange)).T

# USE THIS IF YOU WANT S(THETA) AND S_MIDPLANE

s98 = np.zeros((trange.size,2*R.size+2*z.size))
s99 = np.zeros((trange.size,2*R.size+2*z.size))
s985 = np.zeros((trange.size,2*R.size+2*z.size))
s995 = np.zeros((trange.size,2*R.size+2*z.size))
theta98 = np.zeros((trange.size,2*R.size+2*z.size))
theta99 = np.zeros((trange.size,2*R.size+2*z.size))
theta985 = np.zeros((trange.size,2*R.size+2*z.size))
theta995 = np.zeros((trange.size,2*R.size+2*z.size))
smid985 = np.zeros(trange.size)
smid995 = np.zeros(trange.size)
smid99 = np.zeros(trange.size)
smid98 = np.zeros(trange.size)

#for t in np.arange(6.49,6.5,1e-2):
#for t in trange:

#THIS IS WHERE THE MAGIC HAPPENS

for i,t in enumerate(trange): #
	if True: #'pfm' not in globals():

		if shotfile:
			tidx = np.abs((eq.t_eq-t)).argmin()
			pfm = eq.pfm[:,:,tidx].T #res['pfm']

		else:
			R = np.linspace(-1, 1, 50)
			z = np.linspace(-1.5, 1.5, 50)
			pfm = np.zeros((len(z), len(R)))
			for i, zz in enumerate(z):
				pfm[i,:] = 1-(R**2+(2.*zz/3.)**2)/2.
			B = np.zeros(list(pfm.shape)+[3])
			B[:,:,2] = - np.sqrt(1-0.35*pfm)

		dR = np.average(R[1:]-R[:-1])
		dz = np.average(z[1:]-z[:-1])

		dpdz, dpdR = np.gradient(pfm, dz, dR)

		#plt.contour(R, z, pfm, 50)
		#plt.quiver(R, z, dpdR, dpdz)
		#plt.gca().set_aspect('equal')

		gradPsi = np.zeros(list(pfm.shape)+[3])
		gradPsi[:,:,0] = dpdR
		gradPsi[:,:,1] = dpdz

		if shotfile:
			B = Bl[:,:,i,:]
			#for i, zz in enumerate(z):
			#	for j, RR in enumerate(R):
			#		Bl = eq.Bmesh(t)
			#		B[i, j, 0] = Bl[0]
			#		B[i, j, 1] = Bl[1]
			#		B[i, j, 2] = Bl[2]

		psi_z, psi_R = np.gradient(pfm, dz, dR)
		psi_zz, psi_zR = np.gradient(psi_z, dz, dR)
		psi_Rz, psi_RR = np.gradient(psi_R, dz, dR)

		if not shotfile:
			B[:,:,0] =  -1./(2*np.pi*R)*psi_z
			B[:,:,1] =   1./(2*np.pi*R)*psi_R

		normGradPsi = np.linalg.linalg.norm(gradPsi, axis=2)
		normB = np.linalg.linalg.norm(B, axis=2)
		normGradPsi = np.array([normGradPsi.T]*3).T # 3D again...
		normB = np.array([normB.T]*3).T # 3D again...

		eperp = np.cross(gradPsi/normGradPsi, B/normB)				#COMPUTE PERPENDICULAR UNIT VECTOR
		kappa = np.cross(-B/normB,curl(-B/normB,dR,dz))
		ce = curl(eperp, dR, dz)
		sl =  np.sum(eperp*ce, 2)									#COMPUTE LOCAL SHEAR FOLLOWING MCCARTHEY
		#kn = np.sum(kappa*gradPsi/normGradPsi,2)					#COMPUTE NORMAL CURVATURE
		#kg = np.sum(kappa*np.cross(B/normB,gradPsi/normGradPsi),2)  #COMPUTE GEODAESIC CURVATURE

	#plotting
	#my.PGF(height=20, width=30)
	plt.figure(1)
	plt.clf()
	plt.gca().set_aspect('equal')
	plt.grid('off')
	c = plt.contourf(R, z, sl, levels=[-1,0,1])
	cs = plt.contour(R, z, sl, levels=[0], colors='white', linewidths=2)
	plt.contour(R, z, sl, levels=np.arange(-1,1.01,0.1), colors='white', linestyles=['solid'], linewidths=0.25)

	#plt.clabel(cs, inline=1, fmt=r'$s_\mathrm{local}=%2.1f$', fontsize=8, inline_spacing=10)
	#plt.text(2.12, -0.35, r'$s_\mathrm{local} = 0$', color='white', rotation=16)

	plt.xlabel(r'$\mathrm{R[m]}$')
	plt.ylabel(r'$\mathrm{z[m]}$')

	#plt.colorbar(c)

	if True: #'points' not in globals():
		angles = np.linspace(0, 2*np.pi, 360)
		rhos = np.linspace(0, 1.0, 100)

		points = np.zeros((len(rhos), len(angles), 2))
		#Bs = np.zeros((len(rhos), len(angles), 4))


		#for i in range(len(angles)):
		#	res = eqh.rhopol_to_Rz(t, rhos, angles[i])
		#	points[:, i, :] = np.array([res['R'], res['z']]).T
			#B = eqh.get_B(t, points[:, i, 0], points[:, i, 1])
			#Bs[:, i, :] = np.array([B[k] for k in ['Br', 'Bz', 'Bt', 'Bpol']]).T

		psiAx = eq.psi0
		psiSep = eq.psix
		zmag = eq.ssq['Zmag'][ti]
		Rmag = eq.ssq['Rmag'][ti]
		zidx=np.abs(z-zmag).argmin()
		Ridx=np.abs(R-Rmag).argmin()

		#R1 = np.average(points[:, :, 0], axis=1)
		#z1 = np.average(points[:, :, 1], axis=1)
		#ravg = np.average(np.sqrt((points[:, :, 0] - R1[0])**2 + (points[:, :, 1] - z1[0])**2), axis=1)
		#q = eqh.rhopol_to_q(t, rhos)['q']
		#r2q = interp1d(ravg, q, bounds_error=False, fill_value=np.nan)
		#def dqdr(x, h=0.0001):
		#	return (r2q(x+h)-r2q(x-h))/(2*h)

		#s = ravg/q * dqdr(ravg)

		rhomat = np.sqrt((pfm - psiAx[tidx])/(psiSep[tidx]-psiAx[tidx]))
		if 'rhomid' not in globals():
			rhomid = np.zeros((trange.size,R.size))
			smid = np.zeros((trange.size,R.size))
		tmp = interp2d(R,z,sl)
		Rnew,znew = np.meshgrid(R,z)
		c = _cntr.Cntr(Rnew, znew, rhomat)

		#THIS IS TO CALCUALTE S_L(THETA) FOR DISTINCT RHOS /NOT NEEDED FOR 2D PLOT

		trace = c.trace(0.98)[1]
		s98[ti,:trace.shape[0]] = np.array([tmp(trace[i,0], trace[i,1]) for i in range(len(trace))]).flatten()
		theta98[ti,:trace.shape[0]] = np.array([np.arctan2(trace[i,1] - zmag, trace[i,0] - Rmag) for i in range(len(trace))])/np.pi*180

		trace = c.trace(0.985)[1]
		s985[ti,:trace.shape[0]] = np.array([tmp(trace[i,0], trace[i,1]) for i in range(len(trace))]).flatten()
		theta985[ti,:trace.shape[0]] = np.array([np.arctan2(trace[i,1] - zmag, trace[i,0] - Rmag) for i in range(len(trace))])/np.pi*180

		trace = c.trace(0.99)[1]
		s99[ti,:trace.shape[0]] = np.array([tmp(trace[i,0], trace[i,1]) for i in range(len(trace))]).flatten()
		theta99[ti,:trace.shape[0]] = np.array([np.arctan2(trace[i,1] - zmag, trace[i,0] - Rmag) for i in range(len(trace))])/np.pi*180

		trace = c.trace(0.995)[1]
		s995[ti,:trace.shape[0]] = np.array([tmp(trace[i,0], trace[i,1]) for i in range(len(trace))]).flatten()
		theta995[ti,:trace.shape[0]] = np.array([np.arctan2(trace[i,1] - zmag, trace[i,0] - Rmag) for i in range(len(trace))])/np.pi*180
		rhomid[ti,Ridx:] = rhomat[zidx,Ridx:]
		smid[ti,Ridx:] = sl[zidx,Ridx:]
		smid995[ti]=np.interp(0.995,rhomid[ti],smid[ti])
		smid99[ti]=np.interp(0.99,rhomid[ti],smid[ti])
		smid98[ti]=np.interp(0.98,rhomid[ti],smid[ti])
		smid985[ti]=np.interp(0.985,rhomid[ti],smid[ti])
		#rho2s = interp1d(rhos, s, bounds_error=False, fill_value=np.nan)

	plt.contour(R, z, rhomat, levels=[1], colors='black', linewidths=2)
	#cs = plt.contour(R, z, rho2s(rhomat), colors='black', linewidths=1, levels=[0,0.5,1.5])
	#plt.clabel(cs, inline=1, fmt=r'$s=%2.1f$', fontsize=10, inline_spacing=10)

	plt.ylim(-1.2, 1.)
	plt.xlim(0.9,2.3)
	plt.title(r'$%i \,\,\, @ \,\,\, %3.2f~~\mathrm{s}$'%(shot, t))
	plt.tight_layout()

	print(t)
	ti=ti+1
	plt.savefig(order+r'%i_%3.2f.png'%(shot, t),dpi=300)

plt.show(False)

np.save(order+r's98',s98)
np.save(order+r's985',s985)
np.save(order+r's99',s99)
np.save(order+r's995',s995)
np.save(order+r'theta995',theta995)
np.save(order+r'rhomid',rhomid)
np.save(order+r'smid',smid)
np.save(order+r'smid995',smid995)
np.save(order+r'smid99',smid99)
np.save(order+r'smid98',smid98)
np.save(order+r'smid985',smid985)


try:
	__IPYTHON__
except:
	from IPython import embed
	embed()

#sl = np.sum((-eperp*ce)[:,:,:], axis=2)
