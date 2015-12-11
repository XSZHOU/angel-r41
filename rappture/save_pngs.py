# ===================================================
# ANGEL - An NEGF simulator aimed at LEDs
# auxiliary python script to let matplotlib 
# generate contour plots of simulation results
# ===================================================


import sys
import matplotlib
import numpy as np
#from scipy import interpolate
from scipy import linspace
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'


if len(sys.argv) < 3:
	# no arguments for the simulation to visualize are given
	# --> look for info in result.files
	f = open('result.files', 'r')
	s1 = f.readline()
	assert s1!="", 'invalid result.files'
	s2 = f.readline()
	assert s2!="", 'invalid result.files'
	f.close()
	
	s1.strip()
	s2.strip()
	s1 =  s1[:len(s1)-1] # get rid of '\n'
	s2 =  s2[:len(s2)-1]
	
	basedir=s1
	basename=s2
else:
	basedir  = sys.argv[1]
	basename = sys.argv[2]
	

#basedir='../results/nanohub/'
#basename = 'nanohub_Source0.000V'

print "Starting ANGEL results image generation..."
print "basedir="+basedir
print "basename="+basename

my_cmap = cm.jet
fc = 'white'
shr=0.8
asp = 0.7 # for imshow
xsize=4 # inches
ysize=3 # inches
colbarext='neither' # arrow-like finish of colorbar or not
my_dpi = 200


def plot_it(x, y, f):
	plt.figure(facecolor=fc, figsize=(xsize,ysize))
	
	# C = plt.contourf(x,y,f,100,antialiased=1, cmap=my_cmap)
	
	# newfunc = interpolate.interp2d(x,y,f,kind='linear')
	assert(len(f)==len(y))
	assert(len(f[0])==len(x))
	
	#nx = 100
	#ny = 100
	nx = xsize*my_dpi / 4
	ny = ysize*my_dpi / 4
	
	xnew = linspace(min(x),max(x),nx)
	ynew = linspace(min(y),max(y),ny)
	# fnew = newfunc(xnew,ynew)
	fnew = [];
	
	# do interpolation onto rectilinear grid myself - interp2d is terribly slow
	xidx = range(nx);
	for ii in range(nx):
		xinter = xnew[ii];	
		test = 0;
		while test<len(x)-1 and x[test+1]<xinter:
			test = test+1;
		xidx[ii] = test;
		
	yidx = range(ny);
	for jj in range(ny):
		yinter = ynew[jj];	
		test = 0;
		while test<len(y)-1 and y[test+1]<yinter:
			test = test+1;
		yidx[jj] = test;
	
	for jj in range(ny):
		tmp = range(nx);
		fnew.append(tmp) # by reference --> can't take tmp definition outside for-loop!!!!!
				
		for ii in range(nx):
			
			xxx = xidx[ii];
			yyy = yidx[jj];
			f00 = f[yidx[jj]][xidx[ii]];
			
			if xxx<len(x)-1 and yyy<len(y)-1:
			
				fx = (xnew[ii] - x[xidx[ii]]) / (x[xidx[ii]+1] - x[xidx[ii]]); # fx=0 --> xnew[ii] = x[xidx[ii]]
				fy = (ynew[jj] - y[yidx[jj]]) / (y[yidx[jj]+1] - y[yidx[jj]]);
			
				f01 = f[yidx[jj]  ][xidx[ii]+1];
				f10 = f[yidx[jj]+1][xidx[ii]  ];
				f11 = f[yidx[jj]+1][xidx[ii]+1];
			
				#finter = f00;
				#finter = fy*f00 + (1-fy)*f01;
				finter = (1-fy) * ((1-fx)*f00 + fx*f01) + fy * ((1-fx)*f10 + fx*f11);
			else:
				finter = f00;
					
			fnew[jj][ii] = finter;
		
	#C = plt.pcolor(xnew,ynew,fnew, cmap=my_cmap)
	#C = plt.pcolor(x,y,f, cmap=my_cmap)
	
	#C = plt.imshow(f, interpolation='bilinear', origin='lower',
	C = plt.imshow(fnew, interpolation='bilinear', origin='lower',
          cmap=my_cmap, extent=(min(x), max(x), min(y), max(y)), 
		  aspect=asp*(max(x)-min(x))/(max(y)-min(y)))
	CBI = plt.colorbar(C, orientation='vertical', shrink=shr, extend=colbarext)

	return C




# -----------------------
# n(x,E)
# -----------------------
print "n(x,E)..."
nE_data = np.loadtxt(basedir+basename+'_nE', dtype='float', comments='%')
pE_data = np.loadtxt(basedir+basename+'_pE', dtype='float', comments='%')
# 1st col = x, 1st row = E
x1 = nE_data[ 0, 1:] # python indices start with 0
E1 = nE_data[1:,  0]

nE_data = nE_data[1:, 1:] # take away x,E
pE_data = pE_data[1:, 1:]
npE = nE_data + pE_data
C = plot_it(x1,E1,npE)

# -----------------------
# LDOS(x,E)
# -----------------------
print "LDOS(x,E)..."
LDOS_data = np.loadtxt(basedir+basename+'_LDOS', dtype='float', comments='%')
LDOS_data = LDOS_data[1:, 1:]
C = plot_it(x1,E1,LDOS_data)


# -----------------------
# Je(x,E)
# -----------------------
print "J(x,E)..."
JEe_data = np.loadtxt(basedir+basename+'_JEe', dtype='float', comments='%')
JEh_data = np.loadtxt(basedir+basename+'_JEh', dtype='float', comments='%')
x2 = JEe_data[ 0, 1:] # python indices start with 0
E2 = JEe_data[1:,  0]

JEe_data = JEe_data[1:, 1:]
JEh_data = JEh_data[1:, 1:]
JE = JEe_data + JEh_data
C = plot_it(x2,E2,JE)


# -----------------------
# format figures
# -----------------------
print "Formatting figures..."
#plt.figure(1); plt.title('n(x,E,k=0) [$nm^{-3} eV^{-1}$]')
#plt.figure(2); plt.title('LDOS(x,E,k=0) [$nm^{-3} eV^{-1}$]')
#plt.figure(3); plt.title('J(x,E,k=0) [$ec^{-1} ps^{-1} nm^{-2} eV^{-1}$]')
plt.figure(1); plt.title('n(x,E,k=0)')
plt.figure(2); plt.title('LDOS(x,E,k=0)')
plt.figure(3); plt.title('J(x,E,k=0)')
for ii in range(3):
	plt.figure(ii+1)
	plt.xlabel('x[nm]', fontsize=14, color='k')
	plt.ylabel('E[eV]', fontsize=14, color='k')
	plt.xlim(min(x1), max(x1));
	plt.ylim(min(E1), max(E1));
	
	# left  = 0.125  # the left side of the subplots of the figure
	# right = 0.9    # the right side of the subplots of the figure
	# bottom = 0.1   # the bottom of the subplots of the figure
	# top = 0.9      # the top of the subplots of the figure
	# wspace = 0.2   # the amount of width reserved for blank space between subplots
	# hspace = 0.2   # the amount of height reserved for white space between subplots
	plt.subplots_adjust(bottom=0.1, top = 0.95, left=0.1, right=0.75)



# -----------------------
# save to file
# .png's will be in same directory as the caller of the script is
# -----------------------
print "Saving to file "+basename+"_*.png..."
plt.figure(1); 
plt.savefig(basename+'_npE.png', dpi=my_dpi, facecolor='w', edgecolor='w', format='png', transparent=False)
plt.figure(2); 
plt.savefig(basename+'_LDOS.png', dpi=my_dpi, facecolor='w', edgecolor='w', format='png', transparent=False)
plt.figure(3); 
plt.savefig(basename+'_JE.png', dpi=my_dpi, facecolor='w', edgecolor='w', format='png', transparent=False)
		
		
# -----------------------
# show plots
# script stops until plots have been closed
# -----------------------
#print "Displaying plots..."; plt.show()
