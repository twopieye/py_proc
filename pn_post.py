##this is a ported version of the GNU Radio Post Process PDPs Doppler MATLAB code originally written by Prof. Chris Anderson and CDR Thomas W. Tedesso

##Author: Cristian Mendoza

##MAKE DEFAULT FONT SIZE/WEIGHT 14/Bold

import csv
import numpy as np
import matplotlib.pyplot as plt
#from scipy.io.numpyio import fread
import scipy.signal as sig
import matplotlib.axes as ax
def post_proc(fs):
	#fs = raw_input('enter the sampling freqency in MHz >>')
	fs = fs*1e6
	Ts = 1/fs

	filename = 'pntest.txt'# raw_input('enter PN sequence filename >> ')

	with open(filename) as pnfile:
		readtxt = csv.reader(pnfile)
		j = []
		for row in readtxt:
			for i in row:
				i = int(i)
				j = np.append(j,row[i])
	pnseq = [int(x) for x in j]#convert strings from csv to integers
	pnseq = [2*x-1 for x in pnseq] #antipodal sequence
	pnlength = len(pnseq) #length of PN seq in Chips
	#print pnlength
	pntime = [x*Ts for x in range(pnlength)]

	#let's plot this
	plt.figure(1)
	plt.plot([x*1e6 for x in pntime],pnseq)
	plt.xlabel('Time (us)')
	plt.ylabel('Amplitude (V)')
	plt.title('PN Sequence')
	plt.show()

	#load data
	datafile = '2015.09.28.13.01.31.dat'
	#data = open(datafile,'rb') #read binary data

	#Read/process data in blocks
	numpnseq = int(1e4) #blocks to process
	segsize = int(pnlength*numpnseq) #segment size to process

	count = 0
	
	#time to process.
	with open(datafile,'rb') as dataf:
		z = np.fromfile(dataf,dtype=np.float)
		zd = np.reshape(z,[2,-1])
		data = zd[0,:]+zd[1,:]*1j #complex numbaz
		#data = np.reshape(data,(len(data),1)) ####This part probably needs editing. gotta use segsize in here at some point.

		#gonna downsample this signal:
		ds_data = sig.decimate(data,2)

		#split complex numbaz for inphase/quadrature to gen. Mag/Phase PDPs
		i_data = np.real(ds_data)
		q_data = np.imag(ds_data)

		#run cross-correlation, generate time vector, and plot results
		#time vector starts at an offset of the PN length
		pdp_i = np.correlate(i_data,pnseq)
		pdp_q = np.correlate(q_data,pnseq)
		L = len(pdp_i)
		#throw away 0's
		pdp_i = pdp_i[int(np.floor((L)/2)):L]
		pdp_q = pdp_q[int(np.floor((L)/2)):L]
		#grab mag, phase
		pdp_mag = np.sqrt(pdp_i**2+pdp_q**2)
		pdp_mag_db = 20*np.log10(pdp_mag/max(pdp_mag))#dBFS, decibels relative to full scale
		pdp_phase = [x*180/np.pi for x in np.arctan2(pdp_q,pdp_i)]#convert radians to degrees
		#they always say "you make time":
		pdp_time = [x*Ts for x in range(L/2)]
		#plot mag
		plt.figure(2)
		plt.subplot(211)
		plt.plot([x*1e6 for x in pdp_time],pdp_mag_db,'k-')
		plt.xlabel('Time (us)')
		plt.ylabel('Relative Power (dBFS)')
		#plot phase
		plt.subplot(212)
		plt.plot([x*1e6 for x in pdp_time],pdp_phase,'r-')
		plt.xlabel('Time (us)')
		plt.ylabel('Phase (degrees)')
		plt.show()
		
		
		
		

		
			

	



