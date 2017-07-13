##this is a ported version of the GNU Radio Post Process PDPs Doppler MATLAB code originally written by Prof. Chris Anderson and CDR Thomas W. Tedesso

##Author: Cristian Mendoza

##MAKE DEFAULT FONT SIZE/WEIGHT 14/Bold

import csv
import numpy as np
import matplotlib.pyplot as plt
#from scipy.io.numpyio import fread

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
	numpnseq = 1e4 #blocks to process
	segsize = pnlength*numpnseq #segment size to process

	count = 0
	
	#time to process.
	with open(datafile,'rb') as dataf:
		z = np.fromfile(dataf)
		z = [float(x) for x in z]
		z = np.reshape(z,(2,segsize))
		data = z[0,:]+z[1,:]*1j #complex numbaz
		[r, c] = data.shape
		data = np.reshape(data,c,r)
			

	



