
import sys
import os
import random
import copy
import numpy as np


def readfile(path):
	dsnp = {}
	dpath = {'0':'control/', '1':'disease/'}
	for linfo in lid:
		sid = linfo[0]
		scondition = linfo[1]
		sinfo = '\t'.join(linfo)
		#print linfo
		with open(path+dpath[scondition]+sid+'_variants.vcf','r') as f:
			for l in f:
				if l.startswith('#'):
					continue
				r = l.strip().split()
				npos = int(r[1])
				if l.startswith('chr1\t') and npos < 1500000:
					
					if npos in dsnp.keys():
						if linfo not in dsnp[npos]:
							dsnp[npos].append(linfo)
					else:
						dsnp[npos] = [linfo]


				else:
					break

		
	print 'write to file'
	fout = open(path+'snpMat.txt','w')
	#nrow: num of indviduals, ncol: num of SNPs
	lsnp = sorted(dsnp.keys())
	ltitle = ''
	for k in lsnp:
		ltitle += 'chr1_'+str(k)+' '
	fout.write(ltitle + '\n')
	for linfo in lid:
		lnew = ''
		#print 'row per person', linfo
		for k in lsnp:
			if linfo in dsnp[k]:
				lnew += '1 '
			else:
				lnew += '0 '
		fout.write(lnew + '\n')

	fout.close()




if __name__ == '__main__':

	lid = []
	with open('covariates.csv', 'r') as f:
		for l in f:
			if l.startswith('human'):
				continue
			r = l.strip().split(',')
			lid.append(r[:2])


	path = 'vcf/'
	readfile(path)

