#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: mtormo
"""

import argparse
import multiprocessing
import numpy as np
import pandas as pd
import os
import sys

from collections import OrderedDict as od
from collections import Counter

def InputPar():  #input parameters
	parser = argparse.ArgumentParser(
		description='Extract stats (c10/c30 or c1/c0) per gene in a list of samples (From Bedtools coverage) ')
	parser.add_argument('--list_files', metavar='FILE', required=True,
						help='File with a list of files from bedtools. Each file should have 6 fields (tab-delimited): chrom - start - end - gene.exon - position - coverage (Required)')
	parser.add_argument('--out', metavar='INT', type=int, choices=[1,2,3], default=3,
						help='Output with coverage analysis type: 1 = Table with genes/samples computing c10, c30 and mean; 2 = Table with genes/samples computing c1 and getting 0-coverage positions; 3 = Both [3]')
	parser.add_argument('--base_out', metavar='BASENAME', required=True,
						help='Output basename: <BASE_OUT>_cov1030.tab for out=1, and  <BASE_OUT>_cov0.tab for out=2 (Required)')
	return parser.parse_args()


def MergeC(x):
	l = []
	m = []
	n = []
	# handle NaN values
	if pd.notnull(x['c10']):
		l = str(x['c10'])
	else:
		l = "0"

	if pd.notnull(x['c30']):
		m = str(x['c30'])
	else:
		m = "0"

	if pd.notnull(x['mean']):
		n = str(x['mean'])
	else:
		n = "0"

	# print type(l),type(m),type(n)
	# print x['c10'],x['c30'],x['mean']
	return '/'.join([l,m,n])



def GetCov(filein):
	### get dataframe
	infile = pd.read_csv(filein, sep="\t", header=None)
	### get identifier
	iden = os.path.split(filein)[1].split(".")[0]

	### get gene name from exons' names
	infile["gene"] = infile[3].str.split(".").str[0]

	### get c10,30 and mean
	covs = pd.DataFrame()
	covs["ct"] = infile.groupby("gene")[5].count()
	covs["ct10"] = infile[infile[5] >= 10].groupby("gene")[5].count()
	covs["ct30"] = infile[infile[5] >= 30].groupby("gene")[5].count()
	covs["c10"] = covs.ct10*100 / covs.ct
	covs["c30"] = covs.ct30*100 / covs.ct
	covs["mean"] = infile.groupby("gene")[5].mean()

	### join columns (c10/c30/mean)
	covs = covs.round(2)
	covs[iden] = covs.apply(lambda row: MergeC(row), axis=1)

	### get c1
	covs["ct1"] = infile[infile[5] >= 1].groupby("gene")[5].count()
	covs["c1"] = covs.ct1*100 / covs.ct

	### get exons names with cov0
	covs0 = {}
	for k,v in infile[infile[5] == 0].groupby("gene")[3]:
		covs0[k] = Counter(v)


	### remove columns
	covs = covs.drop(["ct1","ct30","ct10","ct","c10","c30","mean"],axis=1)


	return [iden,covs,covs0]


def main():
	### take arguments
	param = InputPar()

	infiles = []
	with open(param.list_files) as list:
		for line in list:
			line = line.rstrip()
			infiles.append(line)

	pool = multiprocessing.Pool()
	results = pool.map(GetCov, infiles)

	covs_list = [] ### list of coverage per sample to concatenate
	covs0_list = [] ### list of c1/c0
	for each in results:
		### get all results
		iden,covs,covs0 = each

		# ### extract c1 and remove the column
		for row,cols in covs.iterrows():
			if cols[1] < 100:
				# print covs0[row]
				c0 = ",".join([x + ":" + str(covs0[row][x]) for x in covs0[row].keys()])

				covs0_list.append([iden,row,str(round(cols[1],2))+"/"+c0])


		covs = covs.drop(["c1"], axis=1)
		### append coverage
		covs_list.append(covs)

	### concatenate list of coverage
	covs_out = pd.concat(covs_list, axis=1)

	# print type(covs_list)
	if param.out != 2:
		covs_out.to_csv(param.base_out+"_cov1030.tab", sep='\t')

	if param.out != 1:
		with open(param.base_out+"_cov0.tab","w") as out2:
				out2.write("sample\tgene\tc1/c0"+"\n")
				for item in covs0_list:
					out2.write("\t".join(item)+"\n")



if __name__ == "__main__":
	main()
