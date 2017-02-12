#!/usr/bin/env python
import numpy as np
import sys
# takes in input a  methylation counts file
# and a giml output file

methfn=sys.argv[1]
gimlfn=sys.argv[2]

# for the moment assume there is only 
# 1 chromosome ,1 lambda

meth=np.loadtxt(fname=methfn,dtype="S10,int,int,int")
print meth


