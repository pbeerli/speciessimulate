#!/usr/bin/env python
import os
import sys
import random
loci = int(sys.argv[1])
file = sys.argv[2]
seed = int(sys.argv[3])
random.seed(seed)
os.system("python speciessim.py "+str(loci)+" "+file+" "+str(seed))
for i in range(1,loci):
    seed = random.randint(1,2**63)
    os.system("python speciessim.py -1 "+file+" "+str(seed))
