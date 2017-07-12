import matplot
import numpy as np
import sys
import utils


#print utils.configMap(sys.argv[1])

print utils.commaFormat(123456789)
#matplot.boxplot([[3,5,2,6],[4,5,6,7],[6,7,8,9]],['a','b','c'],)
a1=np.random.random(10)*10
a2=np.random.random(10)*10
a3=np.random.random(10)*10
a4=np.random.random(10)*10
a5=np.random.random(10)*10
a6=np.random.random(10)*10
a7=np.random.random(10)*10
a8=np.random.random(10)*10

#for a,b in matplot.styleNum(10):
a,b=matplot.styleNum(10)
#print a 
#print b

import bamio

'''
idx=bamio.Tabix(sys.argv[1])

for item in idx.fetch('1',3100,5000):
        print item
        idx.close()
        a,b,c,d,e,f,g=bamio.mappingstat(sys.argv[1])
        print a,b,c,d,e,f,g
        #matplot.densityplot([a1,a2,a3,a4,a5,a6,a7,a8],['s','e','f','g','h','a','b','c'])
        '''

arr=[['1',42,52],['11',45,78],['2',25,100],['1',23,78],['1',56,89]]

print utils.sortArr(arr,0,1)
from format import *

f=sys.argv[1]
recs = fasta_itr(f)
print dir(recs)
print type(recs)

for rec in recs:
    print rec.id
    print rec.seq
