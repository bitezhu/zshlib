import matplot
import numpy as np
import sys

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
'''

a,b,c,d,e,f,g=bamio.mappingstat(sys.argv[1])
print a,b,c,d,e,f,g
#matplot.densityplot([a1,a2,a3,a4,a5,a6,a7,a8],['s','e','f','g','h','a','b','c'])
