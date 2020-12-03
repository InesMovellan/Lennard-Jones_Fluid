import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc                         
from matplotlib.ticker import FormatStrFormatter  
from scipy.interpolate import interp1d            

font = {'family' : 'serif',   
        'weight' : 'normal',  
        'size' : '16'         
       }                      
plt.rc('font', **font)        

inc = []
n = []
g = []
f = open('g_out','r')
line = f.readline()
while line:
    l = line.split()
    inc.append(float(l[0]))
    n.append(float(l[1]))
    g.append(float(l[4]))
    line = f.readline()
f.close()

# Representation of the energy against the distortion coordinate
fig = plt.figure()
p = fig.add_subplot(111)
#p.plot(d,E,'ro', markersize=4)
#p.plot(cycle,pot,'ro')
#p.plot(inc,n,'ro', markersize=6)
p.plot(inc,g,'ro', markersize=5)
p.plot(inc,g,'r', linewidth=2.0)
axes = plt.gca()
#axes.set_xlim([min(cycle),max(cycle)])
p.set_xlabel('r') 
p.set_ylabel('g(r)')   
#plt.title('K$_2$CuF$_4$ I4/mmm to Cmca', fontsize=16)
#plt.legend()
#plt.hlines(0.0, min(cycle), max(cycle), colors='k', linestyles='dashed',
#        linewidth=1.0)
plt.show()
fig.savefig('gr.png')
