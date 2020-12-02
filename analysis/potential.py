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

cycle = []
pot = []
f = open('v_out','r')
for i in range(0,5000):
    line = f.readline()
for i in range(0,95000):
    line = f.readline()
    l = line.split()
    cycle.append(float(l[0]))
    pot.append(float(l[1]))
f.close()

#line = f.readline()
#while line:
#   l = line.split()
#   cycle.append(float(l[0]))
#   pot.append(float(l[1]))
#   line = f.readline()
#f.close()
#for i in range(0,1000):
#    cycle.pop(i)
#    pot.pop(i)

# Representation of the energy against the distortion coordinate
fig = plt.figure()
p = fig.add_subplot(111)
#p.plot(d,E,'ro', markersize=4)
#p.plot(cycle,pot,'ro')
p.plot(cycle,pot,'r', linewidth=1.0)
axes = plt.gca()
#axes.set_xlim([min(d),max(d)])
p.set_xlabel('Monte Carlo cycle') 
p.set_ylabel('V (reduce units)')   
#plt.title('K$_2$CuF$_4$ I4/mmm to Cmca', fontsize=16)
#plt.legend()
#plt.hlines(0.0, min(cycle), max(cycle), colors='k', linestyles='dashed',
#        linewidth=1.0)
plt.show()
fig.savefig('MC_pot.png')
