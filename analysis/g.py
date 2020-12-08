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
f = open('g.out','r')
for i in range(0,7):
    line = f.readline()

line = f.readline()
while line:
    l = line.split()
    inc.append(float(l[0]))
    #g.append(float(l[1]))
    g.append(float(l[3]))
    line = f.readline()
f.close()

#n = np.asarray(n)
#n = n/max(n)

#fit = np.polyfit(inc,g,10)
#inc1 = np.linspace(min(inc),max(inc))
#g1 = fit[0]*inc1**3+fit[1]*inc1**2+fit[2]*inc1+fit[3]
# Representation of the energy against the distortion coordinate
fig = plt.figure()
p = fig.add_subplot(111)
p.plot(inc,g,'ro', markersize=5)
p.plot(inc,g,'r', linewidth=2.0)
axes = plt.gca()
#axes.set_xlim([min(cycle),max(cycle)])
p.set_xlabel('r/$\sigma$') 
p.set_ylabel('g(r)')   
plt.title('Pair correlation function', fontsize=16)
plt.show()
fig.savefig('gr.png')
