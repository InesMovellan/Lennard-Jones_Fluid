import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc                         
from matplotlib.ticker import FormatStrFormatter  
from scipy.interpolate import interp1d            

font = {'family' : 'serif',   
        'weight' : 'normal',  
        'size' : '14'         
       }                      
plt.rc('font', **font)        

r = []
g = []

f = open('g.out','r')
# Read title and comment lines
for i in range(0,5):
    line = f.readline()
# Store density and temperature to identify the result
line = f.readline()
l = line.split()
rho = float(l[5])
line = f.readline()
l = line.split()
T = float(l[3])
# Read title and comment lines
for i in range(0,6):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    r.append(float(l[0]))
    g.append(float(l[3]))
    line = f.readline()
f.close()

# Representation of g(r) against r
fig = plt.figure()
p = fig.add_subplot(111)
p.plot(r,g,'r', linewidth=2.0, label=r'$\rho$ = ' + str(rho) + ', T = ' +
        str(T))
plt.hlines(1.0,min(r),max(r), colors='k', linestyles='dashed', linewidth=1.5)
axes = plt.gca()
axes.set_xlim([min(r),max(r)])
p.set_xlabel('r/$\sigma$') 
p.set_ylabel('g(r)')   
plt.title('Pair correlation function', fontsize=16)
plt.legend()
plt.show()
fig.savefig('gr.png')
