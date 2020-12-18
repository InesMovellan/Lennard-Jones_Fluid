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

f = open('V.out','r')
# Read the title and comment lines
for i in range(0,11):
    line = f.readline()
# Uncomment to skip the first values of the potential energy
for i in range(0,2):
    line = f.readline()

# Read the potential energy for each MC cycle sampled and store them
line = f.readline()
while line:
    l = line.split()
    cycle.append(float(l[0]))
    pot.append(float(l[1]))
    line = f.readline()
f.close()

# Representation of the potential energy again the MC cycle sampled
fig = plt.figure()
p = fig.add_subplot(111)
#p.plot(cycle,pot,'ro')
p.plot(cycle,pot,'r', linewidth=1.5)
axes = plt.gca()
axes.set_xlim([min(cycle),max(cycle)])
p.set_xlabel('Monte Carlo cycle') 
p.set_ylabel('V (reduce units)')   
plt.title('Potential energy', fontsize=16)
plt.show()
fig.savefig('V_vs_cycleMC.png')
