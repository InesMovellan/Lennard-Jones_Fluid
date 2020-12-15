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
g100 = []
g150 = []
g200 = []
g250 = []

f = open('g_100_particles_4.out','r')
for i in range(0,7):
    line = f.readline()
line = f.readline()
while line:
    l = line.split()
    inc.append(float(l[0]))
    g100.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_150_particles_2.out','r')
for i in range(0,7):
    line = f.readline()
line = f.readline()
while line:
    l = line.split()
    g150.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_200_particles_2.out','r')
for i in range(0,7):
    line = f.readline()
line = f.readline()
while line:
    l = line.split()
    g200.append(float(l[3]))
    line = f.readline()
f.close()

#f = open('g_250_particles.out','r')
#for i in range(0,7):
#    line = f.readline()
#line = f.readline()
#while line:
#    l = line.split()
#    g250.append(float(l[3]))
#    line = f.readline()
#f.close()

#fit = np.polyfit(inc,g,10)
#inc1 = np.linspace(min(inc),max(inc))
#g1 = fit[0]*inc1**3+fit[1]*inc1**2+fit[2]*inc1+fit[3]
# Representation of the energy against the distortion coordinate
fig = plt.figure()
p = fig.add_subplot(111)
#p.plot(inc,g,'ro', markersize=4)
p.plot(inc,g100,'r', linewidth=2.0, label='100 atoms')
p.plot(inc,g150,'b', linewidth=2.0, label='150 atoms')
p.plot(inc,g200,'g', linewidth=2.0, label='200 atoms')
#p.plot(inc,g250,'m', linewidth=2.0, label='250 atoms')
plt.hlines(1.0,min(inc),max(inc), colors='k', linestyles='dashed')
axes = plt.gca()
axes.set_xlim([min(inc),max(inc)])
p.set_xlabel('r/$\sigma$') 
p.set_ylabel('g(r)')   
plt.title('Pair correlation function', fontsize=16)
plt.legend()
plt.show()
fig.savefig('grs.png')
