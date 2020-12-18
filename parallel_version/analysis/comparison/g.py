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
rp2 = []
gp2 = []
rp4 = []
gp4 = []
r1 = []
g1 = []
r1p2 = []
g1p2 = []
r1p4 = []
g1p4 = []
r2 = []
g2 = []
r2p2 = []
g2p2 = []
r2p4 = []
g2p4 = []

f = open('g_100_particles.out','r')
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

f = open('g_100_particles_2.out','r')
# Read title and comment lines
for i in range(0,13):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    rp2.append(float(l[0]))
    gp2.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_100_particles_4.out','r')
# Read title and comment lines
for i in range(0,13):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    rp4.append(float(l[0]))
    gp4.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_150_particles.out','r')
# Read title and comment lines
for i in range(0,5):
    line = f.readline()
# Store density and temperature to identify the result
line = f.readline()
l = line.split()
rho1 = float(l[5])
line = f.readline()
l = line.split()
T1 = float(l[3])
# Read title and comment lines
for i in range(0,6):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    r1.append(float(l[0]))
    g1.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_150_particles_2.out','r')
# Read title and comment lines
for i in range(0,13):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    r1p2.append(float(l[0]))
    g1p2.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_150_particles_4.out','r')
# Read title and comment lines
for i in range(0,13):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    r1p4.append(float(l[0]))
    g1p4.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_200_particles.out','r')
# Read title and comment lines
for i in range(0,5):
    line = f.readline()
# Store density and temperature to identify the result
line = f.readline()
l = line.split()
rho2 = float(l[5])
line = f.readline()
l = line.split()
T2 = float(l[3])
# Read title and comment lines
for i in range(0,6):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    r2.append(float(l[0]))
    g2.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_200_particles_2.out','r')
# Read title and comment lines
for i in range(0,13):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    r2p2.append(float(l[0]))
    g2p2.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_200_particles_4.out','r')
# Read title and comment lines
for i in range(0,13):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    r2p4.append(float(l[0]))
    g2p4.append(float(l[3]))
    line = f.readline()
f.close()

# Representation of g(r) against r
fig = plt.figure()
p = fig.add_subplot(111)
p.plot(r,g,'r', linewidth=1.5, label=r'$\rho$ = ' + str(rho) + ', T = ' +
        str(T))
p.plot(rp2,gp2,'r--', linewidth=1.5, label='2 procs')
p.plot(rp4,gp4,'ro', markersize=2.0, label='4 procs')
p.plot(r1,g1,'b', linewidth=1.5, label=r'$\rho$ = ' + str(rho1) + ', T = ' +
        str(T1))
p.plot(r1p2,g1p2,'b--', linewidth=1.5)
p.plot(r1p4,g1p4,'bo', markersize=2.0)
p.plot(r2,g2,'g', linewidth=1.5, label=r'$\rho$ = ' + str(rho2) + ', T = ' +
        str(T2))
p.plot(r2p2,g2p2,'g--', linewidth=1.5)
p.plot(r2p4,g2p4,'go', markersize=2.0)
plt.hlines(1.0,min(r),max(r), colors='k', linestyles='dashed', linewidth=1.5)
axes = plt.gca()
axes.set_xlim([min(r),max(r)])
p.set_xlabel('r/$\sigma$') 
p.set_ylabel('g(r)')   
plt.title('Pair correlation function', fontsize=16)
plt.legend()
plt.show()
fig.savefig('grs.png')
