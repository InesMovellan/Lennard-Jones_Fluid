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
r1 = []
g1 = []
r2 = []
g2 = []
r3 = []
g3 = []
r4 = []
g4 = []

f = open('g_1.45_T.out','r')
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
f = open('g_1.268_T.out','r')
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

f = open('g_1.15_T.out','r')
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

f = open('g_1.0_T.out','r')
# Read title and comment lines
for i in range(0,5):
    line = f.readline()
# Store density and temperature to identify the result
line = f.readline()
l = line.split()
rho3 = float(l[5])
line = f.readline()
l = line.split()
T3 = float(l[3])
# Read title and comment lines
for i in range(0,6):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    r3.append(float(l[0]))
    g3.append(float(l[3]))
    line = f.readline()
f.close()

f = open('g_0.8_T.out','r')
# Read title and comment lines
for i in range(0,5):
    line = f.readline()
# Store density and temperature to identify the result
line = f.readline()
l = line.split()
rho4 = float(l[5])
line = f.readline()
l = line.split()
T4 = float(l[3])
# Read title and comment lines
for i in range(0,6):
    line = f.readline()
# Read the data from the file: r and g(r)
line = f.readline()
while line:
    l = line.split()
    r4.append(float(l[0]))
    g4.append(float(l[3]))
    line = f.readline()
f.close()


# Representation of g(r) against r
fig = plt.figure()
p = fig.add_subplot(111)
p.plot(r,g,'r', linewidth=2.0, label=r'$\rho$ = ' + str(rho) + ', T = ' +
        str(T))
p.plot(r1,g1,'b', linewidth=2.0, label=r'$\rho$ = ' + str(rho1) + ', T = ' +
        str(T1))
p.plot(r2,g2,'g', linewidth=2.0, label=r'$\rho$ = ' + str(rho2) + ', T = ' +
        str(T2))
p.plot(r3,g3,'m', linewidth=2.0, label=r'$\rho$ = ' + str(rho3) + ', T = ' +
        str(T3))
p.plot(r4,g4,'k', linewidth=2.0, label=r'$\rho$ = ' + str(rho4) + ', T = ' +
        str(T4))
plt.hlines(1.0,min(r),max(r), colors='k', linestyles='dashed', linewidth=1.5)
axes = plt.gca()
axes.set_xlim([min(r),max(r)])
p.set_xlabel('r/$\sigma$') 
p.set_ylabel('g(r)')   
plt.title('Pair correlation function', fontsize=16)
plt.legend()
plt.show()
fig.savefig('grs.png')
