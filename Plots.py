####  This file is for data processing and visualization of simulation output.


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re
from scipy.stats import linregress

####  1. Initialize arrays for initial conditions and pion field data.

param_file = open('inputs.txt','r')

content = param_file.readlines()

for line in content:
    if line.startswith('#'):
        content.remove(line)

c = content.count('\n')
for i in range(c):
    content.remove('\n')

for i in range(len(content)):
    var = content[i].split()[0]
    val = content[i].split()[2]
    if var == 'Omega_m':
        Omega_m = float(val)
    elif var == 'H0':
        H0 = float(val)
    elif var == 'boxsize':
        boxsize = float(val)
    elif var == 'meshpoints':
        meshpoints = int(val)
    elif var == 'sigma':
        sigma = float(val)
    elif var == 'z_init':
        z_init = float(val)
        a_init = 1/(1+z_init)
    elif var == 'a_init':
        a_init = float(val)
        z_init = 1/a_init - 1
    elif var == 'z_final':
        z_final = float(val)
        a_final = 1/(1+z_final)
    elif var == 'a_final':
        a_final = float(val)
        z_final = 1/a_final - 1
    elif var == 'initial_conditions_path':
        initial_conditions_path = val
    elif var == 'pion_output_path':
        pion_output_path = val
    elif var == 'pion_prime_output_path':
        pion_prime_output_path = val
    elif var == 'overdensity_output_path':
        overdensity_output_path = val
    elif var == 'plot_initial_conditions_Pk':
        plot_initial_conditions_Pk = int(val)
    elif var == 'plot_pion_Pk':
        plot_pion_Pk = int(val)
    elif var == 'plot_pion_Pk_NLO':
        plot_pion_Pk_NLO = int(val)
    elif var == 'plot_pion_prime_Pk':
        plot_pion_prime_Pk = int(val)
    elif var == 'plot_overdensity_Pk':
        plot_overdensity_Pk = int(val)
    elif var == 'plot_initial_conditions_3D':
        plot_initial_conditions_3D = int(val)
    elif var == 'plot_pion_3D':
        plot_pion_3D = int(val)
    elif var == 'plot_pion_prime_3D':
        plot_pion_prime_3D = int(val)
    elif var == 'plot_overdensity_3D':
        plot_overdensity_3D = int(val)

param_file.close()

####Simulation parameters -- must match initial and evolved conditions
##Omega_m = 0.31  ##matter density fraction in present epoch
##H0 = 68         ##km/s/Mpc
##z_init = 99     ##initial redshift
##
##boxsize = 600   ##Mpc
##meshpoints = 64 ##runtime scales as N^3 log N
##
##z_final = 0 ## This should be read after running Cosmic_Evolution.py
##
####Input data to be analyzed
##initial_conditions_path = 'initial_conditions.txt'
##pion_output_path = 'PidataC++.txt'
##pion_prime_output_path = 'PiPrimedataC++.txt'
##overdensity_output_path = 'DeltaC++.txt'
##
####2D power spectrum graph options
##plot_initial_conditions_Pk = 1
##plot_pion_Pk = 1
##plot_pion_Pk_NLO = 1
##plot_pion_prime_Pk = 0
##plot_overdensity_Pk = 0
##
####3D plot options
##plot_initial_conditions_3D = 0
##plot_pion_3D = 0
##plot_pion_prime_3D = 0
##plot_overdensity_3D = 0


#### 2. Read data from .txt files into 3D Python arrays.

ic = np.zeros((meshpoints, meshpoints, meshpoints))
pidata = np.zeros((meshpoints, meshpoints, meshpoints))
piprimedata = np.zeros((meshpoints, meshpoints, meshpoints))
overdensity = np.zeros((meshpoints, meshpoints, meshpoints))

file = open(initial_conditions_path,'r')

content = file.readlines()

arr = content[0].split("[[[")

while '' in arr:
    arr.remove('')

arr_i = arr[0].split("[[")

for i in range(len(arr_i)):
    line = arr_i[i].split('[')
    while '' in line:
        line.remove('')
    for j in range(len(line)):
        points = line[j].split(',')
        while '' in points:
            points.remove('')
        for k in range(len(points)):
            ic[i][j][k] = float(points[k].replace(']',''))

file.close()

file = open(pion_output_path,'r')

content = file.readlines()

arr = content[0].split("[[[")

while '' in arr:
    arr.remove('')

arr_i = arr[0].split("[[")

for i in range(len(arr_i)):
    line = arr_i[i].split('[')
    while '' in line:
        line.remove('')
    for j in range(len(line)):
        points = line[j].split(',')
        while '' in points:
            points.remove('')
        for k in range(len(points)):
            pidata[i][j][k] = float(points[k].replace(']',''))
    
file.close()

file = open(pion_prime_output_path,'r')

content = file.readlines()

arr = content[0].split("[[[")

while '' in arr:
    arr.remove('')

arr_i = arr[0].split("[[")

for i in range(len(arr_i)):
    line = arr_i[i].split('[')
    while '' in line:
        line.remove('')
    for j in range(len(line)):
        points = line[j].split(',')
        while '' in points:
            points.remove('')
        for k in range(len(points)):
            piprimedata[i][j][k] = float(points[k].replace(']',''))
    
file.close()

file = open(overdensity_output_path,'r')

content = file.readlines()

arr = content[0].split("[[[")

while '' in arr:
    arr.remove('')

arr_i = arr[0].split("[[")

for i in range(len(arr_i)):
    line = arr_i[i].split('[')
    while '' in line:
        line.remove('')
    for j in range(len(line)):
        points = line[j].split(',')
        while '' in points:
            points.remove('')
        for k in range(len(points)):
            overdensity[i][j][k] = float(points[k].replace(']',''))
    
file.close()

####  3. Calculates the power spectra.

meshsize = boxsize/meshpoints
ic_FT = np.fft.fftn(ic)
pion_FT = np.fft.fftn(pidata)
pion_prime_FT = np.fft.fftn(piprimedata)
overdensity_FT = np.fft.fftn(overdensity)

kvals = []
ic_Pkvals = []
pion_Pkvals = []
pion_prime_Pkvals = []
overdensity_Pkvals = []
bin_counts = []
factor = []
NLO_diff = []

for i in range(1,int(np.floor(meshpoints*np.sqrt(3)/2)+1)):
    kvals.append((6.28318/boxsize)*i)
    ic_Pkvals.append(0)
    pion_Pkvals.append(0)
    pion_prime_Pkvals.append(0)
    overdensity_Pkvals.append(0)
    bin_counts.append(0)
    factor.append(1)
    NLO_diff.append(0)

for i in range(int(np.floor(meshpoints/2))):
    for j in range(int(np.floor(meshpoints/2))):
        for k in range(int(np.floor(meshpoints/2))):
            if i*i +j*j + k*k > 0:
                bin_counts[int(np.sqrt(i*i + j*j + k*k))] += 1
                ic_Pkvals[int(np.sqrt(i*i + j*j + k*k))] += ((299792.458**2)*((meshsize**6)/(boxsize**3))*(np.abs(ic_FT[i][j][k])**2))
                pion_Pkvals[int(np.sqrt(i*i + j*j + k*k))] += ((299792.458**2)*((meshsize**6)/(boxsize**3))*(np.abs(pion_FT[i][j][k])**2))
                pion_prime_Pkvals[int(np.sqrt(i*i + j*j + k*k))] += ((299792.458**2)*((meshsize**6)/(boxsize**3))*(np.abs(pion_prime_FT[i][j][k])**2))
                overdensity_Pkvals[int(np.sqrt(i*i + j*j + k*k))] += (((meshsize**6)/(boxsize**3))*(np.abs(overdensity_FT[i][j][k])**2))

for i in range(len(kvals)):
    if bin_counts[i] != 0:
        ic_Pkvals[i] /= bin_counts[i]
        pion_Pkvals[i] /= bin_counts[i]
        pion_prime_Pkvals[i] /= bin_counts[i]
        overdensity_Pkvals[i] /= bin_counts[i]

del kvals[0]
del ic_Pkvals[0]
del pion_Pkvals[0]
del pion_prime_Pkvals[0]
del overdensity_Pkvals[0]
del factor[0]
del NLO_diff[0]

a_init = 1/(1 + z_init)
a_final = 1/(1+z_final)
lin_pi_growth_times_a2 = (2*np.sqrt(a_init**5))/(5*(H0/299792.458)*np.sqrt(Omega_m))
D_init = lin_pi_growth_times_a2/(a_init**2)
a_scale = a_init

for i in range(0, 1000000):
    if a_scale > a_final:
        break
    
    a_step = a_scale*0.01
    lin_pi_growth_times_a2 += a_step/((H0/299792.458)*np.sqrt(((Omega_m/(a_scale**3))+(1-Omega_m))))
    a_scale += a_step

D_final = lin_pi_growth_times_a2/(a_scale**2)

for i in range(len(kvals)):
    if ic_Pkvals[i] > 1E-12:
        factor[i] = (pion_Pkvals[i]/ic_Pkvals[i])
        NLO_diff[i] = ((pion_Pkvals[i] - ((D_final/D_init)**2)*ic_Pkvals[i])/(((D_final/D_init)**2)*ic_Pkvals[i]))
        #NLO_diff[i] = (pion_Pkvals[i] - (pion_Pkvals[0]/ic_Pkvals[0])*ic_Pkvals[i])/((pion_Pkvals[0]/ic_Pkvals[0])*ic_Pkvals[i])
####  4. 2D power spectra graphs.

if plot_initial_conditions_Pk ==1 or plot_pion_Pk == 1:
    
    if plot_initial_conditions_Pk == 1:

        plt.scatter(kvals, ic_Pkvals, alpha=1)
        plt.title("Initial conditions power spectrum")
        plt.xlabel("k (Mpc^-1)")
        plt.xlim(6.28318/boxsize,(6.28318/boxsize)*(meshpoints/2))
        plt.ylabel("P(k) ((km/s)^2)(Mpc^5)")
        plt.xscale("log")
        plt.yscale("log")

    if plot_pion_Pk == 1:

        plt.scatter(kvals, pion_Pkvals, alpha=1)
        plt.title("Pion field power spectrum")
        plt.xlabel("k (Mpc^-1)")
        plt.xlim(6.28318/boxsize,(6.28318/boxsize)*(meshpoints/2))
        plt.ylim(np.min(pion_Pkvals)/10, np.max(pion_Pkvals)*10)
        plt.ylabel("P(k) ((km/s)^2)(Mpc^5)")
        plt.xscale("log")
        plt.yscale("log")

    plt.show()

if plot_pion_prime_Pk == 1:
    
    plt.scatter(kvals, pion_prime_Pkvals, alpha=1)
    plt.title("Pion prime field power spectrum")
    plt.xlabel("k (Mpc^-1)")
    plt.xlim(6.28318/boxsize,(6.28318/boxsize)*(meshpoints/2))
    plt.ylim(np.min(pion_prime_Pkvals)/10, np.max(pion_prime_Pkvals)*10)
    plt.ylabel("P(k) ((km/s)^2)(Mpc^3)")
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

if plot_overdensity_Pk == 1:
    
    plt.scatter(kvals, overdensity_Pkvals, alpha=1)
    plt.title("Overdensity power spectrum")
    plt.xlabel("k (Mpc^-1)")
    plt.xlim(6.28318/boxsize,(6.28318/boxsize)*(meshpoints/2))
    plt.ylim(np.min(overdensity_Pkvals)/10, np.max(overdensity_Pkvals)*10)
    plt.ylabel("P(k) (Mpc^3)")
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

if plot_pion_Pk_NLO == 1:

    plt.scatter(kvals, NLO_diff, alpha=1)
    plt.title("Fractional size of NLO corrections")
    plt.xlabel("k")
    plt.xlim(6.28318/boxsize,(6.28318/boxsize)*(meshpoints/2))
    plt.ylim(-0.5,0.5)
    plt.ylabel("(P(k)-P_0(k))/P_0(k)")
    plt.xscale("log")
    plt.yscale("linear")
    plt.show()

####  5. 3D plots.


mesh = np.zeros((meshpoints, meshpoints, meshpoints, 3))

for i in range(meshpoints):
    for j in range(meshpoints):
        for k in range(meshpoints):
            mesh[i][j][k] = [meshsize*i, meshsize*j, meshsize*k]

if plot_initial_conditions_3D == 1:

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cax = ax.scatter(mesh[:, :, :, 0].flatten(), mesh[:, :, :, 1].flatten(), mesh[:, :, :, 2].flatten(), c=ic.flatten(), cmap='jet', alpha=1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.suptitle('Pion field')
    plt.title('Initial conditions at a = ' + str(a_init))
    color_bar = plt.colorbar(cax)
    color_bar.set_alpha(1)
    color_bar.draw_all()
    plt.show()

if plot_pion_3D == 1:

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cax = ax.scatter(mesh[:, :, :, 0].flatten(), mesh[:, :, :, 1].flatten(), mesh[:, :, :, 2].flatten(), c=pidata.flatten(), cmap='jet', alpha=1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.suptitle('Pion field')
    plt.title('Final configuration at a = ' + str(a_final))
    color_bar = plt.colorbar(cax)
    color_bar.set_alpha(1)
    color_bar.draw_all()
    plt.show()

if plot_pion_prime_3D == 1:

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cax = ax.scatter(mesh[:, :, :, 0].flatten(), mesh[:, :, :, 1].flatten(), mesh[:, :, :, 2].flatten(), c=piprimedata.flatten(), cmap='jet', alpha=1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.suptitle('Pion primefield')
    plt.title('Final configuration at a = ' + str(a_final))
    color_bar = plt.colorbar(cax)
    color_bar.set_alpha(1)
    color_bar.draw_all()
    plt.show()
    
if plot_overdensity_3D == 1:

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cax = ax.scatter(mesh[:, :, :, 0].flatten(), mesh[:, :, :, 1].flatten(), mesh[:, :, :, 2].flatten(), c=overdensity.flatten(), cmap='jet', alpha=1)#(overdensity.flatten()-np.min(overdensity.flatten()))/(np.max(overdensity.flatten())-np.min(overdensity.flatten())))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.suptitle('Overdensity field')
    plt.title('Final configuration at a = ' + str(a_final))
    color_bar = plt.colorbar(cax)
    color_bar.set_alpha(1)
    color_bar.draw_all()
    plt.show()

