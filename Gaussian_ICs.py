
####This file is for generating Gaussian initial conditions for the pion field
####in a Lambda CDM universe.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re


## 1. Input simulation parameters

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
    elif var == 'Omega_b':
        Omega_b = float(val)
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
    elif var == 'Pk_fitting_option':
        Pk_fitting_option = int(val)
    elif var == 'Pk_input_file':
        Pk_input_file = val
    elif var == 'n_tilt':
        n_tilt = float(val)
    elif var == 'T_CMB':
        T_CMB = float(val)
    elif var == 'initial_conditions_path':
        output_file_path = val
        
param_file.close()

##Uncomment this block to change parameters from within Python

####Cosmic
##Omega_m = 0.31  ##matter density fraction in present epoch
##H0 = 68         ##km/s/Mpc
##z_init = 99     ##initial redshift
##
##boxsize = 600   ##Mpc
##meshpoints = 64 ##runtime scales as N^3 log N
##sigma = 25
##
###Matter power spectrum: set Pk_fitting_option = 0 to use BBKS fitting formula or 1 to use a CLASS output file.
##Pk_fitting_option = 1
##n_tilt = 0.96                               ## if Pk_fitting_option = 0
##Pk_input_file = 'mPk_from_class.dat'            ## if Pk_fitting_option = 1
##
##output_file_path = 'initial_conditions.txt'


print("Generating initial conditions at z = " + str(z_init) + " with Omega_m = " + str(Omega_m) + " and H0 = " + str(H0) + "km/s/Mpc.\n") 
print("Box Size is " + str(boxsize) + " Mpc and grid contains " + str(meshpoints) + "^3 mesh points.\n")


## 2. Calculate the power spectrum.

kvals = []
Pkvals = []
h = H0/100

if Pk_fitting_option == 0:
    
    file = open(Pk_input_file,'r')

    content = file.readlines()

    for line in content:
        if '#' not in line:
            numbers = line.split()
            kvals.append(float(numbers[0])*h)
            Pkvals.append(float(numbers[1])/h**3)
            
    file.close()

    k_max = kvals[len(kvals)-1]
    k_min = kvals[0]

    print('Pk read from k_min = ' + str(k_min) + ' Mpc^{-1} to k_max = ' + str(k_max) + ' Mpc^{-1}.\n') 

elif Pk_fitting_option == 1:
    x_var = 0.328*np.log(431*Omega_m*h**2)*(Omega_b/Omega_m)-0.38*np.log(22.3*Omega_m*h**2)*(Omega_b/Omega_m)**2
    y_var = (0.43*44.5*np.log(9.83/(Omega_m*h**2)))/np.sqrt(1+10*((Omega_m*h**2)**(3/4)))
    kvals.append(6.28318/1E6)
    for i in range(1,int(np.floor(meshpoints*np.sqrt(3)/2))):
        kvals.append(6.28318*i/boxsize)
    for i in range(0,len(kvals)):
        k_mag = kvals[i]
        q = ((T_CMB/2.7)**2)*k_mag*(1+(k_mag*y_var)**4)/((Omega_m*h**2)*(1+(1-x_var)*(k_mag*y_var)**4))
        T = np.log(2*np.e+1.8*q)/(np.log(2*np.e+1.8*q)+(14.2+731/(1+62.5*q))*q**2)
        Pkval = (2*3.14159**2)*((299792.458/H0)**(n_tilt+3))*(k_mag**n_tilt)
        Pkval *= ((2.18E-5)*(Omega_m**(-0.785-0.05*np.log(Omega_m)))*np.exp(-0.95*(n_tilt-1)-0.169*(n_tilt-1)**2))**2
        Pkval *= T**2
        Pkvals.append(Pkval)
    print('Calculating power spectrum from Eisentein and Hu fitting formula.\n')

elif Pk_fitting_option == 2:
    kvals.append(6.28318/1E6)
    for i in range(1,int(np.floor(meshpoints*np.sqrt(3)/2))):
        kvals.append(6.28318*i/boxsize)
    for i in range(0,len(kvals)):
        k_mag = kvals[i]
        q = ((T_CMB/2.7)**2)*k_mag/(Omega_m*(h**2)*np.exp(-Omega_b-Omega_b/Omega_m))
        T = (np.log(1+2.34*q)/(2.34*q))*(1+(8.89*q)+(14.1*q)**2 + (5.46*q)**3 + (6.07*q)**4)**(-1/4)
        Pkval = (2*3.14159**2)*((299792.458/H0)**(n_tilt+3))*(k_mag**n_tilt)
        Pkval *= ((1.94E-5)*(Omega_m**(-0.785-0.05*np.log(Omega_m)))*np.exp(-0.95*(n_tilt-1)-0.169*(n_tilt-1)**2))**2
        Pkval *= T**2
        Pkvals.append(Pkval)
    print('Calculating power spectrum from BBKS fitting formula.\n')

    
## 3. Solve for linear growth factor

a_init = 1/(1 + z_init)
Hubble_init = (H0/299792.458)*a_init*np.sqrt((Omega_m/a_init**3))
time_init = 2*np.sqrt(a_init)*(299792.458)/(H0*np.sqrt(Omega_m))
a_scale = a_init
lin_pi_growth_times_a2 = (2*np.sqrt(a_init**5))/(5*np.sqrt(Omega_m))
lin_delta_init = (2*a_init)/(5*Omega_m)
lin_delta_growth = lin_delta_init

for i in range(0, 1000000):
    if a_scale > 1/(1+z_final):
        break
    
    a_step = a_scale*0.0001
    lin_pi_growth_times_a2 += a_step/np.sqrt(((Omega_m/(a_scale**3))+(1-Omega_m)))
    lin_delta_growth += (a_step*lin_pi_growth_times_a2)/((a_scale**4)*np.sqrt(((Omega_m/(a_scale**3))+(1-Omega_m))))
    a_scale += a_step
    
print("Linear growth factor for overdensity between z_init = " + str(z_init) + " and z =  " + str(z_final) + " is " + str(lin_delta_growth/lin_delta_init) + "\n")


## 4. Generate Gaussian initial conditions for pion field power spectrum initial conditions

pi_Pkvals = []

if Pk_fitting_option == 0:
    for i in range(0,len(kvals)):
        pi_Pkvals.append((((lin_delta_init/lin_delta_growth)*(2/time_init)/(kvals[i]**2))**2)*Pkvals[i])

elif Pk_fitting_option == 1:
    #pi_Pkvals.append(0)
    for i in range(0,len(kvals)):
        k_mag = kvals[i]
        pi_Pkvals.append((((lin_delta_init/lin_delta_growth)*(2/time_init)/(kvals[i]**2))**2)*Pkvals[i])

elif Pk_fitting_option == 2:
    #pi_Pkvals.append(0)
    for i in range(0,len(kvals)):
        k_mag = kvals[i]
        pi_Pkvals.append((((lin_delta_init/lin_delta_growth)*(2/time_init)/(kvals[i]**2))**2)*Pkvals[i])


meshsize = boxsize/meshpoints

mesh = np.zeros((meshpoints, meshpoints, meshpoints, 3))

for i in range(meshpoints):
    for j in range(meshpoints):
        for k in range(meshpoints):
            mesh[i][j][k] = [meshsize*i, meshsize*j, meshsize*k]

white_noise = np.zeros((meshpoints, meshpoints, meshpoints))

for i in range(meshpoints):
    for j in range(meshpoints):
        for k in range(meshpoints):
             white_noise[i][j][k] = np.random.normal(0,1)/meshsize**(3/2)

momentum_ic = np.fft.fftn(white_noise)
momentum_icprime = np.zeros((meshpoints, meshpoints, meshpoints))

for i in range(meshpoints):
    for j in range(meshpoints):
        for k in range(meshpoints):
             k_mag = (6.28318/boxsize)*np.sqrt(min(i, meshpoints-i)**2 + min(j,meshpoints-j)**2 + min(k, meshpoints-k)**2 + 1E-24)
             if k_mag > max(kvals):
                 momentum_ic[i][j][k] *= 0
             else:
                 momentum_ic[i][j][k] *= np.sqrt(np.exp(np.interp(k_mag, kvals, np.log(pi_Pkvals))))
##                if Pk_fitting_option == 1:
##                     momentum_ic[i][j][k] *= np.sqrt(np.interp(k_mag, kvals, pi_Pkvals))
##                elif Pk_fitting_option == 0:
##                    q = k_mag/(Omega_m*h**2)
##                    T = (np.log(1+2.34*q)/(2.34*q))*(1+(8.89*q)+(14.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**(-1/4)
##                    pi_Pkval = (2*3.14159**2)*((299792.458/H0)**(n_tilt+3))*(k_mag**n_tilt)
##                    pi_Pkval *= ((1.94E-5)*(Omega_m**(-0.785-0.05*np.log(Omega_m)))*np.exp(-0.95*(n_tilt-1)-0.169*(n_tilt-1)**2))**2
##                    pi_Pkval *= T**2
##                    pi_Pkval *= (((lin_delta_init/lin_delta_growth)*(2/time_init)/(k_mag**2))**2)
##                    momentum_ic[i][j][k] *= np.sqrt(pi_Pkval)
                                                                                         
momentum_ic[0][0][0] = 0

momentum_icprime = momentum_ic/time_init

##Gaussian smoothing

if sigma > 0:
    k_cutoff = 6.28318/sigma
    for i in range(0, meshpoints):
        for j in range(0, meshpoints):
            for k in range(0, meshpoints):
                k_mag = (6.28318/boxsize)*np.sqrt(min(i, meshpoints-i)**2 + min(j,meshpoints-j)**2 + min(k, meshpoints-k)**2 + 1E-24)
                momentum_ic[i][j][k] *= np.exp(- (k_mag**2)/(2*k_cutoff**2))
                momentum_icprime[i][j][k] *= np.exp(- (k_mag**2)/(2*k_cutoff**2))

ic = np.real(np.fft.ifftn(momentum_ic))
icprime = np.real(np.fft.ifftn(momentum_icprime))


## 5. Write initial conditions to txt file by converting the nested list into a string.

file1 = open(output_file_path,"w")
file1.write("[")
for i in range(0, meshpoints-1):
    file1.write("[")
    for j in range(0, meshpoints-1):
        file1.write("[")
        for k in range(0, meshpoints-1):
            file1.write(str(ic[i][j][k])+",")
        file1.write(str(ic[i][j][meshpoints-1]))
        file1.write("],")
    file1.write("[")
    for k in range(0, meshpoints-1):
        file1.write(str(ic[i][meshpoints-1][k])+",")
    file1.write(str(ic[i][meshpoints-1][meshpoints-1]))
    file1.write("],")
file1.write("[")
for j in range(0, meshpoints-1):
    file1.write("[")
    for k in range(0, meshpoints-1):
        file1.write(str(ic[meshpoints-1][j][k])+",")
    file1.write(str(ic[meshpoints-1][j][meshpoints-1]))
    file1.write("],")
file1.write("[")
for k in range(0, meshpoints-1):
    file1.write(str(ic[meshpoints-1][meshpoints-1][k])+",")
file1.write(str(ic[meshpoints-1][meshpoints-1][meshpoints-1]))
file1.write("]]]")
file1.close()

print("Initial conditions written to file " + str(output_file_path))


### 6. Extras: plots

####Plot the matter power spectrum at present epoch
##
##plt.scatter(kvals, Pkvals, alpha=1)
##plt.title("Matter power spectrum at a = 1")
##plt.xlabel("k (Mpc)")
##plt.ylabel("P(k) (Mpc)^3")
##plt.xscale("log")
##plt.yscale("log")
##plt.show()         

####Plot the initial conditions pion power spectrum
##
##plt.scatter(kvals, pi_Pkvals, alpha=1)
##plt.title("Pion field power spectrum at a = " + str(a_init))
##plt.xlabel("k (Mpc)")
##plt.ylabel("P(k) (Mpc)^5")
##plt.xscale("log")
##plt.yscale("log")
##plt.show()

####Plot the pion field initial conditions in 3D (This can take time for a large number of grid points.)

##fig = plt.figure()
##ax = fig.add_subplot(111, projection='3d')
##cax = ax.scatter(mesh[:, :, :, 0].flatten(), mesh[:, :, :, 1].flatten(), mesh[:, :, :, 2].flatten(), c=ic.flatten(), cmap='jet')
##ax.set_xlabel('X')
##ax.set_ylabel('Y')
##ax.set_zlabel('Z')
##plt.suptitle('Pion field initial conditions in 3D')
##color_bar = plt.colorbar(cax)
##color_bar.set_alpha(1)
##color_bar.draw_all()
##plt.show()
