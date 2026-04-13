####  This file calculates the evolution of the cosmic pion field in a Lambda CDM universe,
####  starting from the initial conditions.

import numpy as np
from numpy import random
import matplotlib.pyplot as plt
import time as myTime
from mpl_toolkits.mplot3d import Axes3D
import pickle

#### 1. Input cosmic and simulation parameters.

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
    elif var == 'cs2_0':
        cs2_0 = float(val)
    elif var == 'cv2_0':
        cv2_0 = float(val)
    elif var == 'Pk_fitting_option':
        Pk_fitting_option = int(val)
    elif var == 'initial_conditions_path':
        initial_conditions_path = val
    elif var == 'Pk_input_file':
        Pk_input_file = val
    elif var == 'n_tilt':
        n_tilt = float(val)
    elif var == 'auto_stop_enabled':
        auto_stop_enabled = int(val)
    elif var == 'max_steps':
        max_steps = int(val)
    elif var == 'step_factor':
        step_factor = float(val)
    elif var == 'show_initial_conditions':
        show_initial_conditions = int(val)
    elif var == 'show_final_configuration':
        show_final_configuration = int(val)
    elif var == 'pion_output_path':
        pion_output_path = val
    elif var == 'pion_prime_output_path':
        pion_prime_output_path = val
    elif var == 'overdensity_output_path':
        overdensity_output_path = val

param_file.close()

####Cosmic
##Omega_m = 0.31  ##matter density fraction in present epoch
##H0 = 68         ##km/s/Mpc
##z_init = 99     ##initial redshift
##
##boxsize = 600   ##Mpc
##meshpoints = 64 ##runtime scales as N^3 log N
##
####Input and ouput
##initial_conditions_path = 'initial_conditions.txt'
##show_initial_conditions = 0  ## Set to 1 to plot the initial conditions in 3D.
##pion_output_path = 'output.txt'
##pion_prime_output_path = 'output_prime.txt'
##overdensity_output_path = 'output_overdensity.txt'
##show_final_configuration = 1
##
####Simulation
##sigma = 0.01       ##Smoothing on small scales
##z_final = 0
##auto_stop_enabled = 1   ##Set to 1 to automatically stop if delta < -1 detected
##step_factor = 0.01      ##a is multiplied by this factor at each step
##max_steps = 6000
##
####Set the values of the EFT parameters in the present day.  Time dependence is assumed to be
####D^2 and D^2/(Omega_M**(1/2)), respectively.
##cs2_0 = -0.5E-6
##cv2_0 = 2.2E-6



print("Boxsize is: " + str(boxsize) + " Mpc and number of grid points is:\n")
print("Number of grid points is:", meshpoints)

#### 2. Initialize initial conditions from the given .txt file.

ic = np.zeros((meshpoints, meshpoints, meshpoints))
icprime = np.zeros((meshpoints, meshpoints, meshpoints))

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

a_init = 1/(1 + z_init)
Hubble_init = (H0/299792.458)*a_init*np.sqrt((Omega_m/a_init**3))
time_init = 2*np.sqrt(a_init)*(299792.458)/(H0*np.sqrt(Omega_m))
time = time_init

for i in range(meshpoints):
    for j in range(meshpoints):
        for k in range(meshpoints):
             icprime[i][j][k] = ic[i][j][k]/time_init

mesh = np.zeros((meshpoints, meshpoints, meshpoints, 3))
meshsize = boxsize/meshpoints

for i in range(meshpoints):
    for j in range(meshpoints):
        for k in range(meshpoints):
            mesh[i][j][k] = [meshsize*i, meshsize*j, meshsize*k]

if show_initial_conditions == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cax = ax.scatter(mesh[:, :, :, 0].flatten(), mesh[:, :, :, 1].flatten(), mesh[:, :, :, 2].flatten(), c=ic.flatten(), cmap='jet', alpha=1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.suptitle('Pion field initial conditions in 3D')
    plt.title('Conformal time = ' + str(time_init) + ', a = ' + str(a_init))
    color_bar = plt.colorbar(cax)
    color_bar.set_alpha(1)
    fig.draw_without_rendering()
    plt.show()

print("Initial conditions at z = " + str(z_init) + " from " + str(initial_conditions_path) + "\n")

####  3. Evolution of the equation of motion for the pion field.

a_scale = a_init
a_final = 1/(1+z_final)
lin_pi_growth_times_a2 = (2*np.sqrt(a_init**5))/(5*(H0/299792.458)*np.sqrt(Omega_m))

for i in range(0, 1000000):
    if a_scale > 1.0:
        break
    
    a_step = a_scale*0.0001
    lin_pi_growth_times_a2 += a_step/((H0/299792.458)*np.sqrt(((Omega_m/(a_scale**3))+(1-Omega_m))))
    a_scale += a_step

D_zero = lin_pi_growth_times_a2

print("Starting the run for the evolution of the pion field.  If the number of mesh points is large, this may take awhile...\n")
pidata = ic.copy()
time = time_init
a_scale = a_init
lin_pi_growth_times_a2 = (2*np.sqrt(a_init**5))/(5*(H0/299792.458)*np.sqrt(Omega_m))
numsteps = 0
pitimedata = []
piprimetimedata = []
overdensity = np.zeros((meshpoints,meshpoints,meshpoints))


piprimedata = icprime.copy()
ic, icprime

prog_start = myTime.time()

for t in range(0, max_steps):

    if a_scale > a_final:
        break
    
    a_step = a_scale*step_factor
    timestep = (a_step*299792.458)/((H0)*(a_scale**2)*np.sqrt((Omega_m/(a_scale**3))+(1-Omega_m)))
    Hubble = (H0/299792.458)*a_scale*np.sqrt((Omega_m/(a_scale**3))+(1-Omega_m))
    Omega_M = (Omega_m/(Omega_m + (1-Omega_m)*a_scale**3))
    Hubble_prime = (1-((3/2)*Omega_M))*(Hubble**2)
    D_growth = lin_pi_growth_times_a2/(a_scale**2)
    tildec2_s = (1.5*Omega_M - 3.5)*((D_growth/D_zero)**2)*np.sqrt(Omega_m/Omega_M)*cv2_0 + (2*D_growth/(Hubble*D_zero**2))*np.sqrt(Omega_m/Omega_M)*cv2_0;
    tildec2_s += ((D_growth/D_zero)**2)*cs2_0 + (2/(3*Omega_M))*((2*D_growth/(Hubble*D_zero**2)) - 3*((D_growth/D_zero)**2))*cs2_0
    tildec2_v_overH = (1/Hubble)*((D_growth/D_zero)**2)*np.sqrt(Omega_m/Omega_M)*cv2_0
    tildec2_v_overH += (1/Hubble)*((2/(3*Omega_M))*((-3*(D_growth/D_zero)**2)+(2*D_growth/(Hubble*D_zero**2))))*cs2_0

    gradsquared = np.zeros((meshpoints,meshpoints,meshpoints))
    gradpipiprime = np.zeros((meshpoints,meshpoints,meshpoints))
    

    pi_FT = np.fft.fftn(pidata)
    pi_prime_FT = np.fft.fftn(piprimedata)
    vx_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)
    vy_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)
    vz_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)
    vx_prime_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)
    vy_prime_FT = np.zeros((meshpoints,meshpoints,meshpoints), dtype = complex)
    vz_prime_FT = np.zeros((meshpoints,meshpoints,meshpoints), dtype = complex)
    

    for i in range(0, meshpoints):
        kx = (2*np.pi/boxsize)*i
        if i > meshpoints/2:
            kx -= (2*np.pi/boxsize)*meshpoints
        for j in range(0, meshpoints):
                ky = (2*np.pi/boxsize)*j
                if j > meshpoints/2:
                    ky -= (2*np.pi/boxsize)*meshpoints
                for k in range(0, meshpoints):
                    kz = (2*np.pi/boxsize)*k
                    if k > meshpoints/2:
                        kz -= (2*np.pi/boxsize)*meshpoints
                    k_mag = (6.28318/boxsize)*np.sqrt(min(i, meshpoints-i)**2 + min(j,meshpoints-j)**2 + min(k, meshpoints-k)**2)
                    vx_FT[i][j][k] = -complex(0,1)*kx*pi_FT[i][j][k]
                    vy_FT[i][j][k] = -complex(0,1)*ky*pi_FT[i][j][k]
                    vz_FT[i][j][k] = -complex(0,1)*kz*pi_FT[i][j][k]
                    vx_prime_FT[i][j][k] = -complex(0,1)*kx*pi_prime_FT[i][j][k]
                    vy_prime_FT[i][j][k] = -complex(0,1)*ky*pi_prime_FT[i][j][k]
                    vz_prime_FT[i][j][k] = -complex(0,1)*kz*pi_prime_FT[i][j][k]
                    
    vx = np.real(np.fft.ifftn(vx_FT))
    vy = np.real(np.fft.ifftn(vy_FT))
    vz = np.real(np.fft.ifftn(vz_FT))
    vx_prime = np.real(np.fft.ifftn(vx_prime_FT))
    vy_prime = np.real(np.fft.ifftn(vy_prime_FT))
    vz_prime = np.real(np.fft.ifftn(vz_prime_FT))
    del vx_prime_FT, vy_prime_FT, vz_prime_FT

    ##Calculating the nonlinear terms
    
    for i in range(0, meshpoints):
        for j in range(0, meshpoints):
            for k in range(0, meshpoints):
                gradsquared[i][j][k] += vx[i][j][k]**2 + vy[i][j][k]**2 + vz[i][j][k]**2
                gradpipiprime[i][j][k] += vx[i][j][k]*vx_prime[i][j][k] + vy[i][j][k]*vy_prime[i][j][k] + vz[i][j][k]*vz_prime[i][j][k]

    NewtonPhi = np.zeros((meshpoints,meshpoints,meshpoints))

    for i in range(0, meshpoints):
        for j in range(0, meshpoints):
            for k in range(0, meshpoints):
                NewtonPhi[i][j][k] += -(piprimedata[i][j][k]+(Hubble)*pidata[i][j][k]+(1/2)*gradsquared[i][j][k])


    NewtonPhi_FT = np.fft.fftn(NewtonPhi)
   
    delta_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)
    theta_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)
    thetaprime_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)


    for i in range(0, meshpoints):
        for j in range(0, meshpoints):
            for k in range(0, meshpoints):
                k_mag = (6.28318/boxsize)*np.sqrt(min(i, meshpoints-i)**2 + min(j,meshpoints-j)**2 + min(k, meshpoints-k)**2)
                delta_FT[i][j][k] = NewtonPhi_FT[i][j][k]*(k_mag**2)
                theta_FT[i][j][k] = -pi_FT[i][j][k]*(k_mag**2)
                thetaprime_FT[i][j][k] = -pi_prime_FT[i][j][k]*(k_mag**2)
                
    delta = np.real(np.fft.ifftn(delta_FT))
    overdensity = delta/(Hubble_prime - Hubble*Hubble)##Note that the variable delta here is proportional but not equal to the overdensity.
    theta = np.real(np.fft.ifftn(theta_FT))
    thetaprime = np.real(np.fft.ifftn(thetaprime_FT))
    #del pi_FT, pi_prime_FT, theta_FT, thetaprime_FT

    if auto_stop_enabled == 1:
        if min(overdensity.flatten()) < -1.0:
            print('Error: underdensity calculated to be less than -1.0 at time ' + str(time) + '\n')
            print('Max overdensity is: ' + str(max(overdensity.flatten())) + '\n')
            print("Stopping the run.\n")
            break

    gradx_delta_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)
    grady_delta_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)
    gradz_delta_FT = np.zeros((meshpoints,meshpoints,meshpoints),dtype = complex)

    for i in range(0, meshpoints):
        kx = (2*np.pi/boxsize)*i
        if i > meshpoints/2:
            kx -= (2*np.pi/boxsize)*meshpoints
        for j in range(0, meshpoints):
                ky = (2*np.pi/boxsize)*j
                if j > meshpoints/2:
                    ky -= (2*np.pi/boxsize)*meshpoints
                for k in range(0, meshpoints):
                    kz = (2*np.pi/boxsize)*k
                    if k > meshpoints/2:
                        kz -= (2*np.pi/boxsize)*meshpoints
                    gradx_delta_FT[i][j][k] = -complex(0,1)*kx*delta_FT[i][j][k]
                    grady_delta_FT[i][j][k] = -complex(0,1)*ky*delta_FT[i][j][k]
                    gradz_delta_FT[i][j][k] = -complex(0,1)*kz*delta_FT[i][j][k]

    gradx_delta = np.real(np.fft.ifftn(gradx_delta_FT))
    grady_delta = np.real(np.fft.ifftn(grady_delta_FT))
    gradz_delta = np.real(np.fft.ifftn(gradz_delta_FT))
        
    integrand = np.zeros((meshpoints,meshpoints,meshpoints))
                
    for i in range(0, meshpoints):
        for j in range(0, meshpoints):
            for k in range(0, meshpoints):

                integrand[i][j][k] = 0

                #integrand[i][j][k] += (Hubble_prime-Hubble**2)*theta[i][j][k]

                integrand[i][j][k] += (theta[i][j][k])*(delta[i][j][k])
                
                integrand[i][j][k] += vx[i][j][k]*gradx_delta[i][j][k] + vy[i][j][k]*grady_delta[i][j][k] + vz[i][j][k]*gradz_delta[i][j][k]
                
    FT = np.fft.fftn(integrand)

    FTinverseLaplacian = np.zeros((meshpoints,meshpoints,meshpoints),dtype = 'complex_')

    for i in range(0, meshpoints):
        for j in range(0, meshpoints):
            for k in range(0, meshpoints):
                if i*i + j*j + k*k > 0:
                    k_mag = (6.28318/boxsize)*np.sqrt(min(i, meshpoints-i)**2 + min(j,meshpoints-j)**2 + min(k, meshpoints-k)**2)
                    FTinverseLaplacian[i][j][k] = -FT[i][j][k]/(k_mag**2)
                #FTinverseLaplacian[i][j][k] += (Hubble_prime-Hubble**2)*pidata[i][j][k]

    integral = np.real(np.fft.ifftn(FTinverseLaplacian))

    del integrand, FT, FTinverseLaplacian

    piprimedatanext = np.zeros((meshpoints,meshpoints,meshpoints))

    for i in range(0, meshpoints):
        for j in range(0, meshpoints):
            for k in range(0, meshpoints):
                piprimedatanext[i][j][k] += -(Hubble/2)*(gradsquared[i][j][k])*(timestep)
                piprimedatanext[i][j][k] += -(gradpipiprime[i][j][k])*timestep
                piprimedatanext[i][j][k] += - (integral[i][j][k])*timestep

    if sigma > 0:
        k_cutoff = 6.28318/sigma
        piprimedatanext_FT = np.fft.fftn(piprimedatanext)
        
        for i in range(0, meshpoints):
            for j in range(0, meshpoints):
                for k in range(0, meshpoints):
                    k_mag = (6.28318/boxsize)*np.sqrt(min(i, meshpoints-i)**2 + min(j,meshpoints-j)**2 + min(k, meshpoints-k)**2 + 1E-24)
                    piprimedatanext_FT[i][j][k] *= np.exp(- (k_mag**2)/(2*k_cutoff**2))
                    
        piprimedatanext = np.real(np.fft.ifftn(piprimedatanext_FT))
        del piprimedatanext_FT

    for i in range(0, meshpoints):
        for j in range(0, meshpoints):
            for k in range(0, meshpoints):
                pidata[i][j][k] += (piprimedata[i][j][k])*timestep
                piprimedatanext[i][j][k] += (-(2*Hubble)*piprimedata[i][j][k] - (2*Hubble_prime)*pidata[i][j][k])*timestep
                piprimedatanext[i][j][k] += (tildec2_s*theta[i][j][k])*timestep
                piprimedatanext[i][j][k] += (tildec2_v_overH*thetaprime[i][j][k])*timestep
                piprimedata[i][j][k] += piprimedatanext[i][j][k]

    time = time + timestep
    print("current scale factor is " + str(a_scale))

    lin_pi_growth_times_a2 += a_step/((H0/299792.458)*np.sqrt(((Omega_m/(a_scale**3))+(1-Omega_m))))
    a_scale += a_step
    numsteps += 1

prog_end = myTime.time()

a_final = a_scale
print("Running time for ", numsteps, " steps is:",(prog_end-prog_start)*(10**3),"ms.\n" )
print("Final time is:", time, " Mpc/c and final scale factor is: ", a_scale, "\n")

####  4. Write the data to text files.

file1 = open(pion_output_path,"w")
file1.write("[")
for i in range(0, meshpoints-1):
    file1.write("[")
    for j in range(0, meshpoints-1):
        file1.write("[")
        for k in range(0, meshpoints-1):
            file1.write(str(pidata[i][j][k])+",")
        file1.write(str(pidata[i][j][meshpoints-1]))
        file1.write("],")
    file1.write("[")
    for k in range(0, meshpoints-1):
        file1.write(str(pidata[i][meshpoints-1][k])+",")
    file1.write(str(pidata[i][meshpoints-1][meshpoints-1]))
    file1.write("],")
file1.write("[")
for j in range(0, meshpoints-1):
    file1.write("[")
    for k in range(0, meshpoints-1):
        file1.write(str(pidata[meshpoints-1][j][k])+",")
    file1.write(str(pidata[meshpoints-1][j][meshpoints-1]))
    file1.write("],")
file1.write("[")
for k in range(0, meshpoints-1):
    file1.write(str(pidata[meshpoints-1][meshpoints-1][k])+",")
file1.write(str(pidata[meshpoints-1][meshpoints-1][meshpoints-1]))
file1.write("]]]")
file1.close()

file2 = open(overdensity_output_path,"w")
file2.write("[")
for i in range(0, meshpoints-1):
    file2.write("[")
    for j in range(0, meshpoints-1):
        file2.write("[")
        for k in range(0, meshpoints-1):
            file2.write(str(overdensity[i][j][k])+",")
        file2.write(str(overdensity[i][j][meshpoints-1]))
        file2.write("],")
    file2.write("[")
    for k in range(0, meshpoints-1):
        file2.write(str(overdensity[i][meshpoints-1][k])+",")
    file2.write(str(overdensity[i][meshpoints-1][meshpoints-1]))
    file2.write("],")
file2.write("[")
for j in range(0, meshpoints-1):
    file2.write("[")
    for k in range(0, meshpoints-1):
        file2.write(str(overdensity[meshpoints-1][j][k])+",")
    file2.write(str(overdensity[meshpoints-1][j][meshpoints-1]))
    file2.write("],")
file2.write("[")
for k in range(0, meshpoints-1):
    file2.write(str(overdensity[meshpoints-1][meshpoints-1][k])+",")
file2.write(str(overdensity[meshpoints-1][meshpoints-1][meshpoints-1]))
file2.write("]]]")
file2.close()

file3 = open(pion_prime_output_path,"w")
file3.write("[")
for i in range(0, meshpoints-1):
    file3.write("[")
    for j in range(0, meshpoints-1):
        file3.write("[")
        for k in range(0, meshpoints-1):
            file3.write(str(piprimedata[i][j][k])+",")
        file3.write(str(piprimedata[i][j][meshpoints-1]))
        file3.write("],")
    file3.write("[")
    for k in range(0, meshpoints-1):
        file3.write(str(piprimedata[i][meshpoints-1][k])+",")
    file3.write(str(piprimedata[i][meshpoints-1][meshpoints-1]))
    file3.write("],")
file3.write("[")
for j in range(0, meshpoints-1):
    file3.write("[")
    for k in range(0, meshpoints-1):
        file3.write(str(piprimedata[meshpoints-1][j][k])+",")
    file3.write(str(piprimedata[meshpoints-1][j][meshpoints-1]))
    file3.write("],")
file3.write("[")
for k in range(0, meshpoints-1):
    file3.write(str(piprimedata[meshpoints-1][meshpoints-1][k])+",")
file3.write(str(piprimedata[meshpoints-1][meshpoints-1][meshpoints-1]))
file3.write("]]]")
file3.close()

print("Output files written: " + str(pion_output_path) + " " + str(pion_prime_output_path) + " " + str(overdensity_output_path))

if show_final_configuration == 1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cax = ax.scatter(mesh[:, :, :, 0].flatten(), mesh[:, :, :, 1].flatten(), mesh[:, :, :, 2].flatten(), c=pidata.flatten(), cmap='jet', alpha=1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.suptitle('Pion field in 3D')
    plt.title('Conformal time = ' + str(time) + ', a = ' + str(a_final))
    color_bar = plt.colorbar(cax)
    color_bar.set_alpha(1)
    fig.draw_without_rendering()
    plt.show()
