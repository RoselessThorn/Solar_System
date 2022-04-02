# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 19:56:09 2021

@author: s1969574
"""

"""
CMod Priject A: Solar System
N-body system interacting through Newtonian gravity.

Produces plots of the position of the Sun and other objects, both as function of time. Also
saves to file. The Sun, the moon, the Earth and Halley's comet

The potential is Gravitational potential, where
all parameters required are hard-coded in the main() method
and passed to the functions that
calculate force and potential energy.
"""

import sys
import math 
import numpy as np
import matplotlib.pyplot as pyplot
from particle3D import Particle3D as p3d
from scipy.signal import argrelextrema

def force_gravity(particle1, particle2, G):
    """
    Method to return the force on two particles at r1 and r2
    in a Morse potential.
    Force is given by
    F(r1,r2) = -G*m1*m2*(r_12/(rm_12)**3)
    G is a constant in unit of AU^3Mass^-1Time^-2
    Mass has unit Earth's mass, time has unit Earth's day
    r_12 = r1 - r2 as the direction vector array
    rm_12 is the length of r_12
    rh_12 is the unit vector array
    
    :param particle1,2: Particle3D instance
    :param r1: particle position r1 in units AU
    :param r2: particle position r2
    :param m1: particle mass m1 in units of Earth's mass
    :param m2: particle mass m2
    :return: force acting on particle as a 3D Numpy array in all directions
    """
    r1 = particle1.pos
    r2 = particle2.pos
    m1 = particle1.mass
    m2 = particle2.mass
    r_12 = r2 - r1
    rm_12 = np.linalg.norm(r_12)
    force = -G*m1*m2*(r_12/((rm_12)**3))
    return force


def pot_energy_gravity(particle1, particle2, G):
    """
    Method to return potential energy in eV
    of particle in Morse potential
    U(r1,r2) = D_e*{[1-exp[-α(r_12-r_e)]]^2-1}

    G is a constant in unit of AU^3Mass^-1Time^-2
    r_12 = r1 - r2
    rm_12 is the length of r_12
    rh_12 is the unit vector array
    :param particle1,2: Particle3D instance
    :param r1: particle position r1 in units of AU
    :param r2: particle position r2
    :param m1: particle mass m1 in units of Earth's mass
    :param m2: particle mass m2
    :return: potential energy of particle as float
    """
    r1 = particle1.pos
    r2 = particle2.pos
    m1 = particle1.mass
    m2 = particle2.mass
    r_12 = r1 - r2
    rm_12 = np.linalg.norm(r_12)
    energy = -G*m1*m2/rm_12
    return energy


# Define a function to find periods
def periods(pos_list,dt):
    # Find the peaks of the position-time diagram, make sure the position index matches the time index
    max1 = argrelextrema(np.array(pos_list),np.greater)
    max1 = [element * dt for element in max1]
    # Period is the average difference between maximum peaks
    # If peaks list is empty or with just 1 peak(the object hasn't finish a whole circle),
    # set period to NaN
    if max1[0].size == 0:
        period = 'NaN'
    elif max1[0].size == 1:
        period = 'NaN'
    # Write periods in 2 decimal places
    else:    
        period = round(np.average(np.diff(max1[0])),2)
    return period 

# Define a function to find apsides
def apsides(distance,numstep,dt):
    aphelions = 0
    perihelions = 0
    # Check whether the object has a proper period
    if (periods(distance,dt)) == 'NaN': 
        aphelions = 'NaN'
        perihelions = 'NaN'
    # Apsis are the max and min position of the position list within a period
    else:
        period = int(periods(distance,dt))
        aphelions = (max(distance[:period]))
        perihelions = (min(distance[:period]))
    
    return aphelions, perihelions

# Define a function to calculate energy deviation, which is the difference between the 
# largest and smallest energy values
def energy_deviation(energy_list):
    energy_deviation = max(energy_list) - min(energy_list)
    
    
    return energy_deviation



# Begin main code
def main():
    # Input names of output files from command line
    outfile_name_1 = input('XYZ Outputfile name: ')
    outfile_name_2 = input('Energy Outputfile name: ')
    
    
    # Open output file
    outfile_1 = open(outfile_name_1, "w")
    outfile_2 = open(outfile_name_2, "w")
    
    
    if outfile_name_1 == outfile_name_2:
        print("Invalid input, please input different names for two files")
        quit()
    # Set an empty list for particles
    particles = []
    
    
    # Read the parameters
    time = 0.
    G = 8.8877e-10
    
    # Input dt and number step, dt should be greater than 0 and numstep should be greater than 1
    dt = float(input('Please input the time between each step '))
    numstep = int(float(input('Please input the total simulation steps ')))
    if dt <= 0:
        print("Invalid input, dt should be greater than 0.")
        quit()
    if numstep <= 1:
        print("Invalid input, numstep should be greater than 1.")
        quit()
        
        
    # Try to read the line from initial condition file until it breaks,
    # append lines to the particle list
    with open(input('Please input the initial condition(type:Solar): ')+'.dat') as a:
        while True:
            try:
                p = p3d.new_particle(a)
                particles.append(p)
            except:
                break
        
    # Record the total number of particles, very important and most used later on
    n = len(particles) 
  
    
    # Find the indices for the Earth and the Moon
    for i in range(n):
        if particles[i].label == 'Earth':
            e = i
        if particles[i].label == 'Moon':
            m = i
            
            
    energy = 0
    # Write out initial energy and earth-sun seperation
    i = 0
    j = 0
    p_energy_total = 0
    k_energy_total = 0
    for i in range(n):
        k_energy = particles[i].kinetic_e()
        k_energy_total += k_energy
        energy += k_energy
        for j in range(i+1,n):
            p_energy = pot_energy_gravity(particles[i], particles[j], G)
            p_energy_total += p_energy
            energy += p_energy
 
    
    # Write to output files
    outfile_1.write("{0:}\nPoint = {1}\n".format(n,time))
    for i in range(n):
        outfile_1.write("{0:}".format(particles[i]))
    outfile_2.write("{0:f}{1:12.8f}{2:12.8f}{3:12.8f}\n".format(time,energy,k_energy_total,p_energy_total))


    # Set up loops to read positions and masses from particles' list
    # Then loop to find seperations and corresponding modules
    r = []
    mass = []
    separation= np.zeros((n,n,3))
    modules =np.zeros((n,n))

    
    for i in range (n):
        r.append(particles[i].pos)
        mass.append(particles[i].mass)


    for i in range (n):
        for j in range(n):
            if i != j:
                separation[i][j]= r[j] - r[i]
            if i == j:
                separation[i][j] = 0

    
    for i in range (n):
        for j in range (n):
            modules[i][j] = (np.linalg.norm(separation[i][j]))
    
    
    # Set a matrix to store the distance form the Sun and Earth-moon distance
    S_distance = np.split(modules[:,0],n)
    E_distance = np.zeros((numstep,1))
    E_distance[0] = modules[e][m]


    # Set an empty martix to store modules and force
    force = np.zeros((n,n,3))


    # Calculate initial forces
    for i in range(n):
        for j in range(n):
            if i != j:              
                force[i][j] = -G*mass[i]*mass[j]*((separation[i][j])/(modules[i][j])**3)     
            if i == j:
                force[i][j] = 0
    

    # Initialise data lists for plotting later
    # One specical position list for earth-moon disance
    time_list = []
    force_new = np.zeros((n,n,3))
    energy_list = np.zeros((numstep,1))
    pos_list = np.zeros((n,numstep,1))
    pos_moon_list = np.zeros((numstep,1))
    for i in range(n):
        pos_list[i][0] = r[i][0]
    pos_moon_list[0] = r[m][0] - r[e][0]


    # Center of mass correction
    momentum = 0
    total_mass = 0
    for i in range(n):
        momentum = momentum + particles[i].mass * particles[i].vel
        total_mass = total_mass + particles[i].mass
    com_velocity = momentum/total_mass


    for i in range(n):
        particles[i].vel = particles[i].vel - com_velocity


    # Start the time integration loop, since intial condition has been recorded, time starts from 1dt
    k = dt
    
    
    # Start loop
    for k in range(numstep):
    # Update object position 
        for i in range(n):
            particles[i].update_pos_2nd(dt,-force.sum(axis=1)[i])

            
    # Update the new seperation and corresponding modulus, record new distances
        for i in range(n):
            for j in range(n):
                if i != j:
                    separation[i][j] = (r[j] - r[i])
                    modules[i][j] = np.linalg.norm(separation[i][j])
                    # Record new distances to the Sun and Earth-Moon distance
                    if i == 0:
                        S_distance[j] = np.append(S_distance[j],modules[i][j])
                    if i == 3:
                        E_distance[k] = modules[e][m]
                if i == j:
                    separation[i][j] = 0
                    modules[i][j] = 0


    # Calculate the updated force
        for i in range(n):
            for j in range(n):
                if i != j:
                    force_new[i][j] = -G*mass[i]*mass[j]*((separation[i][j])/((modules[i][j])**3))
                if i == j:
                    force_new[i][j] = 0
                    
    
    # Update velocity and force
        for i in range(n):
            particles[i].update_vel(dt,0.5*(-force.sum(axis=1)[i]-force_new.sum(axis=1)[i]))              
            force[i] = force_new[i]


    # Record new total energy, kinetic energy and potential energy   
        energy = 0     
        p_energy_total = 0
        k_energy_total = 0
        for i in range(n):
            k_energy = particles[i].kinetic_e()
            k_energy_total += k_energy
            energy += k_energy
            for j in range(i+1,n):
                p_energy = pot_energy_gravity(particles[i], particles[j], G)
                p_energy_total += p_energy
                energy += p_energy
                energy_list[k] = energy


    # Append all information needed to data lists
        time_list.append(time)
        time += dt
        i = 0
        try:
            for i in range(n):
                pos_list[i][k+1] = r[i][0]
        except:
            break
        try:
            pos_moon_list[k+1] = r[m][0] - r[e][0]
        except:
            break


    # Output particle information, write particle positions and time into the first output file,
    # write energies into the second output file
        outfile_1.write("{0:}\nPoint = {1}\n".format(n,time))
        for i in range(n):
            outfile_1.write("{0:}".format(particles[i]))
      
        
        outfile_2.write("{0:f}{1:12.8f}{2:12.8f}{3:12.8f}\n".format(time,energy,k_energy_total,p_energy_total))

    
    # Print periods of objects orbiting aroung the Sun and print the Moon as one special case
    for i in range (1,n):
        if i == m:
            continue
        else:
            print('Period of',str(particles[i]).split()[0],'is',periods(pos_list[i],dt),'days')
    print('Period of',str(particles[m]).split()[0],'is',periods(pos_moon_list,dt),'days')
    
    
    # Print apsides of objects orbiting around the Sun, also print the Moon as a special case
    for i in range (1,n):
        if i == m:
            continue
        else:
            print('Apsides of',str(particles[i]).split()[0],'are',apsides(S_distance[i],numstep,dt)[0],'and',apsides(S_distance[i],numstep,dt)[1],'AU')
    print('Apsides of',str(particles[m]).split()[0],'are',float(apsides(E_distance,numstep,dt)[0]),'and',float(apsides(E_distance,numstep,dt)[1]),'AU')


    # Print the energy deviation of the whole system
    print('Energy deviation is',"%.3g"% float(energy_deviation(energy_list)),'AU^2 M_earth day^-2')
    
    
    # Post-simulation:
    # Close output file
    outfile_1.close()
    outfile_2.close()


    # Plot object trajectory to screen, change the pos_list index to change object, defalut is Mercury
    pyplot.title('Velocity Verlet: x-position of objects seperations vs time')
    pyplot.xlabel('Time/[Days]')
    pyplot.ylabel('x-position of objects seperationvs/AU')
    pyplot.plot(time_list, pos_list[1])
    pyplot.show()
    
    
    # Plot system energy to screen
    pyplot.title('Velocity Verlet: total energy vs time')
    pyplot.xlabel('Time/[Day]')
    pyplot.ylabel('Energy/[AU^2 M_earth day^-2]')
    pyplot.plot(time_list, energy_list)
    pyplot.show()


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()