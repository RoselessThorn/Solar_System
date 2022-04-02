"""
 CompMod Ex3: Particle3D, a class to describe point particles in 3D space

 An instance describes a particle in Euclidean 3D space: 
 velocity and position are [3] arrays

 Includes time integrator methods +...

author: s1969574

"""
import math
import numpy as np


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the particle
    mass: mass of the particle
    pos: position of the particle
    vel: velocity of the particle

        Methods:
    __init__
    __str__
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_p3d - initializes a P3D instance from a file handle
    sys_kinetic - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    """

import math
import numpy as np


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the particle
    mass: mass of the particle
    pos: position of the particle
    vel: velocity of the particle

        Methods:
    __init__
    __str__
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_p3d - initializes a P3D instance from a file handle
    sys_kinetic - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    """

    def __init__(self, label, mass, position, velocity):
        """
        Initialises a particle in 3D space

        :param label: String w/ the name of the particle
        :param mass: float, mass of the particle
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """
        self.label = label
        self.mass = mass
        self.pos = np.array([position[0],position[1],position[2]])
        self.vel = np.array([velocity[0],velocity[1],velocity[2]])


    def __str__(self):
        """
        XYZ-compliant string. The format is
        <label>    <x>  <y>  <z>
        Convert them into strings in suitable format and return
        """
        reStr = str(self.label) + '    ' + str(self.pos[0]) + '   ' + str(self.pos[1]) + '   ' + str(self.pos[2]) + "\n"
        return (reStr)


    def kinetic_e(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2 m v**2
        where v is the overall velocity in 3 directions
        """
        
        return 0.5*self.mass*(np.linalg.norm(self.vel)**2)


    def momentum(self):
        """
        Returns the momentum of a Particle3D instance
        
        :return p:float, mv
        where v is the overall velocity in 3 directions
        """

        return self.mass*self.vel
    

    def update_pos(self, dt):
        """
        To update the position in each direction of the particle after a time dt by velocity
        """
        self.pos += dt * self.vel
        


    def update_pos_2nd(self, dt, force):
        """
        To update the position in each direction of the particle after a time dt
        by velocity and force
        """

        self.pos += dt*self.vel + 0.5*(dt**2)*force/self.mass


    def update_vel(self, dt, force):
        """
        Update the velocity in each direction of the particle caused by a force
        """
        self.vel += dt*force/self.mass

        
        
    @staticmethod
    def new_particle(file_handle):
        """
        Initialises a Particle3D instance given an input file handle.
        
        The input file should contain one line per planet in the following format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>
        
        Note: I changed the name for this method to suit the one used in the checker
        :param inputFile: Readable file handle in the above format

        :return Particle3D instance
        """
        #Read the file_handle line by line
        f = file_handle.readline()
        #Split the line by space
        f = f.split()
        #Convert elements to float apart from the first element which is the label string
        for i in range(1,len(f)):
            f[i] = float(f[i])
        return Particle3D(f[0],f[1],[f[2],f[3],f[4]],[f[5],f[6],f[7]])


    @staticmethod
    def sys_kinetic(p3d_list):
        """
        Computes the total kinetic energy of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return sys_ke: The total kinetic energy of the system 
        """
        i = 0
        #Set the inital total kinetic energy to 0 then add the energy one by one in the list
        sys_ke = 0
        for i in range(len(p3d_list)):
            sys_ke = sys_ke + p3d_list[i].kinetic_e()
        return sys_ke


    @staticmethod
    def com_velocity(p3d_list):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return total_mass: The total mass of the system 
        :return com_vel: Centre-of-mass velocity
        """
        i = 0
        #Set the inital total mass and momentum to 0
        total_mass = 0
        total_momentum = 0
        #Add mass and momentum of each particle one by one from the list
        for i in range(len(p3d_list)):
            total_mass = total_mass + p3d_list[i].mass
            total_momentum = total_momentum + p3d_list[i].momentum()
        #Calculate the center of mass velocity
        com_vel = total_momentum / total_mass
        return total_mass, com_vel