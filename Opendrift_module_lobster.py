# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway
#get_environment(reader,['sea_water_temperature'],startTime,lon,lat,0*lat,None)

#
# CODE to INCLUDE N. norvegicus ecological information
#
# Written by Morane Clavel-Henry in the context of BITER project
#

import os, sys
import numpy as np
from datetime import datetime, timedelta
from astral import LocationInfo # Added
from astral.sun import sun #Added
import pylake

# Working path:
sys.path.append('path_to_folder')

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
import Vertical_module as vp # to call vertical and horizontal diffusivity
from opendrift.elements import LagrangianArray
from opendrift.readers import reader_ROMS_native

# Defining the ecological properties of lobster larvae
class Pelagiclobster(Lagrangian3DArray):
    """ Extending Lagrangian3DArray with specific properties for lobster larval behaviour
    """

    variables = Lagrangian3DArray.add_variables([
        ('age_seconds', {'dtype': np.float32,
                         'units': 's',
                         'default': 0.}),
        ('tau', {'dtype': np.float32,  # See Phelps et al., 10.3354/meps11040
                         'units': 's',
                         'default': 0.}),
        ('particles_velocity',{'dtype':np.float32, # Swimming velocity according to minimum velocity of Chu
                        'units':'',
                        'default':0.}),
        ('stages', {'dtype': np.int16, # Stages of the larvae (e.g.:0 = Z1, 1 = Z2, 2 = Z3, 3 = Decapodid 1, 4 = Decapodid +)
                    'units': '[]',
                    'default': 0})])
                    
	# Function to compute the duration of the stage according to Water temperature. Unit of stage duration is initially "day"
    def updateStageTau(self):
		dt = self.time_step.total_seconds()
        Temperature=self.environment.sea_water_temperature[np.all([self.elements.stages < 3],0)]
        A=np.array([-0.161, -0.175, -0.113])
        B=np.array([4.265, 4.646, 4.188])
        stages = np.int16(self.elements.stages[np.all([self.elements.stages < 3],0)])
        self.elements.tau[np.all([self.elements.stages < 3],0)] += dt/((np.exp(A[stages] * Temperature +B[stages]))*86400)
        
	# Function to change the stages of the larvae if tau >= 1
    def updateIDStages(self):
        self.elements.stages[np.all([self.elements.tau >= 1.],0)] += 1
        self.elements.tau[np.all([self.elements.tau >= 1.],0)] = 0. # return to 0 once the stage has been updated
        self.elements.stages[np.all([self.elements.stages >= 4],0)] = 4
        
    def Thermocline_depth(self, math_function):
        # compute the gradient (diff(T)/diff(z))
        #env,prof,miss = self.get_environment(list(self.required_variables),self.time,self.elements.lon,self.elements.lat,self.elements.z,self.required_profiles)
        prof = self.environment_profiles
        
        if math_function is 'diff':
            l= np.divide(np.diff(prof['sea_water_temperature'],axis=0).T,np.diff(prof['z']))
            # Find the index of the maximum and extract the depth
            depth_cline = np.take(prof['z'], np.argmax(l,axis=1)+1)
        if math_function is 'gradient':
            def my_function(tab):
                return(np.argmax(np.abs(tab)))
            repeated_array = np.tile(prof['z'][:, np.newaxis], (1, prof['sea_water_temperature'].shape[1]))
            Gradient=np.divide(np.gradient(prof['sea_water_temperature'],axis=0),np.gradient(repeated_array,axis=0))
            depth_cline = prof['z'][np.apply_along_axis(my_function, axis=0, arr=Gradient)]
            
        return depth_cline

class PelagicLobsterDrift(OceanDrift, Pelagiclobster):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at ICM-CSIC (SPAIN), from the original script of PelagicEgg developped at MET, Norway

        Generic module for particles that are subject to vertical or horizontal turbulent
        mixing with the possibility for positive or negative buoyancy
        Based on larval development of lobsters

        Under construction.
    """

    ElementType = Pelagiclobster

    required_variables = {'x_sea_water_velocity': {'fallback': 0},
                          'y_sea_water_velocity': {'fallback': 0},
                          'sea_surface_height': {'fallback': 0},
                          'land_binary_mask': {'fallback': None},
                          'sea_floor_depth_below_sea_level': {'fallback': -10},
                          'ocean_vertical_diffusivity': {'fallback': 0.0000001,'profiles':True},
                          'sea_water_temperature': {'fallback': 12,'profiles':True},
                          'sea_water_salinity': {'fallback': 37,'profiles':True},
                          'surface_downward_x_stress': {'fallback': 0},
                          'surface_downward_y_stress': {'fallback': 0},
                          'turbulent_kinetic_energy': {'fallback': 0},
                          'turbulent_generic_length_scale': {'fallback': 2000},
                          'upward_sea_water_velocity': {'fallback': 0},
                          'ocean_horizontal_diffusivity': {'fallback': 10},
                          'surface_boundary_layer': {'fallback': -50}
                          }

    # Vertical profiles of the following parameters will be available in
    # dictionary self.environment.vertical_profiles
    # E.g. self.environment_profiles['x_sea_water_velocity']
    # will be an array of size [vertical_levels, num_elements]
    # The vertical levels are available as
    # self.environment_profiles['z'] or
    # self.environment_profiles['sigma'] (not yet implemented)
    # The depth range (in m) which profiles shall cover
        
    def Vertical_Swimming(self, cline_depth): # check up which direction self.element.z is saved
        dt = self.time_step.total_seconds()
        # vertical swimming for Nephrops norvegicus
        Vel_Z1 = 6.0 # mm/s
        Vel_Z2 = 9.0
        Vel_Z3 = 9.0
        Vel_Z4 = 3.0
        # Define zLower according to day/night criteria. zLower = thermocline_depth - 10 m at night, zLower = - 40 m at day. Zupper = thermocline_depth
        # zUpper = -5.0
        zUpper = cline_depth
        print("Thermocline depth is %s "%zUpper[:])
        
        # Look at sunrise/sunset time based on one location
        loc_larvae = LocationInfo('Barcelona','Spain','Europe/Madrid')
        sun_time = sun(loc_larvae.observer,date=self.time,tzinfo=loc_larvae.timezone)
        dawn_time=sun_time['dawn'].hour
        dusk_time= sun_time['dusk'].hour
        
        if self.time.hour <= dusk_time and self.time.hour >= dawn_time :# Day
            zLower = cline_depth - 40.
        else: # Night
            zLower = cline_depth - 10.
            
        zBottom = -1.*self.environment.sea_floor_depth_below_sea_level
        
        # Pick a Random value "Velocity_Vel_Tamper" between 0 and 1 to tamper velocity of swimming
        Velocity_Vel_Tamper = np.random.uniform(0,1,len(self.elements.lat))
        # Pick a random value "Dir_rand" between 0 and 1 to define direction of the swimming
        Dir_rand = np.random.uniform(0,1,len(self.elements.lat))
        # Create an array of 1 with a size N_Particles
        Direction_z = np.ones(len(self.elements.lat))
        
        #------------ Define the direction of the swimming
        # 1) If being under zLower, 100% of Z1 swim up (no changes) and  10% of Z2 swim down
        Direction_z[np.all([self.elements.z < zLower,self.elements.stages == 1,Dir_rand <= 0.1],0)] = -1.
        # 2) Relatively to being above zUpper, 90% of Z1 and Z2 swim down
        Direction_z[np.all([self.elements.z > zUpper, self.elements.stages < 2, Dir_rand > 0.1],0)]= -1.
        # 3) Relatively to being between zLower and zUpper, 50% of Z1 and Z2 don't move, 25% swim down, 25% swim up
        Direction_z[np.all([self.elements.z > zLower,self.elements.z < zUpper,self.elements.stages < 2,Dir_rand > 0.25],0)] = 0.
        Direction_z[np.all([self.elements.z > zLower,self.elements.z < zUpper,self.elements.stages < 2,Dir_rand > 0.75],0)] = -1.
        # 4) Relatively to being more than 1 m above ZBottom, 90% of Z3 and Z4 swim down
        Direction_z[np.all([self.elements.z > zBottom + 1.0,self.elements.stages > 1,Dir_rand > 0.1],0)] = -1.
        # 5) Relatively to being below 1m above ZBottom, 90% of Z3  don't move, 5% swim down
        Direction_z[np.all([self.elements.z < zBottom + 1.0,self.elements.stages > 1,Dir_rand > 0.05],0)] = 0.
        Direction_z[np.all([self.elements.z < zBottom + 1.0,self.elements.stages > 1,Dir_rand > 0.95],0)] = -1.

       # ------------- Allocate the velocity of the swimming according to the Stage
        self.elements.particles_velocity[np.all([self.elements.stages == 0],0)] = Vel_Z1 * Velocity_Vel_Tamper[np.all([self.elements.stages == 0],0)]* Direction_z[np.all([self.elements.stages == 0],0)]
        self.elements.particles_velocity[np.all([self.elements.stages == 1],0)] = Vel_Z2 * Velocity_Vel_Tamper[np.all([self.elements.stages == 1],0)]*Direction_z[np.all([self.elements.stages == 1],0)]
        self.elements.particles_velocity[np.all([self.elements.stages == 2],0)] = Vel_Z3 * Velocity_Vel_Tamper[np.all([self.elements.stages == 2],0)]*Direction_z[np.all([self.elements.stages == 2],0)]
        self.elements.particles_velocity[np.all([self.elements.stages >= 3],0)] = Vel_Z4 * Velocity_Vel_Tamper[np.all([self.elements.stages >= 3],0)]*Direction_z[np.all([self.elements.stages >= 3],0)]
        
        self.elements.z = self.elements.z + self.elements.particles_velocity*dt*0.001# 0.001: Conversion to m/s.
        
        bottom = np.where(self.elements.z < zBottom)
        if len(bottom[0]) > 0:
            print('%s elements reached seafloor, set to bottom' % len(bottom[0]))
            self.interact_with_seafloor()
            self.bottom_interaction(zBottom)
            
        self.surface_stick()
        
        
    def __init__(self, *args, **kwargs):
        # Calling general constructor of parent class
        super(PelagicLobsterDrift, self).__init__(*args, **kwargs)
          # Configuration
        configspec = '''
            [drift]
                scheme = string(default='runge-kutta')
                vertical_advection = boolean(default=True)
                vertical_mixing = boolean(default=False)
                horizontal_mixing = float(default=0.0)
                vertical_swimming = boolean(default=True)
            [vertical_mixing]
                timestep = float(min=0.1, max=3600, default=1.)
                verticalresolution = float(min=0.01, max=10, default = 1.)
                diffusivitymodel = string(default='environment')
            '''

        self._add_config({'drift:vertical_swimming': {
            'type':'bool',
            'default': False,
            'level': self.CONFIG_LEVEL_ADVANCED,
            'description':'Vertical swimming of the larvae'}
            })
        self.set_config('drift:vertical_swimming',False)



    def update(self):
        """Update positions and properties of particles."""
        #print("Temp Profile 1 is %s "%self.environment_profiles['sea_water_temperature'])
        # Update element age
        # self.elements.age_seconds += self.time_step.total_seconds()

        # Larval development
        self.updateStageTau()
        self.updateIDStages()

        if self.get_config('drift:vertical_mixing') is True: 
            self.vertical_mixing()
            
        # Disable vertical swimming of PL++
        NOW=self.elements.stages
        if np.mean(NOW) > 3.:
            for ind in xrange(len(self.elements.lat)):
                if NOW[ind]==4.:
                    self.elements.particles_velocity[ind] = 0. # m/s


        #Vertical swimming
        if self.get_config('drift:vertical_swimming') is True:
            cline_depth = self.Thermocline_depth(math_function='gradient')
            cline_depth[cline_depth == 0.] = -5.
            cline_depth[cline_depth <= -200.] = -200.
            self.Vertical_Swimming(cline_depth) ### Swimming vertical velocity updated, if drift:vertical_swimming is True
            #print("The velocity of larvae is: %s"%(np.mean(self.elements.particles_velocity)*1000))
    
        # Horizontal advection
        #print "The velocity of current is: %s"%self.environment.x_sea_water_velocity
        #vp.horizontal_mixing(self,'okubo')
        #print "The velocity of current + turbulence is: %s"%self.environment.x_sea_water_velocity
        self.advect_ocean_current()
        
        # Vertical advection
        if self.get_config('drift:vertical_advection') is True:
            self.vertical_advection()
