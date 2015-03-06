# -- start license --
# Highly Optimized Object-oriented Many-particle Dynamics -- Blue Edition
# (HOOMD-blue) Open Source Software License Copyright 2009-2014 The Regents of
# the University of Michigan All rights reserved.

# HOOMD-blue may contain modifications ("Contributions") provided, and to which
# copyright is held, by various Contributors who have granted The Regents of the
# University of Michigan the right to modify and/or distribute such Contributions.

# You may redistribute, use, and create derivate works of HOOMD-blue, in source
# and binary forms, provided you abide by the following conditions:

# * Redistributions of source code must retain the above copyright notice, this
# list of conditions, and the following disclaimer both in the code and
# prominently in any materials provided with the distribution.

# * Redistributions in binary form must reproduce the above copyright notice, this
# list of conditions, and the following disclaimer in the documentation and/or
# other materials provided with the distribution.

# * All publications and presentations based on HOOMD-blue, including any reports
# or published results obtained, in whole or in part, with HOOMD-blue, will
# acknowledge its use according to the terms posted at the time of submission on:
# http://codeblue.umich.edu/hoomd-blue/citations.html

# * Any electronic documents citing HOOMD-Blue will link to the HOOMD-Blue website:
# http://codeblue.umich.edu/hoomd-blue/

# * Apart from the above required attributions, neither the name of the copyright
# holder nor the names of HOOMD-blue's contributors may be used to endorse or
# promote products derived from this software without specific prior written
# permission.

# Disclaimer

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS'' AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND/OR ANY
# WARRANTIES THAT THIS SOFTWARE IS FREE OF INFRINGEMENT ARE DISCLAIMED.

# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -- end license --

# this file exists to mark this directory as a python module

from hoomd_plugins.SurfaceCharge import _SurfaceCharge

# Next, since we are extending an updater, we need to bring in the base class updater and some other parts from 
# hoomd_script
from hoomd_script.analyze import _analyzer;
from hoomd_script import util
from hoomd_script import globals
import hoomd


# need to import all submodules defined in this directory

# NOTE: adjust the import statement to match the name of the template
# (here: plugin_template)


# class to initialize the dumpWriter
#
# Every \a period time steps, particles are checked for evaporation and thermostatting
#
            
class xml(_analyzer):
    def __init__(self, filename="dump", period=None, time_step=None, **params):
        util.print_status_line();

        # initialize base class
        _analyzer.__init__(self);

        # create the c++ mirror class
        self.cpp_analyzer = _SurfaceCharge.HOOMDTypeDumpWriter(globals.system_definition, filename);
        util._disable_status_lines = True;
        self.set_params(**params);
        util._disable_status_lines = False;

        if period is not None:
            self.setupAnalyzer(period);
            self.enabled = True;
            self.prev_period = 1;
        elif filename != "dump":
            util._disable_status_lines = True;
            self.write(filename, time_step);
            util._disable_status_lines = False;
        else:
            self.enabled = False;
    
    def excludeType(self, type):
        self.cpp_analyzer.excludeType(type);
    
    def excludeTypes(self, types):
        for type in types:
            self.cpp_analyzer.excludeType(type);
    
    def includeType(self, type):
        self.cpp_analyzer.includeType(type);
    
    def includeTypes(self, types):
        for type in types:
            self.cpp_analyzer.includeType(type);
        
    def preserveTypes(self, enable):
        self.cpp_analyzer.preserveTypes(enable);

    def set_params(self,
                   all=None,
                   vis=None,
                   position=None,
                   image=None,
                   velocity=None,
                   mass=None,
                   diameter=None,
                   type=None,
                   body=None,
                   wall=None,
                   bond=None,
                   angle=None,
                   dihedral=None,
                   improper=None,
                   acceleration=None,
                   charge=None,
                   orientation=None,
                   vizsigma=None,
                   included=None,
                   excluded=None,
                   preserveTypes=None):
        util.print_status_line();
        self.check_initialization();

        if all:
            position = image = velocity = mass = diameter = type = wall = bond = angle = dihedral = improper = True;
            acceleration = charge = body = orientation = True;

        if vis:
            position = mass = diameter = type = body = bond = angle = dihedral = improper = charge = True;

        if position is not None:
            self.cpp_analyzer.setOutputPosition(position);

        if image is not None:
            self.cpp_analyzer.setOutputImage(image);

        if velocity is not None:
            self.cpp_analyzer.setOutputVelocity(velocity);

        if mass is not None:
            self.cpp_analyzer.setOutputMass(mass);

        if diameter is not None:
            self.cpp_analyzer.setOutputDiameter(diameter);

        if type is not None:
            self.cpp_analyzer.setOutputType(type);

        if body is not None:
            self.cpp_analyzer.setOutputBody(body);

        if wall is not None:
            self.cpp_analyzer.setOutputWall(wall);

        if bond is not None:
            self.cpp_analyzer.setOutputBond(bond);

        if angle is not None:
            self.cpp_analyzer.setOutputAngle(angle);

        if dihedral is not None:
            self.cpp_analyzer.setOutputDihedral(dihedral);

        if improper is not None:
            self.cpp_analyzer.setOutputImproper(improper);

        if acceleration is not None:
            self.cpp_analyzer.setOutputAccel(acceleration);

        if charge is not None:
            self.cpp_analyzer.setOutputCharge(charge);

        if orientation is not None:
            self.cpp_analyzer.setOutputOrientation(orientation);

        if vizsigma is not None:
            self.cpp_analyzer.setVizSigma(vizsigma);
        
        if included is not None:
            if isinstance(included,list):
                self.includeTypes(included)
            else:
                self.includeType(included)
        
        if excluded is not None:
            if isinstance(excluded,list):
                self.excludeTypes(excluded)
            else:
                self.excludeType(excluded)
        
        if preserveTypes is not None:
            self.preserveTypes(preserveTypes)

    def write(self, filename, time_step = None):
        util.print_status_line();
        self.check_initialization();

        if time_step is None:
            time_step = globals.system.getCurrentTimeStep()

        self.cpp_analyzer.writeFile(filename, time_step);
