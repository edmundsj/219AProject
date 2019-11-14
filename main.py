# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
# CONTEXT, APPLICATIONS
# This is a Fourier-Optics based physical optics simulator. Its main purpose is to calculate
# monochromatic and polychromatic wave propagation along an optical axis. Typical applications
# might include:
#   - diffraction through a slit
#   - Gaussian beam propagation through free space or an optical system
#   - Diffraction orders and patterns upon reflecting from a diffraction grating
#   - Spot size/shape of an optical focusing system
#   - Power coupling efficiency into waveguide
#   - Formation of speckle patterns through an inhomogenous scattering medium
#   - Anything that can be adequately described by scalar wave diffraction theory
#   
# Why did I build this simulator? Because I think optics is awesome and Fourier analysis is
# the literally the blood that flows in my veins, I needed an optical simulator for speckle pattern
# formation and diffraction from a transparent cylinder, and I was taking a computational class
# where I had the freedom to actually build said simulator. It has been useful to me, hopefully it will
# be useful to others.
# 
# SIMULATION STRUCTURE
# This is a netlist-based simulator. That means that the only interaction a user needs to have
# with the software is to feed in a netlist into the main program, and run it. For example, you might run
# using the code:
#   $ python3 main.py "my_netlist.txt"
# Everything is defined in the netlist. All your optical components, their location, your excitation source,
# and all the planes of interest at which you want to know the electromagnetic field.
# Why use a netlist? Because it's a pain in the a$$ to read other people's code. I don't like doing it,
# and I have never met anyone who likes doing it. Part of the success of electronic simulators (SPICE)
# also developed at Berkeley, is that they are netlist-based. They require no knowledge of the fundamentals
# of circuit simulation or numerical analysis, and so they dramatically lower the barrier to entry for
# people to start simulating circuits. This was my main goal in making this simulator. I wanted to makke
# the barrier to entry for computational optics as low for other people (and my future self) as possible.
# For that reason, I believe that complexity in actually using the code should only be made apparent on
# an as-needed basis. The simplest simulation you can run, then, should require an extremely minimal
# set of things. Indeed, for the simplest simulation you can run, your netlist only needs three lines:
# one line describing your source (its location and type), one line describing a plane of interest
# for which you want to calculate the new electromagnetic field of that source, and one line describing the
# type of simulation you want to do (the simplest simulation is just a monochromatic simulation at
# a particular wavelength. Such a sample netlist for a Gaussian beam is below:
#   G1 0mm 500um
#   P1 20mm
#   .mono 700nm

# This runs a simulation of a Gaussian beam with a beam waist of 500 microns, and propagates it 20 millimeters.
# The wavelength of the simulation is set to 700nm. The results (the outputs) will be saved as numpy arrays:
# the magnitude and phase of each plane of interest, all sources, and all masks. You can define (and in
# fact it's easier if you do define) analytic amplitude/phase functions, and these will be sampled during
# the simulation appropriately. You can also feed in image files as an amplitude/phase mask.

# INTERNAL STRUCTURE
# For those of you interested in the inner workings of this simulator, it is split up as below:
# -1. You write your netlist
# 0. You run "main.py" using python3 and pass in as an argument your netlist.
# 1. The class NetlistParser parses a netlist and turns everything into a "Mask" object. 
#   These objects just contain the complex amplitude value you are applying to the 
#   electromagnetic field as it propagates through that plane. Masks are fundamentally 
#   two-dimensional. Any three-dimensional objects are sliced up into a bunch of 2D masks 
#   separated by some distance (which the simulator handles). These masks can be defined as sources or
#   'passive' elements, (lenses apertures, stops, etc.). You may have multiple sources at various locations
#   in your simulation. The results will be saved separately and optionally summed.

# 2. A 'Field' object is created for each source the netlist parser finds.
#   This 'Field' object takes the information about the source and creates a discretized version - it 
#   samples the electromagnetic field of the source appropriately (still have to figure out exactly
#   how to do this)
# 3. That 'Field' object is then propagated to each plane and Mask
#   The 'Field' object is then propagated (using Kirchoff/Fresnel/Fraunhoffer diffraction
#   as appropriate to each plane and each mask as they appear along the optical axis. If you first
#   have a plane at z=10mm, and your source is at z=0mm, then the Field gets propagated using
#   the appropriate diffraction formula (and probably a few FFTs). The complex-valued field strength
#   is then saved to a numpy array and written to file at each plane.
# 3a. When appropriate, masks are applied to the field.
#   If you defined any masks (lenses, apertures, and so on), the Field will get propagated to the
#   location of that mask, and then the mask will get applied to the field. The Field will then
#   be saved (both before and after the mask), and then propagated along its merry way.
# 4. When the Field reaches the final plane or mask, the simulation is terminated.
#   Everything is then written to files (which will just be named according to how you named them in
#   the netlist), and the final output field is plotted for you to view. If you would like to use
#   the utilities I have written to load in the results of your simulation and plot them, be my guest.

from netlist.netlist_parser import NetlistParser

parser = NetlistParser('./netlist/sample_netlist.txt');
my_objects = parser.parseNetlist();
print(my_objects);
