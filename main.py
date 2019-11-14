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
# people to start simulating circuits. This was my main goal in making this simulator.



from netlist.netlist_parser import NetlistParser

parser = NetlistParser('./netlist/sample_netlist.txt');
my_objects = parser.parseNetlist();
print(my_objects);
