# This is the core module which contains our field propagation.
# Explanation of overall simulation:
#

### This is the class for essentially everything, it represents a two-dimensional 
# complex-valued mask for our electromagnetic field. This is the general structure used to 
# represent phase masks, sources, etc. Everything is in cartesian coordinates
# and all units are assumed to be millimeters.
# TODO:
# 1. Add Fresnel Number / Error calculation in fresnelPropagate
# 2. Add impulse-response based fresnel propagator
# 3. Add warnings if the source (the majority of the intensity) does not *fit* within
# the aperture (L > 2-3*Dsource). 

import numpy as np

class Mask:
    def __init__(self, ident, location, source, phaseFunction=None, intensityFunction=None,
                realFunction=None, imaginaryFunction=None,
                phaseValues=None, intensityValues=None,
                realValues=None, imaginaryValues=None):

        self.z = location;
        self.ident = ident;
        self.source = source;
        self.phaseFunction = phaseFunction;
        self.intensityFunction = intensityFunction;

        self.realFunction = realFunction;
        self.imaginaryFuunction = imaginaryFunction;
        self.phaseValues = phaseValues;
        self.intensityValues = intensityValues;
        self.realValues = realValues;
        self.imaginaryValues = imaginaryValues;

    source = False; # Is this field a source or not? (does it violate conservation of energy)
    ident = ''; # The unique identifier for this field.
    z = 0; # Location along the optical axis of our field
    phaseFunction = 0; # function describing the phase of our optical field
    intensityFunction = 1; # function describing the intensity of our optical field.
    realFunction = None;
    imaginaryFunction = None;
    phaseValues = None;
    intensityValues = None;
    realValues = None;
    imaginaryValues = None;

# This is just a container for our complex-valued field at a particular location.
class Field:
    def __init__(self, ident, loc, size_x, size_y, size_x_mm, size_y_mm):
        self.ident = ident;
        self.complexValues = np.zeros((size_y, size_x));
        self.z = loc;
        self.Lx = size_x_mm;
        self.Ly = size_y_mm;

    complexValues = None; # The array of complex-valued fields
    z = 0; # Current z-location of the field
    ident = ''; # String identifier for this field.
    Lx = 0; # Physical size along the x dimension
    Ly = 0; # Physical size along the y dimension

    def fresnelPropagate(self, wavelength, distance):
        """Propagate an electromagnetic field from one location to another using the Fresnel
        u1: Electromagnetic field to propagate at some initial plane
        L: side length of both the initial and final plane, in mm
        z: distance to propagate the field, in mm
        """
        k0 = 2*np.pi/wavelength;
        field_shape = self.complexValues.shape;

        new_field = np.zeros(field_shape);
        M = field_shape[0];
        N = field_shape[1];

        dx = L/N;
        dy = L/M;

        # These are the frequencies we can represent on our discretized grid.
        # They are the 'available' frequency bins we will shove stuff into.
        fx = np.linspace(-1/(2*dx),1/(2*dx),M)
        fy = np.linspace(-1/(2*dy),1/(2*dy),N)

        # This turns our frequencies into a 2D grid for 2-dimensional 
        # fourier transformation.
        (FX,FY) = np.meshgrid(fx, fy);

        # Fresnel transfer function
        H = np.exp(-1j*np.pi*wavelength*z*(np.square(FX)+np.square(FY)));

        # Not sure why this 'fftshift' stuff is necessary.
        H = np.fft.fftshift(H); # shift the transfer function in frequency
        U1 = np.fft.fft2(np.fft.fftshift(self.complexValues))
        U2 = H*U1; # Compute the fourier transform of our new field.
        self.complexValues = np.fft.ifftshift(np.fft.ifft2(U2)); # compute our final electromagnetic field
