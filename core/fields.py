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
# 4. Choose the field sampling so we appropriately oversample our stuff and don't alias anything.

import numpy as np
import matplotlib.pyplot as plt

class Mask:
    """
    The main class from which everything inherits. Represents a two-dimensional
    complex-valued 'mask' to be applied to the electromagnetic field.
    """
    def __init__(self, ident, location, phaseFunction=None, amplitudeFunction=None,
            realFunction=None, imaginaryFunction=None, complex_values = None):

        self.z = location;
        self.ident = ident;
        self.phaseFunction = phaseFunction;
        self.amplitudeFunction = amplitudeFunction;
        self.realFunction = realFunction;
        self.imaginaryFuunction = imaginaryFunction;

        self.complex_values =  complex_values;

    ident = '';
    z = 0;
    Lx = 0;
    Ly = 0;
    phaseFunction = None; # The function of x,y,l which determines the phase
    amplitudeFunction = None; # The function of x, y, l which determines the amplitude
    realFunction = None;
    imaginaryFunction = None;

    complex_values = np.array([]); # Raw complex values for the mask (a 2D array)

    def complexFunction(self, x, y, l):
        """ Return a complex-valued function depending on whether we have
        defined a real/imaginary function, amplitude/phase function, or real/imaginary values.
        """
        if(self.realFunction == None and self.imaginaryFunction == None
                and self.amplitudeFunction != None and self.phaseFunction != None):
            return self.amplitudeFunction(x,y,l)*np.exp(1j*self.phaseFunction(x,y,l))
        elif(self.amplitudeFunction == None and self.phaseFunction == None and
                self.realFunction != None and self.imaginaryFunction != None):
            return self.realFunction(x,y,l) + 1j*self.imaginaryFunction(x,y,l);
        else:
            return None;

    def plotMagPhase(self, title=''):
        """
        Plots the magnitude and phase of the electromagnetic field object
        """
        if(self.complex_values.size > 0):
            fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2,figsize=(7,7))
            ax1.set_xlabel('x(mm)')
            ax1.set_ylabel('y(mm)')
            ax1.set_title('Intensity');

            ax2.set_xlabel('x(mm)')
            ax2.set_ylabel('y(mm)')
            ax2.set_title('Phase');


            if(self.Lx != 0 and self.Ly != 0):
                amp = ax1.imshow(np.square(np.square(np.abs(self.complex_values))),
                        extent=[-self.Lx/2,self.Lx/2,-self.Ly/2,self.Ly/2])
                phase = ax2.imshow(np.angle(self.complex_values),
                        extent=[-self.Lx/2,self.Lx/2,-self.Ly/2,self.Ly/2])
            else:
                amp = ax1.imshow(np.square(np.square(np.abs(self.complex_values)))) 
                phase = ax2.imshow(np.angle(self.complex_values));

            fig.colorbar(amp, ax=ax1)
            fig.colorbar(phase, ax=ax2)
            fig.suptitle(title);

    def saveValues(self, filename):
        """
        Prints the complex-valued field to a file with filename "filename".
        """
        np.savetxt(filename, self.complex_values);



# This is just a container for our complex-valued field at a particular location.
# For now, it has a well-defined wavelength.
class Field(Mask):
    def __init__(self, ident, loc, wavelength, phaseFunction, amplitudeFunction):
        self.ident = ident;
        self.z = loc;
        self.wavelength = wavelength;
        self.phaseFunction = phaseFunction;
        self.amplitudeFunction = amplitudeFunction;

    wavelength = 0; # Wavelength of the monochromatic field.


    def fresnelPropagate(self, distance):
        """Propagate an electromagnetic field from one location to another using the Fresnel
        u1: Electromagnetic field to propagate at some initial plane
        L: side length of both the initial and final plane, in mm
        z: distance to propagate the field, in mm
        """
        if(distance > 0):
            field_shape = self.complex_values.shape;
            print(f"field shape is {field_shape}")
            k0 = 2*np.pi / self.wavelength;

            M = field_shape[0];
            N = field_shape[1];
            dx = self.Lx/N;
            dy = self.Ly/M;

            print(f"dx: {dx}, dy: {dy}");

            # This is our oversampling ratio for the transfer function propagator.
            # ideally, we want this to be much greater than one for appropriate sampling
            # of our phase function. We want this number to be much greater than 1.

            # These are the frequencies we can represent on our discretized grid.
            # They are the 'available' frequency bins we will shove stuff into.
            # FX, FY turns our frequencies into a 2D grid for 2-dimensional 
            # fourier transformation.
            fx = np.linspace(-1/(2*dx),1/(2*dx),M)
            fy = np.linspace(-1/(2*dy),1/(2*dy),N)
            (FX,FY) = np.meshgrid(fx, fy);

            # Fresnel transfer function
            H = np.exp(-1j*np.pi*self.wavelength*distance*(np.square(FX) + np.square(FY)));

            # Not sure why this 'fftshift' stuff is necessary.
            H = np.fft.fftshift(H); # shift the transfer function in frequency
            U1 = np.fft.fft2(np.fft.fftshift(self.complex_values))
            U2 = H*U1; # Compute the fourier transform of our new field.

            # update our complex values
            self.complex_values = np.fft.ifftshift(np.fft.ifft2(U2)); # compute our final electromagnetic field
            self.z += distance;

            print(self.wavelength);
            oversampling_ratio = dx * self.Lx / self.wavelength / distance;
            return oversampling_ratio;
        else:
            return 0;

def circ(diameter, x,y):
    """
    Function that returns 1 if the coordinates x, y are within the radius of the circle
    """
    return np.heaviside(diameter/2 - np.sqrt(np.square(x) + np.square(y)),1);

def rect(width, height, x, y):
    """
    Function that returns 1 if the coordinates x, y are within the rectangle.
    """
    return np.heaviside(width/2 - np.abs(x), 1)*np.heaviside(height/2 - np.abs(y), 1);

def sq(width, x, y):
    """
    Function that returns 1 if the coordinates x, y are within the square
    """
    return rect(width, width, x, y);
