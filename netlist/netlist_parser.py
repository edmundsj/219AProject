# This is a parser for our simple netlists.
# TODO:
# 1. Incorporate phase information into Gaussian wave based on the beam waist.
# 2. Switch wavelength to a global property, and only allow people to define a single
#   wavelength per simulation. They can still do multiple wavelengths eventually, but 
#   for now this should be a global setting.


# EXPECTED NETLIST FORMAT
GAUSSIAN_WAIST_POSITION = 3;
WAVELENGTH_POSITION = 2;
PLANE_WAVE_ANGLE_POSITION = 3;
APERTURE_WIDTH = 2;
APERTURE_HEIGHT = 3;

import re
import numpy as np
from core.fields import * # import our core objects that will be sent back to our program.

class NetlistParser:
    filename = '';

    def __init__(self, filename):
        self.filename = filename;

    # Stripts the units after the number
    def stripUnits(self, text):
        if(text[-2:] == 'mm'):
            return float(text[:-2]);
        elif(text[-2:] == 'um'):
            return 0.001*float(text[:-2]);
        elif(text[-2:] == 'nm'):
            return 1e-6*float(text[:-2]);
        else:
            return float(text);

    def parseNetlist(self):
        fields = []; # 
        masks = [];
        f = open(self.filename, 'r');
        sep = '#'; # The separator for our netlist files. Also comment symbol.
        f1 = f.readlines();
        processed_lines = [];
        for line in f1:
            comments_removed = line.split(sep, 1)[0];   # Removes all comments
            stripped = comments_removed.strip();        # Removes trailing and preceding whitespace
            commas_removed = stripped.replace(',', ' ');  # Remove all commas
            whitespace_removed = re.sub('\s+', ' ', commas_removed);
            if(stripped != ''):
                processed_lines.append(whitespace_removed);

        # Loop through all the lines that actually contain something interesting
        for line in processed_lines:
            line_chunks = line.split(' ');
            name = line_chunks[0];
            location_text = line_chunks[1];
            location = self.stripUnits(location_text);


            # If we have a plane wave source, we need to know two pieces of information:
            # the intensity of the plane wave, and its tilt. By default (and for now)
            # we will assume the tilt is zero. If it is not, this just changes the phase
            # function we need to apply to the plane wave.
            # This is expecting a grid of x and y values. It will return a numpy array
            # Even if you feed it with a single value.
            if(line[0] == 'W'): # The line contains a plane-wave source

                wavelength = self.stripUnits(line_chunks[WAVELENGTH_POSITION]);

                phase_function = lambda x, y, l: np.zeros(x.shape);
                amplitude_function = lambda x, y, l: np.ones(x.shape);

                new_planewave = Field(name, location, wavelength, phase_function, amplitude_function);
                fields.append(new_planewave);

            # This needs to be modified to include the Gaussian phase.
            elif(line[0] == 'G'):

                wavelength = self.stripUnits(line_chunks[WAVELENGTH_POSITION]);
                waist = self.stripUnits(line_chunks[GAUSSIAN_WAIST_POSITION]);

                phase_function = lambda x, y, l: np.zeros(x.shape);
                intensity_function = lambda x, y, l: 1* np.exp(-(np.square(x)+np.square(y))/np.square(waist));

                new_gaussian_beam = Field(name, location, wavelength, phase_function, intensity_function);
                fields.append(new_gaussian_beam);

            elif(line[0] == 'L'): # Thin lens - note that the phase mask actually depends on the wavelength.

                focal_length = self.stripUnits(line_chunks[2]);

                phase_function = lambda x, y, l: 2*np.pi / (2*focal_length * l) * (np.square(x) + np.square(y));
                amplitude_function = lambda x, y, l: np.ones(x.shape);

                new_lens = Mask(name, location, phase_function, amplitude_function);
                masks.append(new_lens);

            elif(line[0] == 'P'): # The line contains a plane definition
                new_plane = Mask(name, location);
                masks.append(new_plane);

            elif(line[0] == 'C'): # Circular aperture 
                diameter = self.stripUnits(line_chunks[APERTURE_WIDTH]);

                amplitude_function = lambda x, y, l: circ(diameter, x, y);
                phase_function = lambda x, y, l: np.zeros(x.shape);
                new_aperture = Mask(name, location, phase_function, amplitude_function);
                masks.append(new_aperture);

            elif(line[0] == 'R'): # Rectangular aperture

                width = self.stripUnits(line_chunks[APERTURE_WIDTH]);
                height = self.stripUnits(line_chunks[APERTURE_HEIGHT]);

                amplitude_function = lambda x, y, l: rect(width, height, x, y);
                phase_function = lambda x, y, l: np.zeros(x.shape);
                new_aperture = Mask(name, location, phase_function, amplitude_function);
                masks.append(new_aperture);

            elif(line[0] == 'S'): # Square aperture

                width = self.stripUnits(line_chunks[APERTURE_WIDTH]);

                amplitude_function = lambda x, y, l: sq(width, x, y);
                phase_function = lambda x, y, l: np.zeros(x.shape);
                new_aperture = Mask(name, location, phase_function, amplitude_function);

                masks.append(new_aperture);


            else:
                print("ERROR: Not able to parse line:")
                print(line)
                print("Object type not supported");

        masks = sorted(masks, key=lambda x: x.z);
        fields = sorted(fields, key=lambda x: x.z);

        return [masks, fields];

