# This is a parser for our simple netlists.
# TODO:
# 1. Incorporate phase information into Gaussian wave based on the beam waist.
# 2. Switch wavelength to a global property, and only allow people to define a single
#   wavelength per simulation. They can still do multiple wavelengths eventually, but 
#   for now this should be a global setting.


import re
import numpy as np
from core.fields import Mask # import our core objects that will be sent back to our program.

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
        wavelength = None;
        k0 = None;
        fields_list = []; # 
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

        # First, loop through all the lines to extract the wavelength
        for line in processed_lines:
            line_chunks = line.split(' ');
            name = line_chunks[0];
            if(line[0] == '.'):
                wavelength = self.stripUnits(line_chunks[1]);
                k0 = 2*np.pi / wavelength;


        # Loop through all the lines that actually contain something interesting
        for line in processed_lines:
            line_chunks = line.split(' ');
            name = line_chunks[0];
            location_text = line_chunks[1];
            location = self.stripUnits(location_text);

            print(line);
            print(line_chunks);

            if(line[0] == 'P'): # The line contains a plane definition
                print("Plane! No applied phase mask.");

                phase_function = None; # Do not apply any phase mask
                new_plane = Mask(name, location, False, phase_function);
                fields_list.append(new_plane);

            # If we have a plane wave source, we need to know two pieces of information:
            # the intensity of the plane wave, and its tilt. By default (and for now)
            # we will assume the tilt is zero. If it is not, this just changes the phase
            # function we need to apply to the plane wave.
            elif(line[0] == 'W'): # The line contains a plane-wave source
                print("Plane Wave Source!");

                phase_function = lambda x, y: 0;
                intensity_function = lambda x, y: 1;

                new_planewave = Mask(name, location, True, phase_function, intensity_function);
                fields_list.append(new_planewave);

            elif(line[0] == 'G'):
                print("Gaussian Source!");

                waist = self.stripUnits(line_chunks[2]);

                phase_function = lambda x, y: 0; # FOR NOW, the Gaussian 
                intensity_function = lambda x, y: 1* np.exp(-(np.square(x)+np.square(y))/np.square(waist));

                new_gaussian_beam = Mask(name, location, True, phase_function, intensity_function);
                fields_list.append(new_gaussian_beam);

            elif(line[0] == 'L'): # Thin lens - note that the phase mask actually depends on the wavelength.

                focal_length = self.stripUnits(line_chunks[2]);

                phase_function = lambda x, y: k0 / (2*focal_length) * (np.square(x) + np.square(y));
                intensity_function = lambda x, y: 1;

                new_lens = Mask(name, location, False, phase_function, intensity_function);
                fields_list.append(new_lens);
                print("Lens!");

            elif(line[0] == '.'): # Simulation directive.
                print("Simulation directive!");

            else:
                print("Not able to parse line:")
                print(line)
                print("Object type not supported");

        return fields_list;
