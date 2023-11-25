#!/usr/local/bin/python
#
# This code uses several numerical integration techniques. These are all
# detailed in:
#
# Marron, M. J. (1982). Composite Rules and Romberg Integration.
# Numerical Analysis: A Practical Approach. 
# Macmillan Publishing Company, Incorporated.
import ConfigParser
import datetime
import getopt
import math
import os
import random
import sys

import numarray

class WindowSpacingSequence(object):
    """This knows the next starting offset.
    
    Each window asks this object for the triggering timestep. This way the system always has a uniform window spacing.
    """
    def __init__(self, window_spacing):
        self.spacing = window_spacing
        self.timestep = 0
        
    def get_step(self):
        current_timestep = self.timestep
        self.timestep = self.timestep + self.spacing
        return current_timestep
        
    start = property(get_step)
        
class SlidingWindow(object):
    def __init__(self, initial_length):
        """Sliding window constructor.

        Builds a single sliding window. _offset is the index (0,len(self)) of the internal array. _last_timestep is the last
        timestep that needs to be processed. Once this timestep is processed the data needs to be accumulated and then a new
        start needs to be requested.
        """
        self._offset = 0;
        self._last_timestep = 0
        self._start = 0;
        self.p0xy = 0.0;
        self.p0xz = 0.0;
        self.p0yz = 0.0;
        self.pxy = None
        self.pxz = None
        self.pyz = None
        self.length = initial_length

    def __str__(self):
        """String representation of a sliding window.

        """
        return "current slide %s\noffset %s\nstart %s\nend %s\np0xy %s\np0xz %s\np0yz %s\n" % (self.cur_slide, self.data_offset, self.start, self.end, self.p0xy, self.p0xz, self.p0yz)
    
    def clear(self):
        """Set the arrays to 0.0. Reset the offset.
        """
        self.pxy *= 0.0
        self.pxz *= 0.0
        self.pyz *= 0.0
        self._offset = 0
            
    def get_length(self):
        return len(self.pxy)
        
    def set_length(self, value):
        """Sets the internal arrays by reallocating a zero-initialized numarray.
        
        This wipes out anything you had in there before.
        """
        self.pxy = numarray.zeros(value, numarray.Float64)
        self.pxz = numarray.zeros(value, numarray.Float64)
        self.pyz = numarray.zeros(value, numarray.Float64)
        
    length = property(get_length, set_length)

    def set_start(self, value):
        """Set the starting timestep (t0).
        
        When the starting time step is set the final timestep is also computed using the length of the window.
        """
        self._start = value
        self._last_timestep = self.start + self.length - 1
        
    def get_start(self):
        return self._start
        
    start = property(get_start, set_start)
    
    def process_pressures(self, timestep, pxy, pxz, pyz, accumxy, accumxz, accumyz, number_of_origins):
        """Computes the pressure product, slides the window, and accumulates.
        
        
        """
        if timestep >= self._start:
            # If we're at the beginning then we need to set our t0 values.
            if timestep == self._start:
                self.p0xy = pxy
                self.p0xz = pxz
                self.p0yz = pyz
    
            # Store the value in the arrays.
            self.pxy[self._offset] = pxy
            self.pxz[self._offset] = pxz
            self.pyz[self._offset] = pyz
            self._offset = self._offset + 1
            
            # If we're at the last timestep to be processed then compute the p0xy * pxy product for the whole array and then accumulate
            # it into the overall statistics. Once that's done, clear everything and get a new starting timestep.
            if timestep == self._last_timestep:
                self.pxy *= self.p0xy
                self.pxz *= self.p0xz
                self.pyz *= self.p0yz
                
                accumxy += self.pxy
                accumxz += self.pxz
                accumyz += self.pyz
                
                number_of_origins[0] = number_of_origins[0] + 1
                
                return True
        
        # This is what happens most of the time, we haven't accumulated the entire window.
        return False

def estimate_viscosity(timestep, acf_xy_data, acf_xz_data, acf_yz_data, visc_xy_data, visc_xz_data, visc_yz_data):
    """Finds the viscosity using the asymptotic value of the ACF.
    
    The viscosity for a given Cartesian direction is the value of the integrated ACF function at the timestep which corresponds to
    where the slowest ACF (of the three Cartesian directions) falls below 0.1% of its value at t = 0.
    """
    contribution_limit = 1.0 / 1000.0

    # Time step at which the match is found
    ts_acf_xy = ts_acf_xz = ts_acf_yz = 0.0

    # Find the viscosity
    firstline = True
    for index in range(len(acf_xy_data)):
        if firstline:
            firstline = False
            threshold_acf_xy = acf_xy_data[index] * contribution_limit
            threshold_acf_xz = acf_xz_data[index] * contribution_limit
            threshold_acf_yz = acf_yz_data[index] * contribution_limit

        if acf_xy_data[index] < threshold_acf_xy and ts_acf_xy == 0.0:
            ts_acf_xy = float(index) * timestep
            visc_xy = visc_xy_data[index]
            visc_xz = visc_xz_data[index]
            visc_yz = visc_yz_data[index]
            if ts_acf_xy * ts_acf_xz * ts_acf_yz != 0.0:
                break

        if acf_xz_data[index] < threshold_acf_xz and ts_acf_xz == 0.0:
            ts_acf_xz = float(index) * timestep
            visc_xy = visc_xy_data[index]
            visc_xz = visc_xz_data[index]
            visc_yz = visc_yz_data[index]
            if ts_acf_xy * ts_acf_xz * ts_acf_yz != 0.0:
                break

        if acf_yz_data[index] < threshold_acf_yz and ts_acf_yz == 0.0:
            ts_acf_yz = float(index) * timestep
            visc_xy = visc_xy_data[index]
            visc_xz = visc_xz_data[index]
            visc_yz = visc_yz_data[index]
            if ts_acf_xy * ts_acf_xz * ts_acf_yz != 0.0:
                break
    else:
        visc_xy = 0.0
        visc_xz = 0.0
        visc_yz = 0.0
        
    return visc_xy, visc_xz, visc_yz

if __name__ == "__main__":
    # Check to see if we've passed a config file. Otherwise ask for one.
    optlist, args = getopt.getopt(sys.argv, '') 
    if len(args) != 3:
        print "usage: viscosity.py configfilename inputfilename"
        sys.exit(1)

    # Open the configuration file for reading.
    config = ConfigParser.ConfigParser()
    config.readfp(open(args[1]))

    # If a single boxsize is not present (the typical case) then we MUST have
    # boxsizes in all of the directions.
    try:
        boxsize = config.getfloat("InitialConditions", "boxsize")
        boxsizeX = boxsize
        boxsizeY = boxsize
        boxsizeZ = boxsize
        
    except ConfigParser.NoOptionError:
        try:
            boxsizeX = config.getfloat("InitialConditions", "boxsize_x")
            boxsizeY = config.getfloat("InitialConditions", "boxsize_y")
            boxsizeZ = config.getfloat("InitialConditions", "boxsize_z")
        except ConfigParser.NoOptionError:
            print "No boxsize present, must be a single boxsize or one of each."
            sys.exit(1)
            
    # Get the viscosity computational parameters.
    skip = config.getint("Viscosity", "skip")
    timestep = config.getfloat("Viscosity", "timestep")
    outputFilename = config.get("Viscosity", "filename")
    process_num_steps = config.getint("Viscosity", "process")
    pxy_col = config.getint("Viscosity", "pxy_col")
    pxz_col = config.getint("Viscosity", "pxz_col")
    pyz_col = config.getint("Viscosity", "pyz_col")
    temperature_col = config.getint("Viscosity", "temperature_col")
    window_length = config.getint("Viscosity", "window_length")
    window_spacing = config.getint("Viscosity", "window_spacing")
    convert_to_pa = config.getfloat("Viscosity", "convert_to_pa")
    
    # Open the I/O files.
    infile = file(args[2], "r")
    outfile = file(outputFilename, "w")
    
    # Create a window spacing. This will keep track of where the next offset is located.
    window_spacing_sequence = WindowSpacingSequence(window_spacing)
    
    # Figure out the minimum number of SlidingWindows that we'll need to properly compute the ACF. We use multiple sliding
    # windows because it would waste time to keep rewinding the data file to slide one window. That said, since we're setting
    # our start times independently we are not going to get duplication of coverage even if we have "too many" windows. Note that 
    # if the window spacing is greater than the window width you only need one sliding window since there is no overlap 
    # between windows.
    if window_spacing > window_length:
        number_of_windows = 1
    else:
        number_of_windows = int(window_length / window_spacing) + 1
        
    # Create the sliding window objects. Once created, ask for 
    sw = []
    for i in range(number_of_windows):
        one_window = SlidingWindow(window_length)
        one_window.start = window_spacing_sequence.start
        one_window.clear()
        sw.append(one_window)

    # Create the arrays that will hold the answers. For the scalar we MUST us a one-item list so the function can change the
    # value, a straight number will not work...
    accumxy = numarray.zeros(window_length, numarray.Float64)
    accumxz = numarray.zeros(window_length, numarray.Float64)
    accumyz = numarray.zeros(window_length, numarray.Float64)
    number_of_origins = [0]
     
    # Skip input lines until we get to where we want to go.
    for i in range(skip):
        dummy = infile.readline()
    
    temperature_sum = 0.0
    # Process the right number of input lines. The section also converts the
    # pressure to Pa from the input units.
    for i in range(process_num_steps):
        splitline = infile.readline().split()
        temperature_sum += float(splitline[temperature_col - 1])
        pxy = float(splitline[pxy_col - 1]) * convert_to_pa
        pxz = float(splitline[pxz_col - 1]) * convert_to_pa
        pyz = float(splitline[pyz_col - 1]) * convert_to_pa
        for wsw in sw:
            if wsw.process_pressures(i, pxy, pxz, pyz, accumxy, accumxz, accumyz, number_of_origins):
                wsw.clear()
                wsw.start = window_spacing_sequence.start
            
    # Done accumulating data, print out average temperature and the scaled ACF. Note that the number_of_origins differs from the
    # process_num_steps because of the window_spacing.
    average_temperature = temperature_sum / float(process_num_steps)

    accumxy /= float(number_of_origins[0])
    accumxz /= float(number_of_origins[0])
    accumyz /= float(number_of_origins[0])

    # Compute constant part of integral.
    volume = boxsizeX * boxsizeY * boxsizeZ * 1e-30
    kb = 1.3806505e-23
    uconv = volume / (kb * average_temperature)

    # -------------- RECTANGULAR INTEGRATION --------------
    # Rectangular integration simply multiplies the value of f(x) by dt to get
    # a rectangular area under the curve.
    rect_ppxy = accumxy * timestep
    rect_ppxz = accumxz * timestep
    rect_ppyz = accumyz * timestep
    rect_ave = ((accumxy + accumxz + accumyz) / 3.0) * timestep
    
    for i in range(1, window_length):
        rect_ppxy[i] += rect_ppxy[i - 1] 
        rect_ppxz[i] += rect_ppxz[i - 1] 
        rect_ppyz[i] += rect_ppyz[i - 1] 
        rect_ave[i] += rect_ave[i - 1]

    # -------------- COMPOSITE TRAPEZOIDAL INTEGRATION --------------
    # This is equation 2a in section 7.4. This is applied in an iterative manner
    # by continually moving the termination point of the integral out so a plot
    # can be made of the integration.
    trap_ppxy = numarray.zeros(window_length, numarray.Float64)
    trap_ppxz = numarray.zeros(window_length, numarray.Float64)
    trap_ppyz = numarray.zeros(window_length, numarray.Float64)
    trap_ave = numarray.zeros(window_length, numarray.Float64)
    
    for endpoint in range(window_length):
        endpoint_ave_xy = (accumxy[0] + accumxy[endpoint]) / 2.0
        endpoint_ave_xz = (accumxz[0] + accumxz[endpoint]) / 2.0
        endpoint_ave_yz = (accumyz[0] + accumyz[endpoint]) / 2.0
        interior_sum_xy = 0.0
        interior_sum_xz = 0.0
        interior_sum_yz = 0.0
        if endpoint > 1:
            for interior_point in range(1, endpoint - 1):
                interior_sum_xy += accumxy[interior_point]
                interior_sum_xz += accumxz[interior_point]
                interior_sum_yz += accumyz[interior_point]
        trap_ppxy[endpoint] = (endpoint_ave_xy + interior_sum_xy) * timestep
        trap_ppxz[endpoint] = (endpoint_ave_xz + interior_sum_xz) * timestep
        trap_ppyz[endpoint] = (endpoint_ave_yz + interior_sum_yz) * timestep
        trap_ave[endpoint] = (trap_ppxy[endpoint] + trap_ppxz[endpoint] + trap_ppyz[endpoint]) / 3.0
        
    # -------------- COMPOSITE SIMPSON'S INTEGRATION --------------
    # This is equation 3a in section 7.4. This is applied in an iterative manner
    # by continually moving the termination point of the integral out so a plot
    # can be made of the integration.
    simp_ppxy = numarray.zeros(window_length, numarray.Float64)
    simp_ppxz = numarray.zeros(window_length, numarray.Float64)
    simp_ppyz = numarray.zeros(window_length, numarray.Float64)
    simp_ave = numarray.zeros(window_length, numarray.Float64)
    
    for endpoint in range(window_length):
        even_interior_sum_xy = 0.0
        even_interior_sum_xz = 0.0
        even_interior_sum_yz = 0.0
        odd_interior_sum_xy = 0.0
        odd_interior_sum_xz = 0.0
        odd_interior_sum_yz = 0.0
        if endpoint > 1:
            for interior_point in range(1, endpoint - 1):
                # Accumulate even and odd data
                if interior_point % 2 == 0:
                    even_interior_sum_xy += accumxy[interior_point]
                    even_interior_sum_xz += accumxz[interior_point]
                    even_interior_sum_yz += accumyz[interior_point]
                else:
                    odd_interior_sum_xy += accumxy[interior_point]
                    odd_interior_sum_xz += accumxz[interior_point]
                    odd_interior_sum_yz += accumyz[interior_point]
                    
        simp_ppxy[endpoint] = ((accumxy[0] + accumxy[endpoint] + 2 * even_interior_sum_xy + 4 * odd_interior_sum_xy) * timestep) / 3.0
        simp_ppxz[endpoint] = ((accumxz[0] + accumxz[endpoint] + 2 * even_interior_sum_xz + 4 * odd_interior_sum_xz) * timestep) / 3.0
        simp_ppyz[endpoint] = ((accumyz[0] + accumyz[endpoint] + 2 * even_interior_sum_yz + 4 * odd_interior_sum_yz) * timestep) / 3.0
        simp_ave[endpoint] = (simp_ppxy[endpoint] + simp_ppxz[endpoint] + simp_ppyz[endpoint]) / 3.0
        
    # -------------- ESTIMATE VISCOSITY --------------
    rect_est_visc_xy, rect_est_visc_xz, rect_est_visc_yz = estimate_viscosity(timestep, accumxy, accumxz, accumyz, rect_ppxy, rect_ppxz, rect_ppyz)
    trap_est_visc_xy, trap_est_visc_xz, trap_est_visc_yz = estimate_viscosity(timestep, accumxy, accumxz, accumyz, trap_ppxy, trap_ppxz, trap_ppyz)
    simp_est_visc_xy, simp_est_visc_xz, simp_est_visc_yz = estimate_viscosity(timestep, accumxy, accumxz, accumyz, simp_ppxy, simp_ppxz, simp_ppyz)

    # -------------- OUTPUT RESULTS --------------
    
    # Viscosity calculated with:
    #   3.027996e-27 m^3 Volume.
    #   4.158906e+03 K Average Temperature.
    #   Filler to keep header backward compatible (windows_in_set).
    #   99500 Time Origins.
    #   5000 Timestep wide Data Window.
    #   10 Timesteps between time origins.
    # Estimated viscosity is 3.540758e-03 Pa*s .
    #

    # Output header.
    outfile.write("# Viscosity calculated with:\n")
    outfile.write("#   %e m^3 Volume.\n" % volume)
    outfile.write("#   %e K Average Temperature.\n" % average_temperature)
    outfile.write("#   %d No, Time Origins.\n" % number_of_origins[0])
    outfile.write("#   %d tW, Timestep wide Data Window.\n" % window_length)
    outfile.write("#   %d tS, Timesteps between time origins.\n" % window_spacing)
    outfile.write("# Estimated viscosity:\n")
    
    outfile.write("#   %e Pa*s (Rectangular Integration xy).\n" % (rect_est_visc_xy * uconv,))
    outfile.write("#   %e Pa*s (Rectangular Integration xz).\n" % (rect_est_visc_xz * uconv,))
    outfile.write("#   %e Pa*s (Rectangular Integration yz).\n" % (rect_est_visc_yz * uconv,))
    rect_ave_visc = (rect_est_visc_xy + rect_est_visc_xz + rect_est_visc_yz) / 3.0
    rect_ave_visc_err = math.sqrt(((rect_est_visc_xy - rect_ave_visc)**2 + (rect_est_visc_xz - rect_ave_visc)**2 + (rect_est_visc_yz - rect_ave_visc)**2) / 3.0)
    outfile.write("#   %e +/- %e Pa*s (Rectangular Integration Average).\n" % (rect_ave_visc * uconv, rect_ave_visc_err * uconv))

    outfile.write("#   %e Pa*s (Trapezoidal Integration xy).\n" % (trap_est_visc_xy * uconv,))
    outfile.write("#   %e Pa*s (Trapezoidal Integration xz).\n" % (trap_est_visc_xz * uconv,))
    outfile.write("#   %e Pa*s (Trapezoidal Integration yz).\n" % (trap_est_visc_yz * uconv,))
    trap_ave_visc = (trap_est_visc_xy + trap_est_visc_xz + trap_est_visc_yz) / 3.0
    trap_ave_visc_err = math.sqrt(((trap_est_visc_xy - trap_ave_visc)**2 + (trap_est_visc_xz - trap_ave_visc)**2 + (trap_est_visc_yz - trap_ave_visc)**2) / 3.0)
    outfile.write("#   %e +/- %e Pa*s (Trapezoidal Integration Average).\n" % (trap_ave_visc * uconv, trap_ave_visc_err * uconv))

    outfile.write("#   %e Pa*s (Simpson's Integration xy).\n" % (simp_est_visc_xy * uconv,))
    outfile.write("#   %e Pa*s (Simpson's Integration xz).\n" % (simp_est_visc_xz * uconv,))
    outfile.write("#   %e Pa*s (Simpson's Integration yz).\n" % (simp_est_visc_yz * uconv,))
    simp_ave_visc = (simp_est_visc_xy + simp_est_visc_xz + simp_est_visc_yz) / 3.0
    simp_ave_visc_err = math.sqrt(((simp_est_visc_xy - simp_ave_visc)**2 + (simp_est_visc_xz - simp_ave_visc)**2 + (simp_est_visc_yz - simp_ave_visc)**2) / 3.0)
    outfile.write("#   %e +/- %e Pa*s (Simpson's Integration Average).\n" % (simp_ave_visc * uconv, simp_ave_visc_err * uconv))
    outfile.write("#\n")
    
    outfile.write("timestep acf_xy acf_xz acf_yz rect_ppxy rect_ppxz rect_ppyz rect_ave trap_ppxy trap_ppxz trap_ppyz trap_ave simp_ppxy simp_ppxz simp_ppyz simp_ave\n")
    for i in range(window_length):
        outfile.write("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n" % (timestep * float(i),
        accumxy[i], accumxz[i], accumyz[i],
        rect_ppxy[i] * uconv, rect_ppxz[i] * uconv, rect_ppyz[i] * uconv, rect_ave[i] * uconv,
        trap_ppxy[i] * uconv, trap_ppxz[i] * uconv, trap_ppyz[i] * uconv, trap_ave[i] * uconv,
        simp_ppxy[i] * uconv, simp_ppxz[i] * uconv, simp_ppyz[i] * uconv, simp_ave[i] * uconv
        ))
        
 
