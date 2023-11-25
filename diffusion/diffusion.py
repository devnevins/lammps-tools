#!/usr/bin/python
import ConfigParser
import getopt
import sys
import csv 

class Material:
    """Internal material class which eases getting information about the materials used.
    
    """
    def __init__(self, index=-1, number=0, symbol='', mass=0.0, charge=0.0):
        """material constructor.
        """
        self.index = index
        self.number = number
        self.symbol = symbol
        self.mass = mass
        self.charge = charge

    def __str__(self):
        """String representation of a piece of MD Material.
        """
        return "%s particles of %s. Mass: " % ( self.number, self.symbol, self.mass )

def convertMass(mass):
    """Converts mass from g/mol to kg/particle needed by the velocity calcs.
    """        
    return ((mass * 1e-3 )/ 6.022e23)

def linreg(X, Y):
    """
    Computes linear regression given lists of X and Y.
    
    Taken and modified slightly from William Park who got it from numerical recipies in C. See: http://www.phys.uu.nl/~haque/computing/WPark_recipes_in_python.html

    Returns coefficients to the regression line "y=ax+b" from x[] and
    y[].  Basically, it solves 
        Sxx a + Sx b = Sxy
         Sx a +  N b = Sy
    where Sxy = \sum_i x_i y_i, Sx = \sum_i x_i, and Sy = \sum_i y_i.  The
    solution is
        a = (Sxy N - Sy Sx)/det
        b = (Sxx Sy - Sx Sxy)/det
    where det = Sxx N - Sx^2.  In addition,
        Var|a| = s^2 |Sxx Sx|^-1 = s^2 | N  -Sx| / det
           |b|       |Sx  N |          |-Sx Sxx|
        s^2 = {\sum_i (y_i - \hat{y_i})^2 \over N-2}
            = {\sum_i (y_i - ax_i - b)^2 \over N-2}
            = residual / (N-2)
        R^2 = 1 - {\sum_i (y_i - \hat{y_i})^2 \over \sum_i (y_i - \mean{y})^2}
            = 1 - residual/meanerror
    
    It also prints to <stdout> few other data,
        N, a, b, R^2, s^2,
    which are useful in assessing the confidence of estimation.
    """
    from math import sqrt
    if len(X) != len(Y):  raise ValueError, 'unequal length'

    N = len(X)
    Sx = Sy = Sxx = Syy = Sxy = 0.0
    for x, y in map(None, X, Y):
        Sx = Sx + x
        Sy = Sy + y
        Sxx = Sxx + x*x
        Syy = Syy + y*y
        Sxy = Sxy + x*y

    det = Sxx * N - Sx * Sx
    a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det

    meanerror = residual = 0.0
    for x, y in map(None, X, Y):
        meanerror = meanerror + (y - Sy/N)**2
        residual = residual + (y - a * x - b)**2
    RR = 1 - residual/meanerror
    ss = residual / (N-2)
    Var_a, Var_b = ss * N / det, ss * Sxx / det
    
    #print "y=ax+b"
    #print "N= %d" % N
    #print "a= %g \\pm t_{%d;\\alpha/2} %g" % (a, N-2, sqrt(Var_a))
    #print "b= %g \\pm t_{%d;\\alpha/2} %g" % (b, N-2, sqrt(Var_b))
    #print "R^2= %g" % RR
    #print "s^2= %g" % ss
    
    return a, b

if __name__ == "__main__":
    # Check to see if we've passed a config file. Otherwise ask for one.
    optlist, args = getopt.getopt(sys.argv, '') 
    if len(args) != 5:
        print "usage: diffusion.py msd_filename.txt timestep start_time outfilename.csv"
        sys.exit(1)

    msd = []
    ts_arr = []
    ts = 0
    # Load the entire data set. We only want the total displacement, not any of the cartesian directions.
    for line in file(args[1]):
        msd_arr = line.split()
        try:
            msd.append(float(msd_arr[4]))
            ts_arr.append(int(ts))
            ts = ts + 1
        except:
            continue
        
    # Compute where the diffusion coefficient calculation starts. We do it in terms of time so we need to convert that into timesteps.
    timestep = float(args[2])
    start_time = float(args[3])

    # Compute the linear regression starting from the start_time.
    
    start = int(start_time / timestep)
    #for a in range(len(msd)):
    #    print msd[a]
    
    # Truncate lists to reflect where we want to compute linear regression.
    ts_arr = ts_arr[start:]
    msd = msd[start:]
    (m, b) = linreg(ts_arr, msd)
    
    # Convert slope into D and fix units
    # Initial units: Angstroms^2 fsec^-1
    # Desired units: m^2 sec-1
    d = (m * 1.0e-5) / 6.0    
        
    writer = csv.writer(open(args[4], "wb")) 
    header = [ 'Diffusion (m^2 sec^-1)', 'Slope (Angstroms^2 fsec^-1)', 'Intercept (Angstroms)' ]
    data = [ d , m, b ]
    writer.writerow(header) 
    writer.writerow(data) 
