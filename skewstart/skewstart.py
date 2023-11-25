#!/usr/bin/python
import ConfigParser, os, getopt, sys, math, random, datetime

# Material class used internally by skewstart. This didn't seem to have much utility outside
# of this code so I included it.
class Material:
    def __init__(self, index=-1, number=0, name='', symbol='', mass=0.0, charge=0.0):
        """material constructor.

        """
        self.index = index
        self.number = number
        self.symbol = symbol
        self.name = name
        self.mass = mass
        self.charge = charge

    def __str__(self):
        """String representation of a piece of MD Material.

        """
        return "%s particles of %s. Mass: %s" % ( self.number, self.symbol, self.mass )

class OutType:
    def __init__(self):
        assert self.__class__ is not OutType, \
                'OutType is an abstract class, do not instantiate it.'

    def output(self, boxsizeX, boxsizeY, boxsizeZ, pos, vel, filename="outfile", location="none", material="none"):
        assert self.__class__ is not OutType, \
                'OutType is an abstract class, you must override its methods.'

    def convertMass(self, mass):
        return mass

    def convertCharge(self, charge):
        return charge

    def write(self, msg):
        """Writes the message out to the correct place, a file or the screen.

        """        
        if self.isScreen:
            print msg 
        if self.isFile:
            self.outfile.write(msg + "\n")
            self.outfile.flush()

    def writePotentialParameters(self, numMaterials):
        """Print out the potential parameters in the skewstart generated
        output file. Not always required.
        """
        pass
        
    def setOutputDestination(self, filename, location):
        """Sets the outputs destination, either a file or screen or both.

        """        
        if location == "none":
            print "filename: %s, location: %s" % ( filename, location )

        self.isFile = 0
        self.isScreen = 0
        
        if location == "both" or location == "file":
            self.outfile = open(filename, "w")
            self.isFile = 1

        if location == "both" or location == "screen":
            self.isScreen = 1

class OutTypeLAMMPSG(OutType):
    """Takes MKS pos and vel data and outputs a LAMMPSG compatible file.

    """    
    def __init__(self):
        pass

    def output(self, boxsizeX, boxsizeY, boxsizeZ, pos, vel, filename="outfile", location="none", material="none"):
        self.setOutputDestination(filename, location)

        halfBoxsizeX = boxsizeX / 2.0
        halfBoxsizeY = boxsizeY / 2.0
        halfBoxsizeZ = boxsizeZ / 2.0
        
        # Compute the total number of particles (I know this is repetitive)
        numParticles = 0
        for aMaterial in material:
            numParticles += aMaterial.number

        # Output the config file header.
        # In LAMMPSG (C++) the material starts with 1, not 0.
        self.write("Skew-start Configuration File Generated %s\n" % datetime.datetime(datetime.MINYEAR, 1, 1).now())
        self.write("%s atoms\n0 bonds\n0 angles\n0 dihedrals\n0 impropers\n" % numParticles)
        self.write("%(0)s atom types\n" % {'0' : len(material)})
        self.write("%(1)11.8E %(2)11.8E xlo xhi" %  {'1' : -halfBoxsizeX, '2' : halfBoxsizeX})
        self.write("%(1)11.8E %(2)11.8E ylo yhi" %  {'1' : -halfBoxsizeY, '2' : halfBoxsizeY})
        self.write("%(1)11.8E %(2)11.8E zlo zhi\n" %  {'1' : -halfBoxsizeZ, '2' : halfBoxsizeZ})
        self.write("Masses\n")
        for i in range(len(material)):
            self.write("%d %.3f" % (i + 1, (material[i].mass * 1e3 * 6.022e23)))
        self.writePotentialParameters(len(material))
        self.write("\nAtoms\n")

        # Output the position in Angstroms
        groupNumber = 0 # Only one group for all particles.
        particleNumber = 1
        for p in pos:
            self.write("%d\t%d\t%d\t%.4f\t%.17e %.17e %.17e" % (particleNumber, groupNumber, p[0], material[p[0] - 1].charge, p[1] - halfBoxsizeX, p[2] - halfBoxsizeY, p[3] - halfBoxsizeZ))
            particleNumber += 1

        # Output the velocities
        self.write("\nVelocities\n")
        particleNumber = 1
        for v in vel:
            self.write("%d\t%.17e %.17e %.17e" % (particleNumber, v[1] * 1e-5, v[2] * 1e-5, v[3] * 1e-5))
            particleNumber += 1

    def convertMass(self, mass):
        """Converts mass from g/mol to kg/particle needed by the velocity calcs.

        """        
        return ((mass * 1e-3 )/ 6.022e23)

class OutTypeLAMMPSG_TruncatedBMH(OutTypeLAMMPSG):
    """Create a LAMMPSG styled output with the bmh/truncated cutoff
    style parameters.

    """    
    def __init__(self):
        pass

    def writePotentialParameters(self, numMaterials):
        """Write out the bmh/trucated style parameters.

        There are three blocks of parameters, one for items that are
        associated with a single atom (e.g. gamma_i), one that has
        pair-wise associations (e.g. r_ij), and one that is uniqe for
        three body stets (e.g. thetaA-B-C). All types contain parameters
        and cutoffs.

        The single atom items are (in order, first 3 are params, last 2 cutoffs)
            gamma_i lambda_i rho_i r_i^c 
        The pairwise atom items are (in order, first 2 params, last one cutoff)
            a_ij beta_ij r_ij
        The three body terms are (in order, a single cutoff)
            theta_i^c
        """
        # First write out the single atom coefficients, cutoffs. Check to see
        # if that particular section exists in both the Potential and Cutoff
        # sections and zero out whatever is missing. This allows the creation
        # of files that don't specify every single item which makes for more
        # sensible input files.
        self.write("\nSingle Atom Coeffs\n")
        for i in range(numMaterials):
            if config.has_section("Potential_%s" % (i)):
                gamma_i = config.getfloat("Potential_%s" % (i), "gamma")                
                lambda_i = config.getfloat("Potential_%s" % (i), "lambda")                
                rho_i = config.getfloat("Potential_%s" % (i), "rho")
            else:
                gamma_i = lambda_i = rho_i = 0.0
                
            if config.has_section("Cutoff_%s" % (i)):
                r_ic = config.getfloat("Cutoff_%s" % (i), "r")                
            else:
                r_ic = 0.0
            self.write("%.5E %.5E %.5E %.5E" % (gamma_i, lambda_i, rho_i, r_ic))

        # Now write out the pairwise atom coefficients, cutoffs.
        self.write("\nPairwise Atom Coeffs\n")
        for i in range(numMaterials):
            for j in range(numMaterials):
                if config.has_section("Potential_%s_%s" % (i , j)):
                    aij = config.getfloat("Potential_%s_%s" % (i , j), "aij")                
                    betaij = config.getfloat("Potential_%s_%s" % (i , j), "betaij")
                else:
                    aij = betaij = 0.0
                    
                if config.has_section("Cutoff_%s_%s" % (i , j)):
                    rij = config.getfloat("Cutoff_%s_%s" % (i , j), "rij")
                else:
                    rij = 0.0
                self.write("%.5E %.5E %.5E" % (aij, betaij, rij))

        # Now write out the three body cutoffs.
        self.write("\nThree Atom Coeffs\n")
        for i in range(numMaterials):
            for j in range(numMaterials):
                for k in range(numMaterials):
                    if config.has_section("Cutoff_%s_%s_%s" % (i , j, k)):
                        theta_ic = config.getfloat("Cutoff_%s_%s_%s" % (i , j, k), "theta")
                    else:
                        theta_ic = 0.0
                    self.write("%.5E" % (theta_ic))
        
class OutTypeIODYN(OutType):
    """Takes MKS pos and vel data and outputs an IODYN (cgs) compatible file.

    Outputs an IODYN compatible file for use as an OLDCON type of file. The
    units used are all cgs. The input units are MKS.
    """    
    def __init__(self):
        pass
    
    def output(self, boxsizeX, boxsizeY, boxsizeZ, pos, vel, filename="outfile", location="none", material="none"):
        self.setOutputDestination(filename, location)
        
        # Output the boxsize in cgs.
        self.write(" %.6e" % ( boxsizeX * 100.0 ))

        # Output the position in cgs
        for p in pos:
            self.write("%d %.17e %.17e %.17e" % (p[0], p[1] * 100.0, p[2] * 100.0, p[3] * 100.0))
            
        # Output the velocity in cgs
        for v in vel:
            self.write("%d %.17e %.17e %.17e" % (v[0], v[1] * 100.0, v[2] * 100.0, v[3] * 100.0))

class OutTypeMOLDY(OutType):
    """Creates a (Keith Refson's) MOLDY compatible system specification file.
    
    The format for MOLDY is broken up into the specifications for the particles
    and the potential. The particle format looks like:
    
        species-namei Ni 
        id1 x1 y1 z1 m1 q1 name1 
        id2 x2 y2 z2 m2 q2 name2 
        . 
        . 
        . ... ... ... ... ... ... 
        idn xn yn zn mn qn namen 

    """
    def __init__(self):
        pass
        
    def output(self, boxsizeX, boxsizeY, boxsizeZ, pos, vel, filename="outfile", location="none", material="none"):
        self.setOutputDestination(filename, location)
        
        # Timestamp file
        self.write("# Skew-start Configuration File Generated %s\n#" % datetime.datetime(datetime.MINYEAR, 1, 1).now())

        # Loop over the material and for each species output all generated positions for the species.
        self.write("# Positions in the following format: idn xn yn zn mn qn namen")
        species_offset = 0
        for m in material:
            # species-namei Ni
            self.write("%s %s" % (m.name, m.number))
            mass = m.mass
            charge = m.charge
            name = m.symbol
            
            # id1 x1 y1 z1 m1 q1 name1 to idn xn yn zn mn qn namen
            index = species_offset
            for i in range(m.number):
                p = pos[species_offset + i]
                self.write("%6d %9.6e %9.6e %9.6e %9.6e %9.6e %s" % (index, p[1], p[2], p[3], mass, charge, name))
                index += 1
                
            species_offset = species_offset + m.number
        self.write("end")
        
        self.writePotentialParameters(len(material))
        self.write("end")

    def writePotentialParameters(self, numMaterials):
        """Write out the MOLDY style parameters for a Buckingham potential.

        Write out the potential parameters for the Buckingham potential only.
        NOTE: The parameter names are reversed from our usual order. The first
        parameter is OUR cij though it is listed as MOLDY's aij. 
        """
        # Make sure we're doing a bmh/cutoff potential. If not, this won't
        # work.
        potential = config.get("Format", "potential")

        if potential != "bmh/cutoff":
            print "ERROR: The potential must be bmh/cutoff for this to work."
            sys.exit(1)
            
        # Output header
        self.write("# The MOLDY buckingham potential has a different order than what LAMMPS")
        self.write("# uses. The parameters are given A,B,C but A is the same thing as")
        self.write("# C in LAMMPS, B is like A in LAMMPS, and C is 1/rho.")

        self.write("Buckingham")

        # Read coefficients. Ones that are not specified are set to zero.
        for i in range(numMaterials):
            for j in range(numMaterials):
                if config.has_section("Potential_%s_%s" % (i , j)):
                    aij = config.getfloat("Potential_%s_%s" % (i , j), "aij")                
                    bij = config.getfloat("Potential_%s_%s" % (i , j), "bij")
                    cij = config.getfloat("Potential_%s_%s" % (i , j), "cij")
                elif config.has_section("Potential_%s_%s" % (j , i)):
                    aij = config.getfloat("Potential_%s_%s" % (j , i), "aij")                
                    bij = config.getfloat("Potential_%s_%s" % (j , i), "bij")
                    cij = config.getfloat("Potential_%s_%s" % (j , i), "cij")
                else:
                    aij = bij = cij = 0.0

                self.write("%.5E %.5E %.5E" % (cij, aij, bij))

    
# This code implements a version of Keith Refson's skew start method from his
# MOLDY code (not ours...). See: http://www.earth.ox.ac.uk/~keith/moldy.html
#
# NOTE: Units are MKS
def position(boxsizeX, boxsizeY, boxsizeZ, numberParticles, pos):
    """Returns a list of initial positions computed using Refson's skew start method.

    The skew start method is detailed in Refson's MOLDY code. See above. A list of
    particles is returned. Each list element is a list of X,Y,Z positions in Angstroms.
    The overall idea is to create a line which (when folded back into the simulation cube)
    goes from one corner to the other in "slices". This gives a decent fill of the space in
    the cube (except at the corners). We then go along this line and randomly assign the
    positions to different species.
    """
    # Compute the extents in the X,Y, and Z direction. Refson calls them h,k, and l
    # We round up on these. It gives us a little less interline distance but improves
    # the interparticle distance. This is our line direction vector.
    bx = math.ceil(math.pow(float(numberParticles), 2.0 / 3.0)) * boxsizeX        
    by = math.ceil(math.pow(float(numberParticles), 1.0 / 3.0)) * boxsizeY
    bz = 1 * boxsizeZ

    # The vector equation of the line (see Refson as to why we're doing a line)
    # is r' = a' + t * b' (where ' indicates a vector). We've got our direction
    # vector b'. To generate the coordinates we simply move t from 0 to 1 by an
    # amount of 1 / gTotal_particles. This moves along the vector 1/gTotal_particles
    # each time.
    for wStep in range(numberParticles):
        t = wStep / float(numberParticles)
        t += 0.5  / float(numberParticles)

        # Compute the direction vector r. Note that we're using a fixed
        # a vector and a fixed location from which to start a'. This allows
        # us to add them together to get the total offset. The total offset
        # turns out to be (0,0,0).
        rx = t * bx;
        ry = t * by;
        rz = t * bz;

        # Pick a particle at random. We keep selecting particles until there
        # is an open spot. I know it's not the most efficient but it should
        # work OK since we only generate new positions occasionally and then
        # it's not time critical.
        p = random.choice(pos)
        while(p[0] > 0):
            p = random.choice(pos)

        # Once the coordinate has been computed wrap it so that it lies within
        # the unit simulation cube.
        p.append(math.fmod(rx, float(boxsizeX)))
        p.append(math.fmod(ry, float(boxsizeY)))
        p.append(math.fmod(rz, float(boxsizeZ)))

        # Flip the sign of p[0] back to where it supposed to be (> 0) so that we
        # know it's been processed.
        p[0] *= -1

def velocity(temperature, velocity):
    """Takes a list of lists of particle masses and returns an in-place velocity distribution.

    This uses the technique described in the classic book. M. P. Allen and D. J.
    Tildesley, Computer simulation of liquids, Clarendon Press, Oxford, 1987 
    pg. 170. The formula for the velocity for the kth particle in the Cartesian
    direction i is:
     
           v_ik = sqrt(k_B * T / m_k) * R_ik where:
    
    k_B - Boltzmann's constant, T - Temperature in Kelvin, m_k mass of atom k.
    R_ik - Random number with a Gaussian distribution of unit variance & 0 mean.
    """
    k_b = 1.380658e-23 # in J/K
    
    for v in velocity:
        mass = v[1]
        base_velocity = math.sqrt(k_b * temperature / mass);

        # The first velocity overwrites the mass passed to us with the X component
        # of the velocity, the other two fill out the Y and Z components.
        v[1] = base_velocity * random.gauss(0.0, 1.0)            
        v.append(base_velocity * random.gauss(0.0, 1.0))
        v.append(base_velocity * random.gauss(0.0, 1.0))

def computeMomentum(mass,vel):
    """Takes masses, velocities and computes total momentum.
    """
    mx = my = mz = 0.0
    for v in vel:
        mass = momentum_mass[v[0]]
        mx = mx + v[1] * mass
        my = my + v[2] * mass
        mz = mz + v[3] * mass
    
    return mx, my ,mz

if __name__ == "__main__":
    # Check to see if we've passed a config file. Otherwise ask for one.
    optlist, args = getopt.getopt(sys.argv, '') 
    if len(args) != 2:
        print "usage: skewstart.py configfilename"
        sys.exit(1)

    # Open the configuration file for reading.
    config = ConfigParser.ConfigParser()
    config.readfp(open(args[1]))

    # Get the initial condition parameters
    temperature = config.getint("InitialConditions", "temperature")
    
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
            
    units = config.get("Format", "units")
    potential = config.get("Format", "potential")
    outputFilename = config.get("Format", "filename")
    outputLocation = config.get("Format", "location")

    # Figure out what kind of IO we need to do and pick the correct object to do it.    
    if units == "IODYN":
        unitsStyle = OutTypeIODYN()
    elif units == "LAMMPS":
        if potential == "bmh/truncated":
            unitsStyle = OutTypeLAMMPSG_TruncatedBMH()
        else:
            unitsStyle = OutTypeLAMMPSG()
    elif units == "MOLDY":
        unitsStyle = OutTypeMOLDY()
        
    else:
        print "ERROR: The format \"%s\" is not a recognized output format. It needs to be either \"IODYN\", \"LAMMPS\",  or \"MOLDY\"." % units
        sys.exit(1)

    # Get the material parameters. There is one section called "Material" which states how
    # many materials there are. For n materials there are n sections following. Each is called
    # Material0 .. Materialn - 1.
    numberSpecies = config.getint("Material", "number")
    material = []
    totalParticles = 0
    
    for i in range(numberSpecies):
        number = config.getint("Material%s" % (i), "number")
        totalParticles += number
        name = config.get("Material%s" % (i), "name")
        symbol = config.get("Material%s" % (i), "symbol")
        tmpMass = unitsStyle.convertMass(config.getfloat("Material%s" % (i), "mass"))
        charge = unitsStyle.convertCharge(config.getfloat("Material%s" % (i), "charge"))
        tmpMaterial = Material(i + 1, number, name, symbol, tmpMass, charge)
        material.append(tmpMaterial)

    # Compute the skew start positions. The pos value has all of the species indicies
    # made negative so that we know a position hasn't been generated for them.
    pos = []
    for m in material:
        for i in range(m.number):
            pos.append([-1 * m.index])
    position(boxsizeX, boxsizeY, boxsizeZ, totalParticles, pos)

    # Compute the skew start velocities. To do this we must first construct a list of
    # masses (in MKS) and pass that to skew start. Note that the list is actually a list
    # of single element lists.
    vel = []
    for m in material:
        for i in range(m.number):
            vel.append([m.index, m.mass])
    velocity(temperature, vel)

    momentum_mass = []
    for m in material:
        for i in range(m.number):
            momentum_mass.append(m.mass)
            
    # Show the momentum before removal.
    mv = computeMomentum(momentum_mass, vel)
    print "Before momentum removal x: %s y: %s z: %s" % mv
    
    # Remove the momentum by computing the per particle momentum that needs to
    # be removed and then adjusting the velocity profile.
    mvx = mv[0] / float(totalParticles)
    mvy = mv[1] / float(totalParticles)
    mvz = mv[2] / float(totalParticles)
    
    for v in vel:
        mass = momentum_mass[v[0]]
        v[1] = v[1] - mvx / mass
        v[2] = v[2] - mvy / mass
        v[3] = v[3] - mvz / mass
    
    # Show the momentum after removal.
    print "After momentum removal x: %s y: %s z: %s" % computeMomentum(momentum_mass, vel)

    # Output the POS and VEL files.
    unitsStyle.output(boxsizeX, boxsizeY, boxsizeZ, pos, vel, outputFilename, outputLocation, material)
    
    print "Initial conditions generated."        