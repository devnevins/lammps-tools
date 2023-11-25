#!/usr/bin/python
import ConfigParser, os, os.path, getopt, sys, math, random, re

rdf = []        # Contains the Radial Distribution Function histogram.
material = []   # Material used in simulation.
boxsize = 0.0   # Size of simulation cube.

maxCountRDF = 0.0   # How many runs are we looking at to get g(r)
maxSteps = 0.0      # The number of steps used in the RDF calcs.
countRDF = 0.0      # How many runs so far
countSteps = 0.0        # How many pos files read so far.
processedSteps = 0.0    # How many pos files processed so far.

# Defines what a particle is
class Particle:
    def __init__(self, specie, x, y, z):
        self.specie = specie
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def scale(self):
        """Scale the particle by the input boxsize.

        This is meant to allow the particles (0.0 < x < 1.0 in LAMMPSG)
        to have actual dimensions.
        """
        self.x = self.x * boxsize - (boxsize / 2.0)
        self.y = self.y * boxsize - (boxsize / 2.0)
        self.z = self.z * boxsize - (boxsize / 2.0)

    def __str__(self):
        """What a particles prints like.
        """
        return "%5Ef, %5Ef, %5Ef" % (self.x, self.y, self.z)

    def signR(self, x, y):
        """Returns a sign adjusted value of x depending on y.
        """
        if y >= 0.0:
            return x
        else:
            return -x
        
    def distanceToSquared(self, dst):
        """Computes the adjusted distance from self to dst.

        The distance computed is adjusted so that it always stays within the primary box.        
        """
        half_boxsize = boxsize / 2.0
        
        dx = dst.x - self.x
        if math.fabs(dx) > half_boxsize:
            dx = dx - self.signR(boxsize, dx)
            
        dy = dst.y - self.y
        if math.fabs(dy) > half_boxsize:
            dy = dy - self.signR(boxsize, dy)

        dz = dst.z - self.z
        if math.fabs(dz) > half_boxsize:
            dz = dz - self.signR(boxsize, dz)

        return(dx * dx + dy * dy + dz * dz)

# Material class used internally by skewstart. This didn't seem to have much utility outside
# of this code so I included it.
class Material:
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
        return "%s particles of %s. Mass: " % (self.number, self.symbol, self.mass)



def convertMass(mass):
    """Converts mass from g/mol to kg/particle needed by the velocity calcs.
    """        
    return ((mass * 1e-3)/ 6.022e23)

def computeRDF(pos, sizeHistRdf, rangeRdf):
    """Computes the Radial Distribution Function (RDF) g(r).

    Computes the RDF (g(r)). See Rapaport, D.C. (1995) The Art of Molecular
    Dynamics Simulation, pg.87.
    """
    percentComplete = (countRDF / maxCountRDF) * 100
    print "Computing RDF run %s/%s (%d%%) on configuration %s" % (countRDF, maxCountRDF, percentComplete, countSteps)

    rrRange = rangeRdf * rangeRdf
    deltaR = float(rangeRdf) / float(sizeHistRdf)
    
    for j1index in range(0, len(pos) - 1):
        for j2index in range(j1index + 1, len(pos)):
            j1 = pos[j1index]
            j2 = pos[j2index]
            rr = j1.distanceToSquared(j2)
            if rr < rrRange:
                n = int(math.sqrt(rr) / deltaR)
                #print "j1 %s j2 %s n %s" % (j1.specie, j2.specie, n)
                rdf[j1.specie][j2.specie][n] += 1            
                rdf[j2.specie][j1.specie][n] += 1            

if __name__ == "__main__":
    # Check to see if we've passed a config file. Otherwise ask for one.
    optlist, args = getopt.getopt(sys.argv, '') 
    if len(args) != 3:
        print "usage: rdf.py configfilename posfilename"
        sys.exit(1)

    # Open the configuration file for reading.
    config = ConfigParser.ConfigParser()
    config.readfp(open(args[1]))

    # Get the material parameters. There is one section called "Material" which states how
    # many materials there are. For n materials there are n sections following. Each is called
    # Material0 .. Materialn - 1.
    numberSpecies = config.getint("Material", "number")
    material = []
    nAtom = 0

    # Look to see if you need to skip any of the initial runs. This is the count
    # of the run number divorced from any sense of the timestep.
    try:
        startAtStep = config.getint("RDF", "start_at_step")
    except ConfigParser.NoOptionError:
        startAtStep = 0

    boxsize = config.getfloat("InitialConditions", "boxsize")
    halfboxsize = boxsize / 2.0
    
    for i in range(numberSpecies):
        number = config.getint("Material%s" % (i), "number")
        nAtom += number
        name = config.get("Material%s" % (i), "symbol")
        mass = convertMass(config.getfloat("Material%s" % (i), "mass"))
        charge = config.getfloat("Material%s" % (i), "charge")
        tmpMaterial = Material(i + 1, number, name, mass, charge)
        material.append(tmpMaterial)

    # Get RDF parameters
    sizeHistRdf = config.getint("RDF", "bins")
    rangeRdf = config.getfloat("RDF", "extent")
    rdf = [ [ [ 0 for wbin in range(sizeHistRdf) ] for col in range(numberSpecies) ] for row in range(numberSpecies) ]

    # Number of steps used to compute stats.
    try:
        maxSteps = config.getfloat("RDF", "steps")
    except:
        print 'There must be a "steps" parameter present in the input file.'
        sys.exit(1)

    # You can skip a number of data points within a data file.
    try:
        processEvery = config.getint("RDF", "process_every")
    except ConfigParser.NoOptionError:
        processEvery = 1

    # maxCountRDF needs to be divided by the processEvery number since
    # process every implies a skip in the step. If you want 100 steps
    # but to process_every 2 then you only really have 50 steps.
    maxCountRDF = config.getfloat("RDF", "steps")
    maxCountRDF = maxCountRDF / float(processEvery)
    
    # Parse the position file. "ITEM: ATOMS" indicates that the next
    # nAtom lines are atom positions. This can terminate with
    # the presence of "ITEM: TIMESTEP" or termination of the loop.
    curLine = "nothing"
    pos = ["" for x in range(nAtom)]
    for line in file(args[2]):
        # Checks to make sure the line isn't the beginning of a new
        # block.
        if re.search("ITEM: TIMESTEP", line):
            if curLine == "getAtoms":
                # print "Read configuration number %s" % (countSteps)
                # Ensure that we are supposed to process the step. We may want
                # to skip some of the beginning of the data file.
                if countSteps >= startAtStep:
                    # We've reached the end of the block, process block.
                    if (countSteps - startAtStep) % processEvery == 0:
                        computeRDF(pos, sizeHistRdf, rangeRdf)
                        countRDF += 1.0
                        if countRDF >= maxCountRDF:
                            break
                    if countSteps >= maxSteps:
                        break
                countSteps += 1.0
            curLine = "nothing"

        if curLine == "getAtoms":
            lineAr = line.split()
            pos[int(lineAr[0]) - 1] = Particle(int(lineAr[1]) - 1, lineAr[2], lineAr[3], lineAr[4])
            pos[int(lineAr[0]) - 1].scale()
            
        # These set the stage for the NEXT line so they must come after the processing.            
        if re.search("ITEM: ATOMS", line):
            curLine = "getAtoms"

    # Be sure to process the last block of atoms.
    if curLine == "getAtoms":
        computeRDF(pos, sizeHistRdf, rangeRdf)
        countRDF += 1.0

    # Normalize the RDF by dividing by the shell volume used to create each
    # bin.
    deltaR = rangeRdf / sizeHistRdf

    for j1 in range(numberSpecies):
        for j2 in range(numberSpecies):
            normFac = boxsize**3 / (2.0 * math.pi * deltaR**3 * nAtom**2 * countRDF)
            for n in range(0, sizeHistRdf):
                rdf[j1][j2][n] = rdf[j1][j2][n] * normFac / (float(n) - 0.5)**2
                
    # Get the output filename (if present).    
    try:
        outfilename = config.get("RDF", "filename")
    except:
        outfilename = "rdfout.txt"
        
    # Output RDF [g(r)]
    print "Generating RDF [g(r)] file. Writing to ""%s""." % outfilename

    rdffile = file(outfilename, "w")

    # Output the data in a columnar format to ease cutoff processing and plotting.
    # Header is output first. Note that the radius column is duplicated for each pair.
    for si in range(numberSpecies):
        for sj in range(numberSpecies):
            rdffile.write("r %s_around_%s " % (material[sj].symbol, material[si].symbol))
    rdffile.write("\n")
                
    # Followed by the data.
    for n in range(0, sizeHistRdf):
        rBin = (float(n) + 0.5) * rangeRdf / sizeHistRdf
        for si in range(numberSpecies):
            for sj in range(numberSpecies):
                rdffile.write("%3.3f %3.4E " % (rBin, rdf[si][sj][n]))
        rdffile.write("\n")
        
    rdffile.flush()
    rdffile.close()

    print "Statistics Generated"
