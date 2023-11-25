#!/usr/bin/python
#
# Coordination Statistics cn.py
#
# Code computes the coordination number of one species about another.
import ConfigParser, os, os.path, getopt, sys, math, random, re, time
import numarray

cn = []         # Contains the coordination statistics.
tAboutO = []    # Coordination statistics for the Tetrahedral formers vs. Oxygen.
oAboutT = []    # Coordination statistics for the Oxygen vs. Tetrahedral formers
material = []   # Material used in simulation.
cutoff = []     # Cutoff for coordination stats.
boxsize = 0.0   # Size of simulation cube.
doit = []
numberSpecies = 0
T = []          # Array of tetrahedral formers


maxCountCN = 0.0   # How many runs are we looking at to get CN
maxSteps = 0.0      # The number of steps used in the CN calcs.
countCN = 0.0      # How many runs so far
countSteps = 0.0        # How many pos files read so far.
maxCoordination = 0     # Maximum coordination number expected.

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
        return "%s particles of %s. Mass: " % ( self.number, self.symbol, self.mass )

def convertMass(mass):
    """Converts mass from g/mol to kg/particle needed by the velocity calcs.
    """        
    return ((mass * 1e-3 )/ 6.022e23)

def scale(boxsize, px, py, pz):
    """Scales the positions to physical coordinates.
    """
    halfBoxsize = boxsize / 2.0
    px = px * boxsize - halfBoxsize
    py = py * boxsize - halfBoxsize
    pz = pz * boxsize - halfBoxsize
    return (px, py, pz)
    
def computeCoordinationStats(boxsize, specie, specieStart, specieEnd, doit, cutoff, tFormers, oxygenSpecie):
    """Computes the coordination statistics.
    """
    global px, py, pz, dx, dy, dz, rsq
    half_boxsize = boxsize / 2.0
    negative_half_boxsize = -half_boxsize
    
    percentComplete = (countCN / maxCountCN) * 100
    print "Computing coordination stats run %s/%s (%d%% complete) on configuration %s" % (int(countCN + 1), int(maxCountCN), percentComplete,int(countSteps))
    
    numParticles = len(px)
    
    for center in range(numParticles):
        cpx = px[center]
        cpy = py[center]
        cpz = pz[center]
        centerSpecie = specie[center]
        
        # Skip any species that we don't need to compute the stats for.
        if not doit[centerSpecie]:
            continue

        # Try to speed up the center calculation by using numarrays and doing
        # things in a more numarray oriented approach. This should result in
        # improved processing times at the expense of code clarity.
        dx = px - cpx
        dy = py - cpy
        dz = pz - cpz

        # Compensate for periodic boundary condtions.
        
        # Test the case where dx > boxsize and if that's true it subtracts the boxsize from dx to keep it in the simulation cube.
        # It also needs to check where dx < -boxsize (negative equivalent of the first case). In that case boxsize is added to pull
        # it into the simulation cube.
        numarray.where(dx > half_boxsize, dx - boxsize, dx, dx)
        numarray.where(dx < negative_half_boxsize, dx + boxsize, dx, dx)
        
        # Same for dy and dz
        numarray.where(dy > half_boxsize, dy - boxsize, dy, dy)
        numarray.where(dy < negative_half_boxsize, dy + boxsize, dy, dy)
        numarray.where(dz > half_boxsize, dz - boxsize, dz, dz)
        numarray.where(dz < negative_half_boxsize, dz + boxsize, dz, dz)

        # Now all deltas are within the simulation cube we need to compute the distance squared from the center to all outliers.
        rsq = dx * dx + dy * dy + dz * dz
        
        # For each outlier specie extract the portion of rsq that applies and find the highest coordination.
        for outlierSpecie in range(numberSpecies):
            rsqSpecie = rsq[specieStart[outlierSpecie]:specieEnd[outlierSpecie]]
            scwSpecie = numarray.sum(numarray.logical_and(rsqSpecie > 0.01, rsqSpecie < cutoff[centerSpecie][outlierSpecie]))
            if scwSpecie == 0:
                continue
            # Accumulate the coordination value for each center species, and outlier species to cn.
            try:
                cn[centerSpecie][outlierSpecie][scwSpecie] += 1
            except:
                print "The maximum expected coordination was %d and we have found a coordination of %d. Please raise your coordination number." % (maxCoordination - 1, scwSpecie)
                sys.exit(1)
            
            # Check to see if the outlier specie is Oxygen and the center species are a tetrahedral former.
            if len(tFormers) >= 2 and outlierSpecie == oxygenSpecie and sum(map(lambda x: x == centerSpecie, tFormers)) > 0:
                try:
                    oAboutT[scwSpecie] += 1
                except:
                    print "The maximum expected coordination for O about T was %d and we have found a coordination of %d. Please raise your coordination number." % (maxCoordination - 1, scwSpecie)
                    sys.exit(1)

        # Compute T around O. Since we do everything in a vector manner we need to find the coordination of each member of the tFormers array about O
        # and then combine them to get the total coordination of the T. Note that this is very different from the way that the other coordination
        # stats are computed.
        if centerSpecie == oxygenSpecie:
            coordinationNumber = 0
            for tetrahedralFormerSpecieIndex in range(len(tFormers)):
                outlierSpecie = tFormers[tetrahedralFormerSpecieIndex]
                rsqSpecie = rsq[specieStart[outlierSpecie]:specieEnd[outlierSpecie]]
                coordinationNumber += numarray.sum(numarray.logical_and(rsqSpecie > 0.01, rsqSpecie < cutoff[centerSpecie][outlierSpecie]))

            try:
                tAboutO[coordinationNumber] += 1
            except:
                print "The maximum expected coordination for T about O was %d and we have found a coordination of %d. Please raise your coordination number." % (maxCoordination - 1, scwSpecie)
                sys.exit(1)

if __name__ == "__main__":
    # Check to see if we've passed a config file. Otherwise ask for one.
    optlist, args = getopt.getopt(sys.argv, '') 
    if len(args) != 2:
        print "usage: cn.py configfilename"
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

    boxsize = config.getfloat("InitialConditions", "boxsize")
    halfboxsize = boxsize / 2.0
    
    lastStartIndex = 0
    specieStart = [0 for i in range(numberSpecies)]
    specieEnd = [0 for i in range(numberSpecies)]
    
    for i in range(numberSpecies):
        number = config.getint("Material%s" % (i), "number")
        nAtom += number
        specieStart[i] = lastStartIndex
        specieEnd[i] = specieStart[i] + number - 1
        lastStartIndex = specieEnd[i] + 1
        name = config.get("Material%s" % (i), "symbol")
        mass = convertMass(config.getfloat("Material%s" % (i), "mass"))
        charge = config.getfloat("Material%s" % (i), "charge")
        tmpMaterial = Material(i + 1, number, name, mass, charge)
        material.append(tmpMaterial)

    # Find the O specie
    oxygenSpecies = -1
    for i in range(numberSpecies):
        if material[i].symbol == "O":
            oxygenSpecies = i
    if oxygenSpecies == -1:
        print "Oxygen was not found in the material file. Error."
        sys.exit(1)
        
    # Get Coordination number statistics parameters. We need to add one
    # to the maxCoordination number because we go from 0..n-1 when we
    # want 0..n
    # Number of steps used to compute stats.
    try:
        maxCoordination = config.getint("CN", "max_coordination")
    except:
        print 'There must be a "max_coordination" parameter present in the input file.'
        sys.exit(1)
    maxCoordination = maxCoordination + 1

    # Look to see if you need to skip any of the initial runs. This is the count
    # of the run number divorced from any sense of the timestep.
    try:
        startAtStep = config.getint("CN", "start_at_step")
    except ConfigParser.NoOptionError:
        startAtStep = 0

    # Number of steps used to compute stats.
    try:
        maxSteps = config.getfloat("CN", "steps")
    except:
        print 'There must be a "steps" parameter present in the input file.'
        sys.exit(1)

    # You can skip a number of data points within a data file.
    try:
        processEvery = config.getint("CN", "process_every")
    except ConfigParser.NoOptionError:
        processEvery = 1

    # maxCountRDF needs to be divided by the processEvery number since
    # process every implies a skip in the step. If you want 100 steps
    # but to process_every 2 then you only really have 50 steps.
    maxCountCN = math.floor((maxSteps - float(startAtStep)) / float(processEvery))

    cn = [ [ [ 0 for wbin in range(maxCoordination) ] for col in range(numberSpecies) ] for row in range(numberSpecies) ]
    tAboutO = [ 0 for wbin in range(maxCoordination) ]
    oAboutT = [ 0 for wbin in range(maxCoordination) ]
    cutoff = [ [ 0 for col in range(numberSpecies) ] for row in range(numberSpecies) ]
    doit = [ 0 for wbin in range(numberSpecies) ]

    # Get the Tetrahedral formers if present and load the T array with them. If there are less than 2 present then silently fail.
    for wSpecie in range(numberSpecies):
        if config.has_option("CN", "T_%s" % wSpecie):
            T.append(config.getint("CN", "T_%s" % wSpecie))

    if len(T) >= 2:
        computeTAboutO = True
    else:
        computeTAboutO = False

    # Cutoffs can come from either inside the configuration file or from an
    # external file. Get the source of the cutoffs and then read them.
    try:
        cutoff_filename = config.get("CN", "cutoff_filename")
        cutoff_file = ConfigParser.ConfigParser()
        cutoff_file.readfp(open(cutoff_filename))
    except ConfigParser.NoOptionError:
        cutoff_file = config

    # Get the cutoff range and then examine the range for items we can skip. Something can
    # be skipped if ALL its columnar items are -1.0. We put the cutoff**2 in the cutoff
    # array for speed of checking distances.
    for center in range(numberSpecies):
        for outlier in range(numberSpecies):
            cutoff[center][outlier] = cutoff_file.getfloat("CN_%s_%s" % (center, outlier), "cutoff")
            if cutoff[center][outlier] > 0.0:
                cutoff[center][outlier] = cutoff[center][outlier] * cutoff[center][outlier]

    for center in range(numberSpecies):
        sumcutoffs = 0
        for outlier in range(numberSpecies):
            if cutoff[center][outlier] < 0.0:
                sumcutoffs += 1
        if sumcutoffs == numberSpecies:
            doit[center] = 0
        else:
            doit[center] = 1

    # Look to see if you need to skip any of the initial runs. This is the count
    # of the run number divorced from any sense of the timestep.
    try:
        startAtStep = config.getint("CN", "start_at_step")
    except ConfigParser.NoOptionError:
        startAtStep = 0

    # Setup the input file
    posfile = file(config.get("CN", "position_filename"), "r")
    
    # Parse the position file. "ITEM: ATOMS" indicates that the next
    # nAtom lines are atom positions. This can terminate with
    # the presence of "ITEM: TIMESTEP" or termination of the loop.
    curLine = "nothing"
    pos = ["" for x in range(nAtom)]
    
    # Create parallel arrays for testing new approach.
    px = numarray.zeros(nAtom, numarray.Float)
    py = numarray.zeros(nAtom, numarray.Float)
    pz = numarray.zeros(nAtom, numarray.Float)
    specie = numarray.zeros(nAtom, numarray.Int)
    dx = numarray.zeros(nAtom, numarray.Float)
    dy = numarray.zeros(nAtom, numarray.Float)
    dz = numarray.zeros(nAtom, numarray.Float)
    rsq = numarray.zeros(nAtom, numarray.Float)
    
    for line in posfile:
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
                        (px, py, pz) = scale(boxsize, px, py, pz)
                        computeCoordinationStats(boxsize, specie, specieStart, specieEnd, doit, cutoff, T, oxygenSpecies)
                        countCN += 1.0
                        if countCN >= maxCountCN:
                            break
                    if countSteps >= maxSteps:
                        break
                countSteps += 1.0
            curLine = "nothing"

        if curLine == "getAtoms":
            lineAr = line.split()
            curAtomIndex = int(lineAr[0]) - 1
            specie[curAtomIndex] = int(lineAr[1]) - 1
            px[curAtomIndex] = float(lineAr[2])
            py[curAtomIndex] = float(lineAr[3])
            pz[curAtomIndex] = float(lineAr[4])
            
        # These set the stage for the NEXT line so they must come after the processing.            
        if re.search("ITEM: ATOMS", line):
            curLine = "getAtoms"
    else:
        # Be sure to process the last block of atoms.
        if curLine == "getAtoms" and countCN < maxCountCN and countSteps < maxSteps:
            countSteps += 1.0
            countCN += 1.0
            computeCoordinationStats(pos, doit, cutoff)

    # Compute the number of 0-fold items. This is computed by subtracting the
    # sum of all the other coordinations from the total to be distributed. The
    # total to be distributed is the number of each specie * countCN.
    for center in range(numberSpecies):
        maxToDistribute = material[center].number * countCN
        
        print "material[center].number %s countCN %s" % (material[center].number, countCN)
        
        for outlier in range(numberSpecies):
            if center == outlier:
                continue

            totalDistributed = 0            
            for nfold in range(maxCoordination):
                totalDistributed += cn[center][outlier][nfold]    

            print "center: %s outlier %s maxToD %s totalDist %s" % ( center, outlier, maxToDistribute, totalDistributed)

            # If 0 then nothing distributed. Otherwise compute number of 0-folds
            if totalDistributed == 0:
                continue
            else:
                cn[center][outlier][0] = maxToDistribute - totalDistributed   
    
    # Compute the number of 0-fold items for the O about T.
    maxToDistribute = 0
    totalDistributed = 0            
    for i in range(len(T)):
        center = T[i]
        maxToDistribute += material[center].number * countCN
        
    for nfold in range(maxCoordination):
        totalDistributed += oAboutT[nfold]    

    # If 0 then nothing distributed. Otherwise compute number of 0-folds
    if totalDistributed == 0:
        oAboutT[0] = 0
    else:
        oAboutT[0] = maxToDistribute - totalDistributed   
    
    print "O about T: maxToDistribute %s totalDistributed %s" % (maxToDistribute, totalDistributed)
    
    # Compute number of 0-folds from T about O. The max to distribute is the number of O centers times the number of tetrahedral
    # formers since each set of Os is grouped with each specie.
    maxToDistribute = 0
    totalDistributed = 0            
    center = oxygenSpecies
    maxToDistribute = material[center].number * countCN
    
    for nfold in range(maxCoordination):
        totalDistributed += tAboutO[nfold]    

    # If 0 then nothing distributed. Otherwise compute number of 0-folds
    if totalDistributed == 0:
        tAboutO[0] = 0
    else:
        tAboutO[0] = maxToDistribute - totalDistributed   
    
    print "T about O: maxToDistribute %s totalDistributed %s" % (maxToDistribute, totalDistributed)
    
    # Output Coordination statistics.
    print "\nGenerating coordination statistics file"

    # Get the filename (if present).    
    try:
        outfilename = config.get("CN", "output_filename")
    except:
        outfilename = "cnout.txt"
        
    cnfile = file(outfilename, "w")
    
    # Output the coordination statistics for all single species data.
    for center in range(numberSpecies):
        for outlier in range(numberSpecies):
            cnfile.write("\n%s around %s\n" % ( material[outlier].symbol, material[center].symbol))
            # Output header.
            for n in range(maxCoordination):
                cnfile.write("%d\t" % (n) )
            cnfile.write("\n\n")
            
            # Output coordination data.
            for n in range(maxCoordination):
                cnfile.write("%d\t" % (cn[center][outlier][n]) )
            cnfile.write("\n\n")
                
    # Output the coordination statistics for the T about O data point (if needed)
    if len(T) > 0:
            cnfile.write("\nT around O\n")
            # Output header.
            for n in range(maxCoordination):
                cnfile.write("%d\t" % (n) )
            cnfile.write("\n\n")
            
            # Output coordination data.
            for n in range(maxCoordination):
                cnfile.write("%d\t" % (tAboutO[n]) )
            cnfile.write("\n\n")

    # Output the coordination statistics for the O about T data point (if needed)
    if len(T) > 0:
            cnfile.write("\nO around T\n")
            # Output header.
            for n in range(maxCoordination):
                cnfile.write("%d\t" % (n) )
            cnfile.write("\n\n")
            
            # Output coordination data.
            for n in range(maxCoordination):
                cnfile.write("%d\t" % (oAboutT[n]) )
            cnfile.write("\n\n")

    cnfile.flush()
    cnfile.close()

    print "Coordination Statistics Generated"
