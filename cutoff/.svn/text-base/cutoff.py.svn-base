#!/usr/bin/python
#
# This code computes the cutoffs needed for the coordination number program CN
# (cn.py). The format of the output file is the same as what would be present
# in the material parameter file.
import sys
import datetime
import getopt
import ConfigParser

class IntegrationError(Exception):
    def __init__(self, value):
        self.value = "*** ERROR *** Cannot find %s" % value
    def __str__(self):
        return `self.value`

def extractFloat(file, position, pattern):
    """Extracts a floting point value and returns it.
    
    Position is the position in an array if the matching line was split.
    Pattern is a regular expression of what you're trying to match.
    """
    regpat = re.compile(pattern)
    file.seek(0)
    for aline in file.readlines():
        if regpat.search(aline) != None:
            spaline = aline.split()
            return float(spaline[position])
            
    raise IntegrationError(pattern)
    
def getColumn(file, rcol, grcol):
    """Extract the column for the radius and the corresponding g(r) from a
    LAMMPS formatted run column.
    
    Rewind the file and then loop through it appending the correct column
    data on to the returned array.
    """
    ret_r = []
    ret_gr = []
    
    file.seek(0)
    dummyLine = file.readline()
    dummyLine = file.readline()
    try:
        for aLine in file:
            splitLine = aLine.split()
            ret_r.append(float(splitLine[rcol]))
            ret_gr.append(float(splitLine[grcol]))
    except:
        print "*** ERROR *** Column contains a non-numerical entity. The data format is rigid and must have EXACTLY two lines before the first line of data and NO lines after the last line of data."
        sys.exit(1)
        
    return ret_r, ret_gr
    
if __name__ == "__main__":
    # Check to see if we've passed a config file. Otherwise ask for one.
    optlist, args = getopt.getopt(sys.argv, '') 
    if len(args) != 2:
        print "usage: cutoff.py materialfilename"
        print "    where:"
        print "        materialfilename - Material parameters."
        sys.exit(1)

    # Open the configuration file for reading.
    config = ConfigParser.ConfigParser()
    config.readfp(open(args[1]))

    # Open the I/O files.
    rdffile = file(config.get("Cutoff", "input_filename"), "r")
    
    # If the optional output filename is present then output to this file
    # otherwise just send it out the usual way (through stdout).
    try:
        outfilename = config.get("Cutoff", "output_filename")
        outfile = file(outfilename, "w")
    except ConfigParser.NoOptionError:
        outfile = sys.stdout

    # Figure out how many materials there are, their names, and compute what 
    # you're going to need to know.
    numberSpecies = config.getint("Material", "number")
    
    # Get the names of the material.
    materialSymbols = []
    for materialNumber in range(numberSpecies):
        materialSymbols.append(config.get("Material%d" % materialNumber, "symbol"))

    # Create cutoff array
    cutoff = [ [ 0.0  for i in range (numberSpecies) ] for j in range(numberSpecies) ]
    
    # Now we know what we need to do. There should be numberSpecies^2 combinations.
    # Iterate over each one separately.
    rcol = config.getint("Cutoff", "r_col")

    for center in range(numberSpecies):
        for outlier in range(numberSpecies):
            grcol = config.getint("Cutoff", "g_%s_%s_col" % (center, outlier))
            r, gr = getColumn(rdffile, rcol, grcol)
            
            # Find the first minima. For g(r) functions the first minima always follows
            # the first maxima. So find the first maxima since that is easy.
            grMaxIndex = 0
            for wgr in range(len(gr)):
                if gr[wgr] > gr[grMaxIndex]:
                    grMaxIndex = wgr
                    
            # Starting from the maximum point search for the first minima. This is
            # defined as a point which is the lowest point for five iterations.
            grMinIndex = grMaxIndex
            grMinDuration = 0
            grMin = gr[grMaxIndex]
            
            for wgr in range(grMaxIndex, len(gr)):
                if gr[wgr] < gr[grMinIndex]:
                    grMinIndex = wgr
                    grMinDuration = 0
                else:
                    if grMinDuration == 5:
                        break
                    else:
                        grMinDuration += 1
                        
            cutoff[center][outlier] = float(r[grMinIndex])

    # Dump cutoffs. The format needs to be CN (program cn.py) compatible. The
    # format is:
    #       # Mg surrounded by Mg
    #       [CN_0_0]
    #       cutoff: 4.15
    outfile.write("# Cutoff File Generated %s\n#\n" % datetime.datetime(datetime.MINYEAR, 1, 1).now())
    outfile.write("# Coordination Number (CN) cutoffs .\n")
    outfile.write("# Specify what the cutoffs are for each material pair (like the potential parameters). There should\n")
    outfile.write("# be one cutoff for each pair. Units are in Angstroms. Note that for Silicon surrounded by Oxygen\n")
    outfile.write("# (O about Si) you list the center atom FIRST and the surrounding atom LAST (e.g. CN_2_3). Items\n")
    outfile.write("# generic to all cutoffs are under the [CN] section. Note that cutoffs < 0.0 indicate that no \n")
    outfile.write("# coordination statistics should be calculated.\n")
    outfile.write("#\n")

    for i in range(numberSpecies):
        for j in range(numberSpecies):
            outfile.write("\n")
            outfile.write("# %s surrounded by %s\n" % (materialSymbols[i], materialSymbols[j]))
            outfile.write("[CN_%s_%s]\n" % (i, j))
            outfile.write("cutoff: %4.2f\n" % (cutoff[i][j]))
