import sys, getopt
import re

from ebdtFunctions import *

try:
	opts, args = getopt.getopt(sys.argv[1:],"p:r:a:c:")
except getopt.GetoptError:
	sys.exit()

inhibitionThreshold = float(0.5)
ratioThreshold = float(0.5)
probabilityThreshold = float(0.75)
cellLinesFiles = ["MCF7.xlsm", "HL60.xlsm", "NTERA2.xlsm"]

for opt, arg in opts:
	if opt == "-a":
		inhibitionThreshold = float(arg)
	if opt == "-c":
		cellLinesFiles = [cellLineFile for cellLineFile in arg.split(',')]
	if opt == "-r":
		ratioThreshold = float(arg)
	if opt == "-p":
		probabilityThreshold = float(arg)

for cellLineFile in cellLinesFiles:
	if not re.match(r"^.+\.+xl.+$", cellLineFile):
		print(cellLineFile+" is not a correcly formatted file name. File name must match regex ^.+\\.+xl.+$")
		quit()

print("Inhibition threshold = %.2f" % inhibitionThreshold)
print("Ratio threshold = %.2f" % ratioThreshold)
print("Probability threshold = %.2f" % probabilityThreshold)
print("Cell lines = %s" % cellLinesFiles)

GetExpectancyOfBeingDownstreamTarget(inhibitionThreshold, ratioThreshold, probabilityThreshold, cellLinesFiles)