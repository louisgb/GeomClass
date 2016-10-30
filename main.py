
from geomclass import GeomClass
import sys

#-- Set up path of IO files
ioPath = None
if len(sys.argv) > 1:
    ioPath = sys.argv[1]
    if ioPath[-1] != '/': ioPath += '/'

#-- Create GeomClass object from input files
geom = GeomClass.createFromAllInputs(ioPath)
#-- Print initial geometry
print '=== Initial Cartesian geometry ==='
geom.printCarts()
#-- Print initial internal coordinates
print '\n=== Initial internal coordinates ==='
geom.printAllRedundIntCoords()
#-- Displace Cartesian coordinates
print '\n=== Algorithm statistics ==='
geom.displaceCarts(updateSelf=True, printOutput=True)
#-- Print final geometry
print '\n=== Final internal coordinates ==='
geom.printAllRedundIntCoords()

