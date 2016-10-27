
from geomclass import GeomClass
import sys

#--- Set up path of IO files
ioPath = None
if len(sys.argv) > 1:
    ioPath = sys.argv[1]
    if ioPath[-1] != '/': ioPath += '/'

geom = GeomClass.createFromAllInputs(ioPath)
geom.displaceCarts(updateSelf=True)

geom.printAllRedundIntCoords()
