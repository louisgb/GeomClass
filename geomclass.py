import numpy as np
import os
import sys
import geomclass_aux as aux
import copy

class GeomClass(object):
    """A class for molecular geometry.

    """
    def __init__(self, ioPath=None):
        """Constructor.

        ioPath: path containing all input files.
        """
        #-- atomic symbols
        self.atoms = None
        #-- atomic masses
        self.masses = None
        #-- cartesian coordinates
        self.carts = None
        #-- definition of redundant int. coord.
        #--   as a list of lists, each sub-list
        #--   being atomic indexes (1-based)
        #-- bl, ba, to, ob, od =
        #-- bond length, bond angle, torsion, oop bend, oop dist.
        self.iblrs = None
        self.ibars = None
        self.itors = None
        self.iobrs = None
        self.iodrs = None
        #-- # of nonredund. int. coord.
        self.nblnr = None
        self.nbanr = None
        self.ntonr = None
        self.nobnr = None
        self.nodnr = None
        self.nintnr = None
        #-- A matrix = dCart/dNRIC [3N * 3N-6]
        #-- B matrix = dNRIC/dCart [3N-6 * 3N]
        #-- T matrix = dNRIC/dRIC  [3N-6 * nRIC]
        #--   where N = non-; R = redundant;
        #--        IC = internal coordinate
        self.amat = None
        self.bmat = None
        self.tmat = None
        #-- path of code
        self.codePath = os.path.dirname(__file__)
        if not self.codePath:
            self.codePath = '.'
        self.codePath += '/'
        #-- path of IO files; default to current path
        if ioPath is None:
            self.ioPath = os.getcwd()+'/'
        else:
            self.ioPath = ioPath
        #-- full path of individual files
        self.filenameCart = None
        self.filenameR = None
        self.filenameNR = None
        self.filenameDisplNR = None

        
    @classmethod
    def createFromAllInputs(cls, ioPath=None):
        """Instantiate by reading in all input files.

        ioPath: path containing all input files.
        """
        newObject = cls(ioPath)
        newObject.readAllInputs()
        return newObject


    def readAllInputs(self):
        """Wrapper to read in all input files. 

        """
        self.readCarts()
        self.readRedundIntCoords()
        self.readNonRedundIntCoords()
        self.calcABmat()
        

    def setIOpath(self, path):
        """Set path of IO files.

        """
        self.ioPath = path
        if self.ioPath[-1] != '/': self.ioPath += '/'


    def readCarts(self, filename=None):
        """Read in Cartesian coordinates from a file.

        Default filename = self.ioPath+'cart0.txt'
        Note: self.carts will be a 1-D array with 3*N elements.
              To get a N-by-3 matrix, use np.reshape(self.carts, (N, 3))
        """
        if filename is None:
            filename = self.filenameCart
        if filename is None:
            filename = self.ioPath+'cart0.txt'
            self.filenameCart = filename
        #-- Read in initial Cartesian coordinates
        with open(filename, 'r') as fcart:
            line = aux.skipComment(fcart, '#').split()
            natom = int(line[0])
            atoms = []; masses = []; cartslist = []
            for i in xrange(natom):
                line = aux.skipComment(fcart, '#').split()
                atoms.append(line[0])
                masses.append(float(line[1]))
                for j in xrange(3):
                    #-- Setting the zero coord. to 1e-10 
                    #--  to prevent some numerical problems
                    if abs(float(line[j+2]))<1e-10:
                        cartslist.append(1e-10)
                    else:
                        cartslist.append(float(line[j+2]))
            self.carts = np.array(cartslist)
            self.atoms = atoms
            self.masses = masses
                

    def readRedundIntCoords(self, filename=None):
        """Read in definition of redund. int. coord. from a file.

        Default filename = self.ioPath+'rdef0.txt'
        """
        if filename is None:
            filename = self.filenameR
        if filename is None:
            filename = self.ioPath+'rdef0.txt'
            self.filenameR = filename

        iblrs, ibars, itors, iobrs, iodrs = [], [], [], [], []
        with open(filename, 'r') as f:
            #-- bond lengths
            line = aux.skipComment(f, '#')
            nblr = int(line)
            for i in xrange(nblr):
                line = aux.skipComment(f, '#')
                iblrs.append(line.split())
            #-- bends
            line = aux.skipComment(f, '#')
            nbar = int(line)
            for i in xrange(nbar):
                line = aux.skipComment(f, '#')
                ibars.append(line.split())
            #-- torsions
            line = aux.skipComment(f, '#')
            ntor = int(line)
            for i in xrange(ntor):
                line = aux.skipComment(f, '#')
                itors.append(line.split())
            #-- oop bends
            line = aux.skipComment(f, '#')
            nobr = int(line)
            for i in xrange(nobr):
                line = aux.skipComment(f, '#')
                iobrs.append(line.split())
            #-- oop distances
            line = aux.skipComment(f, '#')
            nodr = int(line)
            for i in xrange(nodr):
                line = aux.skipComment(f, '#')
                iodrs.append(line.split())
            iblrs = np.array(iblrs, dtype=int)
            ibars = np.array(ibars, dtype=int)
            itors = np.array(itors, dtype=int)
            iobrs = np.array(iobrs, dtype=int)
            iodrs = np.array(iodrs, dtype=int)
        #-- account for the array index convention
        self.iblrs = iblrs-1
        self.ibars = ibars-1
        self.itors = itors-1
        self.iobrs = iobrs-1
        self.iodrs = iodrs-1


    def readNonRedundIntCoords(self, filename=None):
        """Read in definition of nonredund. int. coord. from a file.

        Also set T matrix and prepare input for abmat.exe
          for calculating A and B matrices. 
        Default filename = self.ioPath+'nrdef0.txt'
        """
        if filename is None:
            filename = self.filenameNR
        if filename is None:
            filename = self.ioPath+'nrdef0.txt'
            self.filenameNR = filename

        nintr = (len(self.iblrs)+len(self.ibars)+
                 len(self.itors)+len(self.iobrs)+len(self.iodrs))
        with open(filename, 'r') as f:
            with open(self.ioPath+'intnrdef.txt', 'w') as fout:
                line = aux.skipComment(f, '#')
                nintnr = int(line)
                line = aux.skipComment(f, '#')
                line = np.array(line.split(), dtype=int)
                [nblnr, nbanr, ntonr, nobnr, nodnr] = line[:]
                print >>fout, nintnr
                ncount = 0
                tmat = np.zeros((nintnr, nintr), dtype=float)
                while ncount < nintnr:
                    line = aux.skipComment(f, '#')
                    if '=' in line:
                        line = line.split()
                        inrstart = int(line[0])
                        inrend = int(line[2])
                        irstart = int(line[4])
                        irend = int(line[6])
                        if irend-irstart != inrend-inrstart:
                            raise RuntimeError(
                                  '*** Problem in definition of nonredundant internals. Stop')
                        #-- array index convention needs to be handled VERY CAREFULLY
                        for i in xrange(inrend-inrstart+1):
                            print >>fout, inrstart+i, 1
                            print >>fout, ' ', irstart+i
                            print >>fout, '  1.0'
                            tmat[inrstart+i-1][irstart+i-1] = 1.0
                        ncount = inrend
                    else:
                        print >>fout, line,
                        inr = int(line.split()[0])
                        line = aux.skipComment(f, '#')
                        line = np.array(line.split(), dtype=int)
                        line2 = aux.skipComment(f, '#')
                        line2 = np.array(line2.split(), dtype=float)
                        line2 = line2/np.linalg.norm(line2)
                        #-- array index convention needs to be handled VERY CAREFULLY
                        for i in xrange(len(line)):
                            tmat[inr-1][line[i]-1] = line2[i]
                        print >>fout, ' ',
                        np.savetxt(fout, line, fmt='%d', newline=' ')
                        print >>fout
        
                        print >>fout, ' ',
                        np.savetxt(fout, line2, 
                                   fmt='%15.10f', newline=' '); 
                        print >>fout
        
                        ncount += 1
        self.tmat = tmat
        self.nblnr = nblnr
        self.nbanr = nbanr
        self.ntonr = ntonr
        self.nobnr = nobnr
        self.nodnr = nodnr
        self.nintnr = nintnr


    def calcABmat(self):
        """Calculate A and B matrices for the current geometry.

        Must have read in the Cartesian, redundant int., 
          and nonredundant int.
        """
        if self.carts is None:
            raise RuntimeError(
                  '*** Cartesian coordinates not set. Run readCarts first. ***')
        if self.iblrs is None:
            raise RuntimeError(
                  '*** Redundant internal coordinates not set. Run readRedundIntCoords first. ***')
        if self.tmat is None:
            raise RuntimeError(
                  '*** Nonredundant internal coordinates not set. Run readNonRedundIntCoords first. ***')

        #== Prepare data
        #-- indexes should be converted from 0-based to 1-based
        #--   which will be printed to files used by abmat.exe
        iblrs = self.iblrs+1; nblr = len(iblrs)
        ibars = self.ibars+1; nbar = len(ibars)
        itors = self.itors+1; ntor = len(itors)
        iobrs = self.iobrs+1; nobr = len(iobrs)
        iodrs = self.iodrs+1; nodr = len(iodrs)
        nintr = nblr+nbar+ntor+nobr+nodr
        ioPath = self.ioPath
        atoms = self.atoms
        masses = self.masses
        carts = self.carts
        natom = len(atoms)

        #== Prepare files for abmat.exe
        with open(ioPath+'ibl.txt', 'w') as f:
            print >>f, nblr
            np.savetxt(f, iblrs, fmt='%d')
        with open(ioPath+'iba.txt', 'w') as f:
            print >>f, nbar
            np.savetxt(f, ibars, fmt='%d')
        with open(ioPath+'ito.txt', 'w') as f:
            print >>f, ntor
            np.savetxt(f, itors, fmt='%d')
        with open(ioPath+'iob.txt', 'w') as f:
            print >>f, nobr
            np.savetxt(f, iobrs, fmt='%d')
        with open(ioPath+'iod.txt', 'w') as f:
            print >>f, nodr
            np.savetxt(f, iodrs, fmt='%d')

        with open(ioPath+'cart.txt', 'w') as f:
            print >>f, natom
            for i in xrange(natom):
                print >>f, '%2s %8.3f %15.10f %15.10f %15.10f'%(
                           atoms[i], masses[i], carts[i*3], 
                           carts[i*3+1], carts[i*3+2])
            
        #-- Calculate initial A and B matrices
        os.system(self.codePath+'abmat.exe '+ioPath)

        #-- Read in A matrix
        with open(ioPath+'amat.txt', 'r') as famat:
            line = famat.readline().split()
            ncart, nintnr = int(line[0]), int(line[1])
            if ncart != natom*3:
                raise RuntimeError(
                      '*** Number of Cartesian coordinates != Number of atoms * 3. Something is wrong here. ***')
            self.amat = np.loadtxt(famat)
        #-- Read in B matrix
        with open(ioPath+'bmat.txt', 'r') as fbmat:
            line = fbmat.readline().split()
            self.bmat = np.loadtxt(fbmat)


    def displaceCarts(self, filenameDisplNR=None, updateSelf=False,
                      printOutput=False, fileHandleOut=sys.stdout):
        """Displace Cartesian coordinates.

        filenameDisplNR: Name of input file of displacement of nonredund. int. coord.
                         Default to self.ioPath+'nrdispl0.txt'
        updateSelf: If True, update self.carts to the displaced Carts;
                    If False, return an instance of GeomClass with updated Carts.
                    Default to False
        printOutput: Whether or not print the displaced Carts and convergence info.
                     Default to False
        fileHandleOut: File handle to which the output is printed.
                       Default to standard output
        """
        if filenameDisplNR is None:
            filenameDisplNR = self.filenameDisplNR
        if filenameDisplNR is None:
           filenameDisplNR = self.ioPath+'nrdispl0.txt'
           self.filenameDisplNR = filenameDisplNR

        #-- calculate A and B matrices if not done
        if self.amat is None:
            self.calcABmat()
        #== Prepare data
        iblrs = self.iblrs; nblr = len(iblrs)
        ibars = self.ibars; nbar = len(ibars)
        itors = self.itors; ntor = len(itors)
        iobrs = self.iobrs; nobr = len(iobrs)
        iodrs = self.iodrs; nodr = len(iodrs)
        nintr = nblr+nbar+ntor+nobr+nodr
        ioPath = self.ioPath
        atoms = self.atoms
        masses = self.masses
        carts = self.carts
        natom = len(atoms)
        nblnr = self.nblnr
        nbanr = self.nbanr
        ntonr = self.ntonr
        nobnr = self.nobnr
        nodnr = self.nodnr
        nintnr = self.nintnr
        tmat = self.tmat
        amat = self.amat
        
        #-- Read in nonredundant internal displacements
        with open(filenameDisplNR, 'r') as fdispl:
            dqsnr = np.zeros(nintnr, dtype=float)
            while True:
                line = aux.skipComment(fdispl, '#')
                if len(line)==0:
                    break
                line = line.split()
                dqsnr[int(line[0])-1] = float(line[1])
                #-- convert degree to rad for angles
                if nblnr<int(line[0])<=nblnr+nbanr+ntonr+nobnr:
                    dqsnr[int(line[0])-1] = float(line[1])/180*np.pi

        #-- Determine internal coordinates of the initial geometry
        qsr = aux.calcInt(nblr, nbar, ntor, nobr, nodr, 
                          iblrs, ibars, itors, iobrs, iodrs, 
                          np.reshape(carts, (natom, 3)))
        qsnr = np.dot(tmat, qsr)

        #== Iteratively determine displaced Cartesian coordinates 
        #==   from initial Cart. coord. and internal displacements
        #== See Pulay et al. JACS 1979, 101, 2550. 
        
        #-- x1 = x0+A*dq
        if dqsnr.size==1:
            #-- np.dot would not work if there's only one internal 
            #--   coord. and dqsnr has only one element
            cartsnew = carts+amat*dqsnr
        else:
            cartsnew = carts+np.dot(amat, dqsnr)
        
        cycle = 0
        while True:
            cycle += 1
            if cycle > 10:
                print >>sys.stderr, '*** Exceeds 10 iterations. Stop.                         ***'
                print >>sys.stderr, '*** Unconvergence may be due to too large displacements. ***'
                raise RuntimeError
            #-- Calculate new internals from displaced Cartesians
            qsrnew = aux.calcInt(nblr, nbar, ntor, nobr, nodr, 
                                 iblrs, ibars, itors, iobrs, iodrs, 
                                 np.reshape(cartsnew, (natom, 3)))
            dqsrnew = qsrnew-qsr
            #-- Deal with "stepping over -pi|pi" problem for torsions
            for i in xrange(nblr+nbar, nblr+nbar+ntor):
                if abs(dqsrnew[i]) > 2*np.pi-abs(dqsrnew[i]):
                    dqsrnew[i] = -np.sign(dqsrnew[i])*(
                                  2*np.pi-abs(dqsrnew[i]))
            dqsnrnew = np.dot(tmat, dqsrnew)
            
            #-- Update A matrix
            with open(ioPath+'cart.txt', 'w') as f:
                print >>f, natom
                for i in xrange(natom):
                    print >>f, '%2s %8.3f %15.10f %15.10f %15.10f'%(
                               atoms[i], masses[i], cartsnew[i*3], 
                               cartsnew[i*3+1], cartsnew[i*3+2])
            os.system(self.codePath+'abmat.exe '+ioPath)
            with open(ioPath+'amat.txt', 'r') as famat:
                famat.readline()
                amat = np.loadtxt(famat)
            #-- x2 = x1+A*(dq-dq1) 
            if dqsnr.size==1:
                #-- np.dot would not work if there's only one internal 
                #--   coord. and dqsnr has only one element
                cartsnewer = cartsnew+amat*(dqsnr-dqsnrnew)
            else:
                cartsnewer = cartsnew+np.dot(amat, dqsnr-dqsnrnew)
            #-- Check convergence in terms of max change in Cartesians
            if np.amax(np.absolute(cartsnewer-cartsnew))<1e-6:
                if printOutput:
                    print >>fileHandleOut, 'Algorithm converged in', cycle, 'iterations. Max Cartesian change between the last two iterations below 1e-6.'
                break
            cartsnew = cartsnewer

        #== Finalize
        #-- new object with updated Carts
        newObject = copy.deepcopy(self)
        newObject.carts = cartsnew
        #-- calculate T matrix
        newObject.readNonRedundIntCoords() 
        #-- calculate A and B matrices
        newObject.calcABmat()
        if updateSelf:
            #-- copy so that self and returned are different objects
            self.carts = cartsnew
            self.readNonRedundIntCoords() 
            self.calcABmat()

        #-- Print final info and Cartesian coordinates
        if printOutput:
            print >>fileHandleOut, 'Intended change in NR internals:'
            for i in xrange(nintnr):
                if abs(dqsnr[i])>1e-10:
                    if i<nblnr or i>=nblnr+nbanr+ntonr+nobnr:
                        unit = 'Angstroms'
                        num = dqsnr[i]
                    else:
                        unit = 'degrees'
                        num = dqsnr[i]/np.pi*180
                    print >>fileHandleOut, '%12.4f %s in internal # %d'%(
                                           num, unit, i+1)
            rmserrlen = np.sqrt(np.linalg.norm(
                        (dqsnrnew-dqsnr)[:nblnr]
                         +(dqsnrnew-dqsnr)[nblnr+nbanr+ntonr+nobnr:]
                        )**2/(nblnr+nodnr))
            rmserrang = np.sqrt(np.linalg.norm(
                        (dqsnrnew-dqsnr)[nblnr:nblnr+nbanr+ntonr+nobnr]
                        )**2/(nintnr-nblnr-nodnr))/np.pi*180
            maxerr = np.amax(np.absolute(dqsnrnew-dqsnr))
            imaxerr = np.argmax(np.absolute(dqsnrnew-dqsnr))+1
            if nblnr < imaxerr <= nblnr+nbanr+ntonr+nobnr:
                maxerr = maxerr/np.pi*180
                uniterr = 'degrees'
            else:
                uniterr = 'Angstroms'
            print >>fileHandleOut, 'RMS error in NR lengths       = %.1e Angstroms'%(rmserrlen)
            print >>fileHandleOut, 'RMS error in NR angles        = %.1e degrees'%(rmserrang)
            print >>fileHandleOut, 'Max abs error in NR internals = %.1e %s in NR internal # %d'%(
                  maxerr, uniterr, imaxerr)
            print >>fileHandleOut
            print >>fileHandleOut, 'Final Cartesian coordinates:'
            newObject.printCarts(fileHandleOut)

        return newObject


    def calcIntCoord(self, ctype, indexstr):
        """Calculate an internal coordinate.
    
        ctype: type of internal coord. to calculate 
               1 or 'bl': bond length
               2 or 'ba': bond angle
               3 or 'to': torsion
               4 or 'ob': out-of-plane bend
               5 or 'od': out-of-plane distance
        indexstr: list of atom indexes (1-based) 
        Return: internal coord., unit
                where unit is 'Angstroms' for lengths and 'degrees' for angles
        """
        x = np.reshape(self.carts, (len(self.atoms), 3))
        # try:
        indexstr = np.array(indexstr, dtype=int)-1
        #-- primary 'torsional' coordinate for thioanisole
        if ctype in [-1,'phi']:
            i,j,k,l,m = 1,2,3,11,12
            e12 = (x[j,:]-x[i,:])/aux.calcBL([i, j], x)
            e23 = (x[k,:]-x[j,:])/aux.calcBL([j, k], x)
            e13 = (x[k,:]-x[i,:])/aux.calcBL([i, k], x)
            e45 = (x[m,:]-x[l,:])/aux.calcBL([l, m], x)
            e24 = (x[l,:]-x[j,:])/aux.calcBL([j, l], x)
            etmp = np.cross(e12, e23)
            eZ = etmp/np.linalg.norm(etmp)
            etmp = np.cross(e13, eZ) 
            eX = etmp/np.linalg.norm(etmp)
            etmp = e45-eX*np.dot(e45, eX)
            ep = etmp/np.linalg.norm(etmp)
            cosSP1 = np.dot(ep, e13)
        
            #-- determine sign
            etmp = np.cross(e13, ep)
            if np.dot(etmp, eX)>0:
                sgn = 1
            else:
                sgn = -1
            #-- avoid numerical problem
            if cosSP1>1:
                evalSP1 = 0e0
            elif cosSP1<-1:
                evalSP1 = sgn*np.pi
            else:
                evalSP1 = sgn*np.arccos(cosSP1)
        
            return evalSP1/np.pi*180, 'degrees'
            
        #-- bond length
        elif ctype in [1,'bl']:
            return aux.calcBL(indexstr, x), 'Angstroms'
        #-- bond angle
        elif ctype in [2,'ba']:
            return aux.calcBA(indexstr, x)/np.pi*180, 'degrees'
        #-- torsion
        elif ctype in [3,'to']:
            return aux.calcTO(indexstr, x)/np.pi*180, 'degrees'
        #-- out-of-plane bend 1-2-3-4: angle b/ 1-4 & 2-4-3
        elif ctype in [4,'ob']:
            return aux.calcOB(indexstr, x)/np.pi*180, 'degrees'
        #-- out-of-plane distance 1-2-3-4: distance b/ 4 & 2-4-3
        elif ctype in [5,'od']:
            return aux.calcOD(indexstr, x), 'Angstroms'
        # except:
        #     print 'Something wrong in calcIntCoord'
        #     return None


    def printCarts(self, fileHandleOut=sys.stdout):
        """Print Cartesian coordinates.

        fileHandleOut: File handle to which the output is printed.
                       Default to standard output
        """
        atoms = self.atoms
        carts = self.carts
        for i in xrange(len(atoms)):
            print >>fileHandleOut, '%2s %12.7f %12.7f %12.7f'%(
                  atoms[i], carts[i*3], carts[i*3+1], 
                  carts[i*3+2])

    
    def printAllRedundIntCoords(self, fh=sys.stdout):
        """Print all defined redundant internal coordinates. 

        fh: file handle for output
        """
        names = ['Bond lengths', 'Bond angles', 'Torsions', 
                 'OOP bends','OOP distances']
        coords = [self.iblrs, self.ibars, self.itors, 
                  self.iobrs, self.iodrs]
        for i in xrange(len(names)):
            print >>fh, names[i]+':'
            for idxes in coords[i]:
                print >>fh, '  ', idxes+1, ': ',
                print >>fh, '% .6f %s'%self.calcIntCoord(i+1, idxes+1)
