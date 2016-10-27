import pdb

import numpy as np
import displace_int2cart_aux as aux
import os
import sys

#================================
#=== Main program starts here ===
#================================
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#--- Naming conventions:
#--- ibl, iba, ito, iob = bond length, bond angle, torsion, oop bend
#--- cart, int = Cartesian, internal
#--- -nr, -r = nonredundant, redundant
#--- -s = array
#--- n- = number of
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#--- Set up paths
if len(sys.argv) > 1:
    ioPath = sys.argv[1]
    if ioPath[-1] != '/': ioPath += '/'
else:
    ioPath = os.getcwd()+'/'

codePath = os.path.dirname(__file__)
if not codePath:
    codePath = '.'
codePath += '/'
#--- Read in initial Cartesian coordinates
with open(ioPath+'cart0.txt', 'r') as fcart:
    line = aux.skipcomment(fcart, '#').split()
    natom, irun = int(line[0]), int(line[1])
    atoms = []; masses = []; cartslist = []
    for i in range(natom):
        line = aux.skipcomment(fcart, '#').split()
        atoms.append(line[0])
        masses.append(float(line[1]))
        for j in range(3):
            #--- Setting the zero coord. to 1e-10 to prevent some numerical problems
            if abs(float(line[j+2]))<1e-10:
                cartslist.append(1e-10)
            else:
                cartslist.append(float(line[j+2]))
    carts = np.array(cartslist)
    del cartslist
#--- Read in definition of (redundant) internal coordinates 
#--- and write files for calc. A matrix
iblrs, ibars, itors, iobrs, iodrs = [], [], [], [], []
with open(ioPath+'rdef0.txt', 'r') as f:
    #--- bond lengths
    line = aux.skipcomment(f, '#')
    niblr = int(line)
    for i in range(niblr):
        line = aux.skipcomment(f, '#')
        iblrs.append(line.split())
    #--- bends
    line = aux.skipcomment(f, '#')
    nibar = int(line)
    for i in range(nibar):
        line = aux.skipcomment(f, '#')
        ibars.append(line.split())
    #--- torsions
    line = aux.skipcomment(f, '#')
    nitor = int(line)
    for i in range(nitor):
        line = aux.skipcomment(f, '#')
        itors.append(line.split())
    #--- oop bends
    line = aux.skipcomment(f, '#')
    niobr = int(line)
    for i in range(niobr):
        line = aux.skipcomment(f, '#')
        iobrs.append(line.split())
    #--- oop distances
    line = aux.skipcomment(f, '#')
    niodr = int(line)
    for i in range(niodr):
        line = aux.skipcomment(f, '#')
        iodrs.append(line.split())
    iblrs = np.array(iblrs, dtype=int)
    ibars = np.array(ibars, dtype=int)
    itors = np.array(itors, dtype=int)
    iobrs = np.array(iobrs, dtype=int)
    iodrs = np.array(iodrs, dtype=int)

with open(ioPath+'ibl.txt', 'w') as f:
    print >>f, niblr
    np.savetxt(f, iblrs, fmt='%d')
with open(ioPath+'iba.txt', 'w') as f:
    print >>f, nibar
    np.savetxt(f, ibars, fmt='%d')
with open(ioPath+'ito.txt', 'w') as f:
    print >>f, nitor
    np.savetxt(f, itors, fmt='%d')
with open(ioPath+'iob.txt', 'w') as f:
    print >>f, niobr
    np.savetxt(f, iobrs, fmt='%d')
with open(ioPath+'iod.txt', 'w') as f:
    print >>f, niodr
    np.savetxt(f, iodrs, fmt='%d')
nintr = niblr+nibar+nitor+niobr+niodr
#--- account for the array index convention
iblrs, ibars, itors, iobrs, iodrs = iblrs-1, ibars-1, itors-1, iobrs-1, iodrs-1

#--- Read in definition of nonredundant coordinates
#--- and prepare input for abmat.exe
with open(ioPath+'nrdef0.txt', 'r') as f:
    with open(ioPath+'intnrdef.txt', 'w') as fout:
        line = aux.skipcomment(f, '#')
        nintnr = int(line)
        line = aux.skipcomment(f, '#')
        line = np.array(line.split(), dtype=int)
        [niblnr, nibanr, nitonr, niobnr, niodnr] = line[:]
        print >>fout, nintnr
        ncount = 0
        tmat = np.zeros((nintnr, nintr), dtype=float)
        while ncount<nintnr:
            line = aux.skipcomment(f, '#')
            if '=' in line:
                line = line.split()
                inrstart = int(line[0])
                inrend = int(line[2])
                irstart = int(line[4])
                irend = int(line[6])
                if irend-irstart!=inrend-inrstart:
                    print '*** Problem in definition of nonredundant internals. Stop'
                    quit()
                #--- array index convention needs to be handled VERY CAREFULLY
                for i in range(inrend-inrstart+1):
                    print >>fout, inrstart+i, 1
                    print >>fout, ' ', irstart+i
                    print >>fout, '  1.0'
                    tmat[inrstart+i-1][irstart+i-1] = 1.0
                ncount = inrend
            else:
                print >>fout, line,
                inr = int(line.split()[0])
                line = aux.skipcomment(f, '#')
                line = np.array(line.split(), dtype=int)
                line2 = aux.skipcomment(f, '#')
                line2 = np.array(line2.split(), dtype=float)
                line2 = line2/np.linalg.norm(line2)
                #--- array index convention needs to be handled VERY CAREFULLY
                for i in range(len(line)):
                    tmat[inr-1][line[i]-1] = line2[i]
                print >>fout, ' ',
                np.savetxt(fout, line, fmt='%d', newline=' ')
                print >>fout

                print >>fout, ' ',
                np.savetxt(fout, line2, fmt='%15.10f', newline=' '); 
                print >>fout

                ncount += 1
# for i in range(nintnr):
#     for j in range(nint):
#         print tmat[i][j],
#     print

#--- Calculate initial A matrix
with open(ioPath+'cart0.txt', 'r') as fin:
    with open(ioPath+'cart.txt', 'w') as fout:
        while True:
            line = aux.skipcomment(fin, '#')
            if(len(line)==0):
                break
            print >>fout, line

#os.system('cat cart0.txt >cart.txt')
os.system(codePath+'abmat.exe '+ioPath)

#--- Stop here if it is a test run
if irun==0:
    quit()    

#--- Read in A matrix
with open(ioPath+'amat.txt', 'r') as famat:
    line = famat.readline().split()
    ncart, nintnr = int(line[0]), int(line[1])
    if ncart!=natom*3:
        print 'Number of Cartesian coordinates != Number of atoms * 3. Something is wrong here.'
        quit()
    amat = np.loadtxt(famat)
#--- Read in B matrix
# with open(path+'bmat.txt', 'r') as fbmat:
#     line = fbmat.readline().split()
#     bmat = np.loadtxt(fbmat)


#--- Read in nonredundant internal displacements
with open(ioPath+'nrdispl0.txt', 'r') as fdispl:
    dqsnr = np.zeros(nintnr, dtype=float)
    while True:
        line = aux.skipcomment(fdispl, '#')
        if len(line)==0:
            break
        line = line.split()
        dqsnr[int(line[0])-1] = float(line[1])
        #--- convert degree to rad for angles
        if niblnr<int(line[0])<=niblnr+nibanr+nitonr+niobnr:
            dqsnr[int(line[0])-1] = float(line[1])/180*np.pi

#--- Determine internal coordinates of the initial geometry
qsr = aux.calcInt(niblr, nibar, nitor, niobr, niodr, iblrs, ibars, itors, iobrs, iodrs, np.reshape(carts, (natom, 3)))
qsnr = np.dot(tmat, qsr)

#=== Iteratively determine displaced Cartesian coordinates from initial Cart. coord. and internal displacements
#=== See Pulay et al. JACS 1979, 101, 2550. 

#--- x1 = x0+A*dq
if dqsnr.size==1:
    #--- np.dot would not work if there's only one internal coord. and dqsnr has only one element
    cartsnew = carts+amat*dqsnr
else:
    cartsnew = carts+np.dot(amat, dqsnr)

cycle = 0
#print 'iter',
while True:
    cycle += 1
#    print cycle,
    if cycle > 10:
        print '*** Exceeds 10 iterations. Stop.                         ***'
        print '*** Unconvergence may be due to too large displacements. ***'
        break
    #--- Calculate new internals from displaced Cartesians
    qsrnew = aux.calcInt(niblr, nibar, nitor, niobr, niodr, 
                         iblrs, ibars, itors, iobrs, iodrs, 
                         np.reshape(cartsnew, (natom, 3)))
    dqsrnew = qsrnew-qsr
    #--- Deal with "stepping over -pi|pi" problem for torsions
    for i in range(niblr+nibar, niblr+nibar+nitor):
        if abs(dqsrnew[i]) > 2*np.pi-abs(dqsrnew[i]):
            dqsrnew[i] = -np.sign(dqsrnew[i])*(2*np.pi-abs(dqsrnew[i]))
    dqsnrnew = np.dot(tmat, dqsrnew)
    
    #--- Update A matrix
    with open(ioPath+'cart.txt', 'w') as f:
        print >>f, natom
        for i in range(natom):
            print >>f, '%2s %8.3f %15.10f %15.10f %15.10f'%(
                       atoms[i], masses[i], cartsnew[i*3], 
                       cartsnew[i*3+1], cartsnew[i*3+2])
    os.system(codePath+'abmat.exe '+ioPath)
    with open(ioPath+'amat.txt', 'r') as famat:
        famat.readline()
        amat = np.loadtxt(famat)
    #--- x2 = x1+A*(dq-dq1) 
    if dqsnr.size==1:
        #--- np.dot would not work if there's only one internal coord. and dqsnr has only one element
        cartsnewer = cartsnew+amat*(dqsnr-dqsnrnew)
    else:
        cartsnewer = cartsnew+np.dot(amat, dqsnr-dqsnrnew)
    #--- Check convergence in terms of max change in Cartesians
    if np.amax(np.absolute(cartsnewer-cartsnew))<1e-6:
        print 'Converged in', cycle, 'iterations. Max Cartesian change below 1e-6.'
        break
    cartsnew = cartsnewer
    

#--- Print final info and Cartesian coordinates
print 'Intended change in NR internals:'
for i in range(nintnr):
    if abs(dqsnr[i])>1e-10:
        if i<niblnr or i>=niblnr+nibanr+nitonr+niobnr:
            unit = 'Angstroms'
            num = dqsnr[i]
        else:
            unit = 'degrees'
            num = dqsnr[i]/np.pi*180
        print '%12.4f %s in internal # %d'%(num, unit, i+1)
rmserrlen = np.sqrt(np.linalg.norm(
                    (dqsnrnew-dqsnr)[:niblnr]
                   +(dqsnrnew-dqsnr)[niblnr+nibanr+nitonr+niobnr:]
                                  )**2/(niblnr+niodnr))
rmserrang = np.sqrt(np.linalg.norm(
                    (dqsnrnew-dqsnr)[niblnr:niblnr+nibanr+nitonr+niobnr]
                                  )**2/(nintnr-niblnr-niodnr))/np.pi*180
maxerr = np.amax(np.absolute(dqsnrnew-dqsnr))
imaxerr = np.argmax(np.absolute(dqsnrnew-dqsnr))+1
if niblnr < imaxerr <= niblnr+nibanr+nitonr+niobnr:
    maxerr = maxerr/np.pi*180
    uniterr = 'degrees'
else:
    uniterr = 'Angstroms'
print 'RMS error in NR lengths       = %.1e Angstroms'%(rmserrlen)
print 'RMS error in NR angles        = %.1e degrees'%(rmserrang)
print 'Max abs error in NR internals = %.1e %s in NR internal # %d'%(
      maxerr, uniterr, imaxerr)
print 'Final Cartesian coordinates:'

for i in range(natom):
    print '%2s %12.7f %12.7f %12.7f'%(
          atoms[i], cartsnew[i*3], cartsnew[i*3+1], cartsnew[i*3+2])

