import numpy as np

#== Aux. functions
def calcBL(iatomlist, carts):
    '''Calculate bond length.

    iatomlist: (integer list) index of atoms
    carts: (natom-by-3 numpy array) Cartesian coordinates
    '''
    return np.linalg.norm(carts[iatomlist[0]]-carts[iatomlist[1]])


def calcBA(iatomlist, carts):
    '''Calculate bond angle 1-2-3 in radian.
    
    iatomlist: (integer list) index of atoms
    carts: (natom-by-3 numpy array) Cartesian coordinates
    '''
    unit21 = carts[iatomlist[0]]-carts[iatomlist[1]]
    unit21 = unit21/np.linalg.norm(unit21)
    unit23 = carts[iatomlist[2]]-carts[iatomlist[1]]
    unit23 = unit23/np.linalg.norm(unit23)
    cosBA = np.dot(unit21, unit23)
    #-- avoid numerical problem
    if cosBA>1:
        return 0.0
    elif cosBA<-1:
        return np.pi
    else:
        return np.arccos(cosBA)
        

def calcTO(iatomlist, carts):
    '''Calculate torsion 1-2-3-4 in radian.

    iatomlist: (integer list) index of atoms
    carts: (natom-by-3 numpy array) Cartesian coordinates
    '''
    unit12 = carts[iatomlist[1]]-carts[iatomlist[0]]
    unit12 = unit12/np.linalg.norm(unit12)
    unit23 = carts[iatomlist[2]]-carts[iatomlist[1]]
    unit23 = unit23/np.linalg.norm(unit23)
    unit34 = carts[iatomlist[3]]-carts[iatomlist[2]]
    unit34 = unit34/np.linalg.norm(unit34)
    sin123 = np.sin(calcBA(iatomlist[0:3], carts))
    sin234 = np.sin(calcBA(iatomlist[1:4], carts))
    cosTO =  np.dot( np.cross(unit12, unit23), (np.cross(unit23, unit34)) )\
        /sin123/sin234 

    # if 15 in iatomlist:
    #     pdb.set_trace()

    #-- determine sign of torsion
    # if np.dot(np.cross(-unit12, unit34), unit23)>0:
    if np.dot(np.cross( np.cross(unit12, unit23), np.cross(unit23, unit34) ), unit23)>0:
        sign = 1
    else:
        sign = -1
    #-- avoid numerical problem
    if cosTO>1:
        return 0.0
    elif cosTO<-1:
        return sign*np.pi
    else:
        return sign*np.arccos(cosTO)


def calcOB(iatomlist, carts):
    '''Calculate out-of-plane bend 1-2-3-4 in radian.
    (definition: angle between bond 1-4 and plane 2-4-3)

    iatomlist: (integer list) index of atoms
    carts: (natom-by-3 numpy array) Cartesian coordinates
    '''
    unit41 = carts[iatomlist[0]]-carts[iatomlist[3]]
    unit41 = unit41/np.linalg.norm(unit41)
    unit42 = carts[iatomlist[1]]-carts[iatomlist[3]]
    unit42 = unit42/np.linalg.norm(unit42)
    unit43 = carts[iatomlist[2]]-carts[iatomlist[3]]
    unit43 = unit43/np.linalg.norm(unit43)
    atom243 = [iatomlist[1], iatomlist[3], iatomlist[2]]
    sin243 = np.sin(calcBA(atom243, carts))
    sinOB = np.dot( np.cross(unit42, unit43), unit41 )/sin243
    return np.arcsin(sinOB)

def calcOD(iatomlist, carts):
    '''Calculate out-of-plane distance 1-2-3-4 in radian.
    (definition: distance of center atom 4 from plane 1-2-3)

    iatomlist: (integer list) index of atoms
    carts: (natom-by-3 numpy array) Cartesian coordinates
    '''
    unit12 = carts[iatomlist[1]]-carts[iatomlist[0]]
    unit12 = unit12/np.linalg.norm(unit12)
    unit13 = carts[iatomlist[2]]-carts[iatomlist[0]]
    unit13 = unit13/np.linalg.norm(unit13)
    v14    = carts[iatomlist[3]]-carts[iatomlist[0]]
    atom213 = [iatomlist[1], iatomlist[0], iatomlist[2]]
    sin213 = np.sin(calcBA(atom213, carts))
    return abs( np.dot(v14, np.cross(unit12, unit13))/sin213 )

def calcInt(nibl, niba, nito, niob, niod, ibls, ibas, itos, iobs, iods, carts):
    '''Calculate all primitive internals from Cartesians.

    nibl, niba, nito, niob: (integer) number of bond lengths, bond angles, torsions, oop bends
    ibls, ibas, itos, iobs: (2D numpy arrays) each row is a list of indexes defining the internal
    carts: (natom-by-3 numpy array) Cartesian coordinates
    return: (1D numpy array) calculated internals
    '''
    qs = []
    for i in range(nibl):
        qs.append(calcBL(ibls[i], carts))
    for i in range(niba):
        qs.append(calcBA(ibas[i], carts))
    for i in range(nito):
        qs.append(calcTO(itos[i], carts))
    for i in range(niob):
        qs.append(calcOB(iobs[i], carts))
    for i in range(niod):
        qs.append(calcOD(iods[i], carts))
    return np.array(qs)

def skipComment(f, char='#'):
    """Skip comment lines and empty lines when reading a file
    and return the first non-comment, non-empty line.

    f: file handle
    char: character denoting a comment line, e.g. '#'
    """
    while True:
        line = f.readline()
        linesplit = line.split()
        #-- EOF
        if not line:
            break
        #-- comment or empty line
        elif line.startswith(char) or not linesplit:
            continue
        else:
            break
    return line
    
