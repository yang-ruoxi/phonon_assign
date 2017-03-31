#!/usr/bin/python

import numpy as np
import math as m

# Extract index of last line in list at which string occurs. Within optional min-max range.
def strindex(list,string,nmin=0,nmax=0):
    if nmax == 0:
        nmax = len(list)
    i = 0
    index = 0 
    for item in list:
        if string in item:
            if i >= nmin and i <= nmax:
                index = i
        i = i+1
    return index

# Extract list of indicies of lines in list at which string occurs. Within optional min-max range.
def strindicies(list,string,nmin=0,nmax=0):
    if nmax == 0:
        nmax = len(list)
    indicies = []
    i = 0
    for item in list:
        if string in item:
            if i >= nmin and i <= nmax:
                indicies.append(i)
        i = i+1
    return indicies

# Extract number of ions from .yaml file3
def get_no_ions(yamllines):
    nions = strindex(yamllines,'natom')
    N = int(yamllines[nions].split()[1])
    return N

# Extract lattice vectors from POSCAR file
def get_lat_vecs(flines):
    lat = float(flines[1].split()[0])
    avec = [float(flines[2].split()[0])*lat,float(flines[2].split()[1])*lat,float(flines[2].split()[2])*lat]
    bvec = [float(flines[3].split()[0])*lat,float(flines[3].split()[1])*lat,float(flines[3].split()[2])*lat]
    cvec = [float(flines[4].split()[0])*lat,float(flines[4].split()[1])*lat,float(flines[4].split()[2])*lat]
    latvec = [avec, bvec, cvec]
    return latvec

# Forms a list of all q-points in the yaml file
def get_qpts(yamllines):
    nqpts = strindicies(yamllines,'q-position')
    qpts = []
    for nqpt in nqpts:
        qpt = [float(yamllines[nqpt].split()[3].replace(',','')),float(yamllines[nqpt].split()[4].replace(',','')),float(yamllines[nqpt].split()[5].replace(',',''))]
        qpts.append(qpt)
    return qpts

# Extract the index of a given qpt in the .yaml file and the index of when that qpt ends
def qpt_index_range(yamllines,qpt):
    nqpts = strindicies(yamllines,'q-position')
    nqpoint = 0
    for nqpt in nqpts:
        if float(yamllines[nqpt].split()[3].replace(',',''))==qpt[0] and float(yamllines[nqpt].split()[4].replace(',',''))==qpt[1] and float(yamllines[nqpt].split()[5].replace(',',''))==qpt[2]:
            nqpoint = nqpt
    if nqpoint == 0:
        print "Requested q-point: " + ' '.join([str(q) for q in qpt]) + " does not appear in output file.\nPlease request a different set of q-points."
        exit()
    if nqpts.index(nqpoint) == len(nqpts)-1:
        nqmax = len(yamllines)
    else:
        nqmax = nqpts[nqpts.index(nqpoint)+1]
    return (nqpoint, nqmax)

# Extract either absolute or fractional ion positions from POSCAR file
def get_ion_posns(flines,Nions,frac=0):
    latvec = get_lat_vecs(flines)
    npos = strindex(flines,'Direct') + 1
    fracposns = []
    posns = []
    for i in range(0,Nions):
        line = flines[npos+i]
        pos = [float(line.split()[0]),float(line.split()[1]),float(line.split()[2])]
        fracposns.append(pos)
    for fpos in fracposns:
        if frac == 0:                               # If frac option not selected, output absolute positions
            posns.append(np.dot(latvec,fpos))
        elif frac == 1:                             # If frac option is selected, output fractional positions
            posns.append(fpos)
    return np.array(posns)

# Extract elements from POSCAR file
def get_elements(flines,Nions):
    elements = []
    symbol = flines[5].split()
    number = flines[6].split()
    for i in range(0,len(symbol)):
        if number[i] == '1':
            elements.append(symbol[i])              # append element once if there's only one in POSCAR
        else:
            elements.extend(symbol[i]*int(number[i]))  # otherwise append to the element list as many as it occurs in POSCAR
    return elements

# Extract masses from .phonon file  ##############################THERE IS NO MASS###########################
def get_masses(flines,Nions):
    npos = strindex(flines,'Fractional Co-ordinates') + 1
    masses = []
    for i in range(0,Nions):
        line = flines[npos+i]
        masses.append(float(line.split()[5]))
    return masses

# Extracts a list of the frequency of each branch at a given qpt for a list of branches from mesh.yaml
def get_freqencies(yamllines,branches,qpt):
    (nqpoint, nqmax) = qpt_index_range(yamllines,qpt)
    freqs = []
    nfrequencies = strindicies(yamllines, 'frequency:')   # get the line numbers for "frequency"
    for nfrequency in nfrequencies:                                 
        if nfrequency < (nqpoint + ((4*Nions+3)*(3*Nions)+3)):          # when the frequency index is within the limit for this qpoint upto next qpoint
            freqs.append(float(yamllines[nfrequency].split()[1]))       
    #for b in branches:
    #    freqs.append(float(yamllines[nqpoint + b].split()[1]))
    return freqs

# Extracts real and imaginary components of phonon eigenvectors as separate lists from yaml file
def get_eigvects(yamllines,Nions,branch,qpt):
    #masses = get_masses(flines,Nions)
    (nqpoint, nqmax) = qpt_index_range(yamllines,qpt)
    nvectsstart = strindex(yamllines,'eigenvector',nmin=nqpoint,nmax=nqmax)         # eigenvector starts from the last branch
    nbranch = nvectsstart - int(15-branch)*(3+int(Nions)*4)                         # so we count back from the last branch 
    i = 0
    Reigvects = []
    Ieigvects = []
    for i in range(0,Nions):        
        x_vec = yamllines[nbranch+i*4+2]
        y_vec = yamllines[nbranch+i*4+3]
        z_vec = yamllines[nbranch+i*4+4]                                        # get the x_vec y_vec z_vec on each atom
        #sqrtmass = m.sqrt(masses[i])
        rvect = [float(x_vec.split()[2].replace(',','')),float(y_vec.split()[2].replace(',','')),float(z_vec.split()[2].replace(',',''))]
        Reigvects.append(rvect)
        ivect = [float(x_vec.split()[3].replace(',','')),float(y_vec.split()[3].replace(',','')),float(z_vec.split()[3].replace(',',''))]
        Ieigvects.append(ivect)                                                 # form the x,y,z on either real/imaginary eigenvector for each atom. 
    return (Reigvects,Ieigvects)                                            

# Forms a vectors of -1s, 0s and 1s from a yaml real eigenvector
def nonzerovector(eigvects):
    nonzeros = []
    for vec in eigvects:
        nonzero = [round(float(vec[0]),8)/abs(divzero(round(float(vec[0]),8))),round(float(vec[1]),8)/abs(divzero(round(float(vec[1]),8))),round(float(vec[2]),8)/abs(divzero(round(float(vec[2]),8)))]
        nonzeros = nonzeros + nonzeros          # WHY NOT APPEND????? Plus will lose the array inside the array. 
    return nonzeros

# Produces a number safe to divide by without ZeroDivisionError (change tol to app value for non-doubles)
def divzero(number,tol=1e-15):
    if abs(number) <= tol:
        number = tol
    return number

# Produces a list of maximum displacements from complex phonon eigenvectors
def modvects(rvects,ivects):
    pos_digit = 16               # Number of decimal places in atom posns
    string = "{0:." + str(pos_digit) + "f}"
    i = 0
    recipsq2 = 1.0/m.sqrt(2.0)
    modeigvects = []
    phshifts = []
    for rvect in rvects:
        ivect = ivects[i]
        X02 = rvect[0]**2+ivect[0]**2
        Y02 = rvect[1]**2+ivect[1]**2
        Z02 = rvect[2]**2+ivect[2]**2
        psix = m.atan(ivect[0]/divzero(rvect[0]))
        psiy = m.atan(ivect[1]/divzero(rvect[1]))
        psiz = m.atan(ivect[2]/divzero(rvect[2]))
        A = X02*m.cos(2*psix)+Y02*m.cos(2*psiy)+Z02*m.cos(2*psiz)
        B = X02*m.sin(2*psix)+Y02*m.sin(2*psiy)+Z02*m.sin(2*psiz)
        R = m.sqrt(A**2+B**2)
        atan = m.atan(abs(B/divzero(A)))
        if A>=0 and B>=0:
            alpha = atan
        elif A>=0 and B<0:
            alpha = -atan
        elif A<0 and B>=0:
            alpha = m.pi-atan
        elif A<0 and B<0:
            alpha = atan-m.pi
        rmax = recipsq2*m.sqrt(X02+Y02+Z02+R)
        theta = -0.5*alpha
        modeigvects.append(rmax)
        phshifts.append(theta)
        i = i + 1
    return (modeigvects,phshifts)

# Outputs the sign of a number only
def sign(number):
    sign = number/abs(divzero(number))
    return sign

# Shifts phase of each ion displacement in an eigenvector
def shifteigvect(rvects,ivects,shift):
    pos_digit = 16               # Number of decimal places in atom posns
    string = "{0:." + str(pos_digit) + "f}"
    i = 0
    recipsq2 = 1.0/m.sqrt(2.0)
    newrvects = []
    newivects = []
    for rvect in rvects:
        ivect = ivects[i]
        X0 = sign(rvect[0])*m.sqrt(rvect[0]**2+ivect[0]**2)
        Y0 = sign(rvect[1])*m.sqrt(rvect[1]**2+ivect[1]**2)
        Z0 = sign(rvect[2])*m.sqrt(rvect[2]**2+ivect[2]**2)
        psix = m.atan(ivect[0]/divzero(rvect[0]))
        psiy = m.atan(ivect[1]/divzero(rvect[1]))
        psiz = m.atan(ivect[2]/divzero(rvect[2]))
        newrvect = [X0*m.cos(psix+shift),Y0*m.cos(psiy+shift),Z0*m.cos(psiz+shift)]
        newivect = [X0*m.sin(psix+shift),Y0*m.sin(psiy+shift),Z0*m.sin(psiz+shift)]
        newrvects.append(newrvect)
        newivects.append(newivect)
        i = i + 1
    return (newrvects,newivects)

# Calculates the point in the phase at which the ion with max amplitude displays that max amplitude
def phaseshift(eigvects,ivects):
    (modeigvects,phshifts) = modvects(eigvects,ivects)     # Returns list of abs val of each complex eigenvector and phase shift
    eigmax = max(modeigvects)                              # Largest rmax of eigenvectors
    maxphaseshift = phshifts[modeigvects.index(eigmax)]    # Phase shift required to see largest amplitude displacement
    (eigvects,ivects) = shifteigvect(eigvects,ivects,maxphaseshift)
    return (eigvects,ivects) 

# Appends vertically one numpy array (array2) to a first (array1) even if array1 is empty
def npvstack(array1,array2):
    if len(array1) == 0:
        array1 = array2
    else:
        array1 = np.vstack((array1,array2))
    return array1

# Prints a castep .cell file with absolute atomic positions
def createcell(printcell,latvecs,elements,sumposns):
    cell = []
    pos_digit = 16               # Number of decimal places in atom posns
    string = "{0:." + str(pos_digit) + "f}"
    cell.append('%BLOCK LATTICE_CART\n')
    cell.append(str(string.format(latvecs[0][0]))+'\t'+str(string.format(latvecs[0][1]))+'\t'+str(string.format(latvecs[0][2]))+'\n')
    cell.append(str(string.format(latvecs[1][0]))+'\t'+str(string.format(latvecs[1][1]))+'\t'+str(string.format(latvecs[1][2]))+'\n')
    cell.append(str(string.format(latvecs[2][0]))+'\t'+str(string.format(latvecs[2][1]))+'\t'+str(string.format(latvecs[2][2]))+'\n')
    cell.append('%ENDBLOCK LATTICE_CART\n')
    cell.append('\n')
    cell.append('%BLOCK POSITIONS_ABS\n')
    i = 0
    for element in elements:
        cell.append(element+'\t'+str(string.format(sumposns[i][0]))+'\t'+str(string.format(sumposns[i][1]))+'\t'+str(string.format(sumposns[i][2]))+'\n')
        i = i+1
    cell.append('%ENDBLOCK POSITIONS_ABS\n')
    f = open(printcell,'w')
    f.writelines( cell )
    f.close()
    return

# Forms a vector from a selection of entries in a larger array
def makevect(array,arrayinds):
    vector = []
    for i in range(0,len(arrayinds)):
        ai = int(arrayinds[i])
        for j in range(0,3):
            vector.append(int(array[ai*3+j]))
    return vector

# Find the index of a particular instance of a particular character in a string (note 1st instance = 1)
def findchar(string,specchar,instance=1):
    i=0
    for j in range(0,len(string)):
        character = string[j]
        if character == specchar:
            i=i+1
            if i==int(instance):
                break
    return j

# Confine float to a 0.0<=x<1.0 range
def confine(x):
    x = float(x)
    if x >= 1.000:
        y = x - int(x)
    elif x < 0.000:
        y = 1.0 - (abs(x) - int(abs(x)))
    else:
        y = x
    return y

# Confine 3D vector in fractional coordinates to a positive 1x1x1 unit cell
def confinevec(vec):
    newvec = []
    for element in vec:
        newvec.append(confine(element))
    return newvec

# Calculate the distance between two 3D fractional vectors accounting for periodic boundary conditions
def pbcdist(vec1,vec2):
    difvec = []
    for j in range(0,len(vec1)):
        d = abs(float(vec1[j]) - float(vec2[j]))
        if d > 0.5 and d <= 1.0: 
            d = 1.0 - d
        elif d > 1.0 and d <= 1.5:
            d = d - 1.0
        elif d > 1.5 and d <= 2.0:
            d = 2.0 - d
        elif d > 2.0 and d <= 2.5:
            d = d - 2.0
        elif d > 2.5: # Hopefully this would never happen e.g. everything should be in a range -1 to 1.5 
            print "Warning: distance between two coordinates greater than 2.5 lattice parameters"
        difvec.append(d)
    dist = np.linalg.norm(np.array(difvec))
    return dist

# Calculate the distance between two 3D fractional vectors accounting for periodic boundary conditions matching ion labels
def pbcdist2(vec1,vec2,label1,label2):
    difvec = []
    if label1 in label2:
        for j in range(0,len(vec1)):
            d = abs(float(vec1[j]) - float(vec2[j]))
            if d > 0.5 and d <= 1.0: 
                d = 1.0 - d
            elif d > 1.0 and d <= 1.5:
                d = d - 1.0
            elif d > 1.5 and d <= 2.0:
                d = 2.0 - d
            elif d > 2.0 and d <= 2.5:
                d = d - 2.0
            elif d > 2.5: # Hopefully this would never happen e.g. everything should be in a range -1 to 1.5 
                print "Warning: distance between two coordinates greater than 2.5 lattice parameters"
            difvec.append(d)
        dist = np.linalg.norm(np.array(difvec))
    else:
        dist = 10
    return dist
