# -*- coding: utf-8 -*-
import numpy as np
import subprocess,os

def checkCEA(compiler='gfortran'):
    """
    Check that all the fortran CEA codes have their 
    corresponding executables:
    """
    cwd = os.getcwd()
    os.chdir('CEA+FORTRAN')
    if not os.path.exists('FCEA2'):
        subprocess.Popen(compiler+' cea2.f -o FCEA2',\
                         shell = True).wait()
    if not os.path.exists('b1b2b3'):
        subprocess.Popen(compiler+' b1b2b3.f -o b1b2b3',\
                         shell = True).wait()
    if not os.path.exists('syntax'):
        subprocess.Popen(compiler+' syntax.f -o syntax',\
                         shell = True).wait()
    if not os.path.exists('thermo.lib'):
        p = subprocess.Popen(['./FCEA2'], stdout=subprocess.PIPE, \
                             stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_data = p.communicate(input='thermo')[0]
    if not os.path.exists('trans.lib'):
        p = subprocess.Popen(['./FCEA2'], stdout=subprocess.PIPE, \
                             stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_data = p.communicate(input='trans')[0]

    os.chdir(cwd)
    #print('CEA up and running!')

def runCEA(filename):
    filename_no_extension = (filename.split('.inp')[0]).split('/')[1]
    subprocess.Popen('cp '+filename+' CEA+FORTRAN/.',shell = True).wait()
    cwd = os.getcwd()
    os.chdir('CEA+FORTRAN')
    p = subprocess.Popen(['./FCEA2'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_data = p.communicate(input=filename_no_extension.encode())[0]
    subprocess.Popen('mv '+filename_no_extension+'.out ../CEAoutput/.',shell = True).wait() 
    subprocess.Popen('rm '+filename_no_extension+'.inp',shell = True).wait()
    os.chdir(cwd)

def read_abundances(filename):
    """
    This code reads abundances from a file containing Z, name, A(X) and error on this 
    number, where A(X) = log N(X)/N(H) + 12, and N(X) is the number of atoms of element 
    X. The code returns:

       Z            Atomic number of the element X
       Name         Name of the element X
       10**A        10 to the (log N(X)/N(H) - 12), where N(X) is the number of atoms of element X,
                    and N(H) is the number of hydrogen atoms (thus, this takes 1e12 for hydrogen 
                    by definition).
    """
    Z = np.array([])
    name = np.array([])
    logN = np.array([])

    f = open(filename,'r')
    while True:
        line = f.readline()
        if line != '':
            if line[0] != '#':
                data = line.split()
                Z = np.append(Z,np.double(data[0]))
                name = np.append(name,data[1])
                logN = np.append(logN,np.double(data[2]))
        else:
            break
    f.close()
    return Z,name,10**(logN)
    
def calcCEA(T,P,name,N,only_consider_these,prefix='benchmark',ions = False):
    """
    This function generates input files for CEA, runs them, and 
    spits out the equilibrium compositions:
    """

    if not os.path.exists('CEAoutput'):
        os.mkdir('CEAoutput')
    f = open('CEAoutput/'+prefix+'_'+str(T)+'_'+str(P)+'.inp','w')
    f.write('# problem dataset: transpec\n')
    if ions:
        f.write('problem tp ions\n')
    else:
        f.write('problem tp\n')
    f.write('  p(bar) = '+str(P)+'\n')
    f.write('  t(k) = '+str(T)+'\n')
    f.write('reac\n')
    for i in range(len(name)):
        f.write('  na '+name[i]+' moles = '+str(N[i]).upper()+'\n')
    f.write('only '+only_consider_these)
    f.write('output calories siunits trace 1e-20\n')
    f.write('end')
    f.close()
    runCEA('CEAoutput/'+prefix+'_'+str(T)+'_'+str(P)+'.inp')

def readCEA(T,P,prefix='benchmark'):
    f = open('CEAoutput/'+prefix+'_'+str(T)+'_'+str(P)+'.out','r')
    species = []
    moles = []
    while True:
        line = f.readline()
        if 'MOLE FRACTIONS' in line:
            f.readline()
            break
        if line[:21] == ' CALCULATIONS STOPPED' or line == '':
            f.close()
            return species,moles
    while True:
        line = f.readline()
        vec = line.split()
        if len(vec) == 0:
            break
        elif line[:21] == ' CALCULATIONS STOPPED':
            break
        else:
            if '*' not in vec[1]:
                species.append(vec[0])
                if '-' in vec[1]:
                    num,exp = vec[1].split('-')
                    moles.append(np.double(num+'e-'+exp))
                else:
                    moles.append(np.double(vec[1]))
    f.close()
    return species,moles

def check_elements(name,only_consider_these):
    idx = []
    for i in range(len(name)):
        c_name = name[i]
        length = len(c_name)
        for j in range(len(only_consider_these)-length):
            if only_consider_these[j:j+length] == c_name and only_consider_these[j-1]!='(' \
                            and only_consider_these[j-1]!='I':
                if length == 1:
                    if not only_consider_these[j+1] == only_consider_these[j+1].lower():
                        if c_name == 'C':
                            if only_consider_these[j+1] == 'L':
                                if only_consider_these[j+1:j+3] == 'Li' or\
                                   only_consider_these[j+1:j+3] == 'La' or\
                                   only_consider_these[j+1:j+3] == 'Lu':
                                       idx.append(i)
                                       break
                            else:
                                idx.append(i)
                                break
                        else:
                            idx.append(i)
                            break
                    else:
                        if only_consider_these[j+1] == ' ':
                            idx.append(i)
                            break
                else:
                    idx.append(i)
                    break
    return idx

def readPropCEA(prop,t,p,fuel,oxid,of,prefix='benchmark'):
    filename = 'CEAoutput/'+prefix+'_'+str(t)+'_'+str(p)+'_'+str(fuel)+'_'+oxid+'_'+str(of)+'.out'
    filename = filename.replace('(','').replace(')','')
    f = open(filename,'r')

    if prop == 'gama':
        while True:
            line = f.readline()
            if 'GAMMAs' in line:
                gama = float(line.split(' ')[-1][:-1])
                break
        f.close()
        return gama
    if prop == 'cp':
        while True:
            line = f.readline()
            if 'WITH FROZEN REACTIONS' in line:
                f.readline()
                line = f.readline()
                cp = round(float(line.split(' ')[-1][:-1])*1000,1)
                break
        f.close()
        return cp

def calcPropCEA(t,p,fuel,oxid,of,prefix='benchmark'):
    filename = 'CEAoutput/'+prefix+'_'+str(t)+'_'+str(p)+'_'+str(fuel)+'_'+oxid+'_'+str(of)+'.inp'
    filename = filename.replace('(','').replace(')','')

    if os.path.isfile(filename):
        return
        
    if not os.path.exists('CEAoutput'):
        os.mkdir('CEAoutput')
    f = open(filename,'w')
    f.write('problem case=1 tp\n')
    f.write('  p(bar) = '+str(p)+'\n')
    f.write('  t(k) = '+str(t)+'\n')
    f.write('reac\n')
    f.write('  fuel='+str(fuel)+' wt=1\n')
    f.write('  oxid='+str(oxid)+' wt='+str(of)+'\n')
    f.write('output\n')
    f.write('  siunits short, transport\n')
    f.write('  t gam cp\n')
    f.write('end\n')
    f.close()
    runCEA(filename)

def calcPropStagnationCEA(p,fuel,oxid,of,prefix='benchmark'):
    filename = 'CEAoutput/'+prefix+'_'+str(p)+'_'+str(fuel)+'_'+oxid+'_'+str(of)+'_rkt.inp'
    filename = filename.replace('(','').replace(')','')

    if os.path.isfile(filename):
        return

    if not os.path.exists('CEAoutput'):
        os.mkdir('CEAoutput')
    f = open(filename,'w')
    f.write('problem    o/f='+str(of)+',\n')
    f.write('    rocket  fac   ac/at=1  tcest,k=3800\n')
    f.write('  p,bar='+str(p)+',\n')
    f.write('react  \n')
    f.write('  fuel='+str(fuel)+' \n')
    f.write('  oxid='+str(oxid)+' \n')
    f.write('output  transport \n')
    f.write('end\n')
    f.close()
    runCEA(filename)

def readPropStagnationCEA(prop,p,fuel,oxid,of,prefix='benchmark'):
    filename = 'CEAoutput/'+prefix+'_'+str(p)+'_'+str(fuel)+'_'+oxid+'_'+str(of)+'_rkt.out'
    filename = filename.replace('(','').replace(')','')
    f = open(filename,'r')

    if prop == 't':
        while True:
            line = f.readline()
            if 'T, K' in line:
                t = float(line.split('  ')[-3])
                break
        f.close()
        return t
    if prop == 'cp':
        while True:
            line = f.readline()
            if 'WITH FROZEN REACTIONS' in line:
                f.readline()
                line = f.readline()
                cp = round(float(line.split('  ')[-3])*1000,1)
                break
        f.close()
        return cp
    if prop == 'pr':
        while True:
            line = f.readline()
            if 'WITH FROZEN REACTIONS' in line:
                f.readline()
                f.readline()
                f.readline()
                line = f.readline()
                pr = float(line.split('  ')[-3])
                break
        f.close()
        return pr
    if prop == 'mi':
        while True:
            line = f.readline()
            if 'VISC,MILLIPOISE' in line:
                mi = float(line.split('  ')[-3])/1000
                break
        f.close()
        return mi

checkCEA()
