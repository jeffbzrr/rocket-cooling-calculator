'''
Biblioteca para calculo de refrigeracao regenerativa em motores foguetes bi propelentes
Jefferson Bezerra
https://github.com/jeffersonmsb/rocket-cooling-calculator
'''

import csv
import numpy as np
import math
import pyCEA
from scipy import optimize
import os
import subprocess

def geometry(data_in, data_out):
    with open(data_in['geometry_path']) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        data_out['geometry'] = list(csv_reader)

        data_out['size'] = len(data_out['geometry'])

        #Buscar geometria da garganta
        Rt = float(data_out['geometry'][0][1])
        zt = float(data_out['geometry'][0][0])
        for row in data_out['geometry']:
            if float(row[1]) < Rt:
                Rt = float(row[1])
                zt = float(row[0])
        data_out['Rt'] = Rt
        data_out['zt'] = zt
        data_out['At'] = np.pi*np.power(Rt,2)

        #Cálculo das razões de área
        data_out['r1'] = []
        data_out['r2'] = []
        data_out['r3'] = []
        data_out['Ae'] = []
        data_out['Ae/At'] = []
        data_out['z'] = []
        data_out['N'] = []
        data_out['CCH'] = []
        data_out['CCW'] = []
        data_out['FT'] = []
        n = 0

        for row in data_out['geometry']:
            A = np.pi*np.power(float(row[1]),2) 
            data_out['r1'].append(float(row[1]))
            r2 = float(row[1]) + data_in['IWT']
            data_out['r2'].append(r2)
            data_out['r3'].append(float(row[1]) + data_in['IWT'] + data_in['CCH'])
            data_out['Ae'].append(A)
            data_out['Ae/At'].append(A/data_out['At'])
            data_out['z'].append(float(row[0]))

            if float(row[0]) > data_in['channel_number'][n][0]:
                n = n + 1

            N = data_in['channel_number'][n][1]
            data_out['N'].append(N)

            data_out['CCH'].append(data_in['CCH'])
            
            if data_in['dim_constant'] == 'FT':
                data_out['FT'].append(data_in['FT'])
                aux = (2*np.pi*r2)/N - data_in['FT']
                if aux <= 0:
                    data_out['error_code'] = 1
                    return
                data_out['CCW'].append(aux)
            else:
                data_out['CCW'].append(data_in['CCW'])
                aux = (2*np.pi*r2)/N - data_in['CCW']
                data_out['FT'].append(aux)

        data_out['L'] = []
        for i in range(0, data_out['size']):
            if(i==0):
                A = 0.5*(data_out['z'][i+1]+data_out['z'][i]) - data_out['z'][i]
                B = 0.5*(data_out['r1'][i+1]+data_out['r1'][i]) - data_out['r1'][i]
                data_out['L'].append(math.sqrt(A**2 + B**2))
            else:
                if(i!=(data_out['size']-1)):
                    A = 0.5*(data_out['z'][i+1]+data_out['z'][i]) - 0.5*(data_out['z'][i]+data_out['z'][i-1])
                    B = 0.5*(data_out['r1'][i+1]+data_out['r1'][i]) - 0.5*(data_out['r1'][i]+data_out['r1'][i-1])
                    data_out['L'].append(math.sqrt(A**2 + B**2))
                else:
                    A = data_out['z'][i] - 0.5*(data_out['z'][i]+data_out['z'][i-1])
                    B = data_out['r1'][i] - 0.5*(data_out['r1'][i]+data_out['r1'][i-1])
                    data_out['L'].append(math.sqrt(A**2 + B**2))
    data_out['error_code'] = 0

def coolant_prop(coolant_name, prop_name, temperature):
    if coolant_name == 'RP-1':
        if temperature > 800:
            temperature = 800
        if temperature < 300:
            temperature = 300

        if prop_name == 'ro':
            return 820
        if prop_name == 'cp':
            return -2.82649e-3*temperature**2.0 + 6.77751e0*temperature - 2.45234e1 #BOYSAN
        if prop_name == 'k':
            return 9.64e-8*temperature**2-2.95e-4*temperature+0.261 #BOYSAN
        if prop_name == 'mi':
            return -1.46e-11*temperature**3+3.22e-8*temperature**2-2.39e-5*temperature+6E-3 #BOYSAN
    if coolant_name == 'C2H5OH(L)':
        if prop_name == 'ro':
            return 785.3
        if prop_name == 'cp':
            return 2570
        if prop_name == 'k':
            return 0.167
        if prop_name == 'mi':
            return 1.36e-3
    else:
        print('Coolant proprieties not found')
        return -1

def create_prop(data_in, data_out):
    data_out['Tc'] = data_out['size']*[data_in['Tc_primary']]
    data_out['Twg'] = data_out['size']*[data_in['Twg_primary']]
    data_out['Twc'] = data_out['size']*[data_in['Twc_primary']]
    data_out['Taw'] = data_out['size']*[data_in['Taw_primary']]

    data_out['cp_c'] = data_out['size']*[None]
    data_out['k_c'] = data_out['size']*[None]
    data_out['mi_c'] = data_out['size']*[None]
    data_out['Pr_c'] = data_out['size']*[None]
    data_out['gama'] = data_out['size']*[None]
    data_out['M'] = data_out['size']*[None]
    data_out['cp'] = data_out['size']*[None]
    data_out['R'] = data_out['size']*[None]
    data_out['h_g'] = data_out['size']*[None]
    data_out['Re_c'] = data_out['size']*[None]
    data_out['D_h'] = data_out['size']*[None]
    data_out['mi_s'] = data_out['size']*[None]
    data_out['h_c'] = data_out['size']*[None]
    data_out['Aa'] = data_out['size']*[None]
    data_out['Atotal'] = data_out['size']*[None]
    data_out['m'] = data_out['size']*[None]
    data_out['eta_f'] = data_out['size']*[None]
    data_out['eta_o'] = data_out['size']*[None]
    data_out['R_c'] = data_out['size']*[None]
    data_out['R_g'] = data_out['size']*[None]
    data_out['R_w'] = data_out['size']*[None]
    data_out['q'] = data_out['size']*[None]
    data_out['Q'] = data_out['size']*[None]
    data_out['f'] = data_out['size']*[None]
    data_out['ro'] = data_out['size']*[None]
    data_out['V_c'] = data_out['size']*[None]
    data_out['hl'] = data_out['size']*[None]
    data_out['deltap'] = data_out['size']*[None]
    data_out['T_static'] = data_out['size']*[None]
    data_out['p_static'] = data_out['size']*[6000000]

def calc_prop(data_in, data_out):

    data_in['p0_pyCEA'] = data_in['p0']/1e5 #Conversão de [Pa] para [bar]
    pyCEA.calcPropStagnationCEA(data_in['p0_pyCEA'],data_in['fuel'], data_in['oxidizer'],data_in['of'], data_in['motor_name'])
    T0 = pyCEA.readPropStagnationCEA('t', data_in['p0_pyCEA'], data_in['fuel'], data_in['oxidizer'], data_in['of'], data_in['motor_name'])
    cp0 = pyCEA.readPropStagnationCEA('cp', data_in['p0_pyCEA'], data_in['fuel'], data_in['oxidizer'], data_in['of'], data_in['motor_name']) 
    Pr0 = pyCEA.readPropStagnationCEA('pr', data_in['p0_pyCEA'], data_in['fuel'], data_in['oxidizer'], data_in['of'], data_in['motor_name'])
    mi0 = pyCEA.readPropStagnationCEA('mi', data_in['p0_pyCEA'], data_in['fuel'], data_in['oxidizer'], data_in['of'], data_in['motor_name'])

    Tc1 = data_in['Tc_primary']
    IWT = data_in['IWT']
    k_w = data_in['k_w']
    mponto_c = data_in['m._c']
    e = data_in['e']
    p0 = data_in['p0']

    Re_c = data_out['Re_c']
    N = data_out['N']
    mi_c = data_out['mi_c']
    CCW = data_out['CCW']
    FT = data_out['FT']
    D_h = data_out['D_h']
    mi_s = data_out['mi_s']
    Tc = data_out['Tc'] 
    Twg = data_out['Twg']
    Twc = data_out['Twc']
    Taw = data_out['Taw']
    h_c = data_out['h_c']
    k_c = data_out['k_c']
    Pr_c = data_out['Pr_c']
    Aa = data_out['Aa']
    L = data_out['L']
    r1 = data_out['r1']
    r2 = data_out['r2']
    r3 = data_out['r3']
    Atotal = data_out['Atotal']
    m = data_out['m']
    eta_f = data_out['eta_f']
    eta_o = data_out['eta_o']
    R_c = data_out['R_c']
    R_g = data_out['R_g']
    R_w = data_out['R_w']
    h_g = data_out['h_g']
    q = data_out['q']
    Q = data_out['Q']
    cp_c = data_out['cp_c']
    k_c = data_out['k_c']
    f = data_out['f']
    ro = data_out['ro']
    V_c = data_out['V_c']
    hl = data_out['hl']
    deltap = data_out['deltap']
    T_static = data_out['T_static']
    p_static = data_out['p_static']
    gama = data_out['gama']
    M = data_out['M']
    CCH = data_out['CCH']

    data_out['p_drop'] = 0
    
    def f_mach(M):
        A = 2/(data_out['gama'][i]+1)
        B = 1+(((data_out['gama'][i]-1)/2)*(M**2))
        C = (data_out['gama'][i]+1)/(data_out['gama'][i]-1)
        D = (data_out['Ae/At'][i]*M)**2
        return ( (A*B)**C-D )

    def f_coolebrook(f):
        return (1/(-2*math.log(((e/D_h[i])/3.7)+(2.51/(Re_c[i]*f**0.5)), 10))**2-f)

    for i in reversed(range(0,data_out['size'])):
        cp_c[i] = coolant_prop(data_in['coolant'], 'cp', Tc[i])
        k_c[i] = coolant_prop(data_in['coolant'], 'k', data_out['Tc'][i])
        data_out['mi_c'][i] = coolant_prop(data_in['coolant'], 'mi', data_out['Tc'][i])
        data_out['Pr_c'][i] = data_out['cp_c'][i]*data_out['mi_c'][i]/data_out['k_c'][i]

        pyCEA.calcPropCEA(data_out['Taw'][i] , data_in['p0_pyCEA'], data_in['fuel'], data_in['oxidizer'], data_in['of'], data_in['motor_name'])
        data_out['cp'][i] = pyCEA.readPropCEA('cp', data_out['Taw'][i], data_in['p0_pyCEA'], data_in['fuel'], data_in['oxidizer'], data_in['of'], data_in['motor_name'])
        #data_out['cp'][i] = -5.84399e-05*data_out['Taw'][i]**2.0 + 4.23454e-01*data_out['Taw'][i] + 1.29256e+03
        data_out['gama'][i] = 1.23854e-8*data_out['Taw'][i]**2 - 8.09028e-5*data_out['Taw'][i] + 1.34563
        #Gama para o L-75
        #data_out['gama'][i] = pyCEA.readPropCEA('gama', data_out['Taw'][i], data_in['p0_pyCEA'], data_in['fuel'], data_in['oxidizer'], data_in['of'], data_in['motor_name'])
        data_out['R'][i] = (data_out['cp'][i]*(1 - 1/data_out['gama'][i]))
        mponto = data_in['p0']*data_out['At']*((data_out['gama'][i]/(data_out['R'][i]*T0))*(2/(data_out['gama'][i]+1))**((data_out['gama'][i]+1)/(data_out['gama'][i]-1)))**0.5
        c = (data_in['p0']*data_out['At'])/mponto
        if(data_out['z'][i] > data_out['zt']):
            a = 1
            b = 25
        else:
            a = 0
            b = 1
        data_out['M'][i] = optimize.bisect(f_mach, a, b, rtol=8.881784197001252e-16)
        aux1 = 1 + ((data_out['gama'][i]-1)/2)*data_out['M'][i]**2
        sigma = ((data_out['Twg'][i]/(2*T0))*aux1+0.5  )**-0.68  *  aux1**-0.12
        data_out['h_g'][i] = (  0.026  *  ((mi0/(2*data_out['Rt']))**0.2)  *  (cp0/(Pr0**0.6))  *  (data_in['p0']/c)**0.8  *  (data_out['At']/data_out['Ae'][i])**0.9  *  sigma  )

        D_h[i] = (4*CCW[i]*CCH[i])/(2*(CCW[i]+CCH[i]))
        Re_c[i] = (4*mponto)/(N[i]*mi_c[i]*2*(CCW[i]+CCH[i]))

        mi_s[i] = coolant_prop(data_in['coolant'], 'mi', Twc[i])
        h_c[i] = ((k_c[i]/D_h[i]) * 0.027 * Re_c[i]**0.8 * Pr_c[i]**(1/3) * (mi_c[i]/mi_s[i])**0.14 )

        Aa[i] = (2*CCH[i]*L[i])
        Atotal[i] = (N[i]*Aa[i] + L[i]*(2*math.pi*r2[i]-N[i]*FT[i]))
        m[i] = math.sqrt((2*h_c[i])/(k_c[i]*FT[i]))
        eta_f[i] = (math.tanh(m[i]*CCH[i])/(m[i]*CCH[i]))
        eta_o[i] = 1-((N[i]*Aa[i]*(1-eta_f[i])) / Atotal[i])

        R_g[i] = (1/(2*math.pi*r1[i]*L[i]*h_g[i]))

        R_w[i] = (math.log(r2[i]/r1[i]) / (2*math.pi*L[i]*k_w))

        R_c[i] = (1 / (eta_o[i]*h_c[i]*Atotal[i]))

        q[i] = ((Taw[i] - Tc[i]) / (R_g[i] + R_w[i] + R_c[i]))
        Q[i] =  ( q[i]/(2*math.pi*r1[i]*L[i])/1000000 )

        aux = 0.5*(data_out['gama'][i] - 1)*data_out['M'][i]**2
        Taw[i] = (T0 * ((1 + Pr0**(1/3)*aux) / (1 + aux)))        
        Twg[i] = -R_g[i]*q[i]+Taw[i]
        Twc[i] = -q[i]*(R_g[i]+R_w[i])+Taw[i]

        lista = reversed(range( i,data_out['size']))
        Tc1 = 303
        for j in lista:
            Tc2 = (q[j] / (mponto_c*cp_c[j])) + Tc1
            Tc[j] = (Tc2+Tc1)/2
            Tc1 = Tc2


        p_static[i] = p0*(1+((gama[i]-1)/2)*M[i]**2)**-(gama[i]/(gama[i]-1))

        #Cálculo da perda de carga
        f[i] = optimize.bisect(f_coolebrook, 0.00001, 2, rtol=8.881784197001252e-16)
        ro[i] = coolant_prop(data_in['coolant'], 'ro', Tc[i])
        V_c[i] = mponto_c/(ro[i]*CCH[i]*CCW[i]*N[i])
        hl[i] = f[i]*((L[i]/D_h[i])/(V_c[i]**2/2))
        deltap[i] = ro[i]*hl[i]*N[i]
        data_out['p_drop'] += deltap[i]

        #Cálculo da temperatura estática e pressão estática
        T_static[i] = T0*(1+((gama[i]-1)/2)*M[i]**2)**-1
        
def iteration(data_in , data_out):
    geometry(data_in, data_out)
    if data_out['error_code'] != 0:
        print('CCW <= 0')
        return
    create_prop(data_in, data_out)
    for i in range(0,data_in['max_iterations']):
        print('Iteration {}'.format(i+1))
        calc_prop(data_in, data_out)
        if i==0:
            Tc_0 = sum(data_out['Q'])
            Twg_0 = sum(data_out['Twg'])
            Twc_0 = sum(data_out['Twc'])
            Taw_0 = sum(data_out['Taw'])

            Tc_prev = Tc_0
            Twg_prev = Twg_0
            Twc_prev = Twc_0
            Taw_prev = Taw_0    
        else:
            Tc = sum(data_out['Q'])
            Twg = sum(data_out['Twg'])
            Twc = sum(data_out['Twc'])
            Taw = sum(data_out['Taw'])

            Tc_L1 = abs(Tc-Tc_prev)/Tc_0
            Twg_L1 = abs(Twg-Twg_prev)/Twg_0
            Twc_L1 = abs(Twc-Twc_prev)/Twc_0
            Taw_L1 = abs(Taw-Taw_prev)/Taw_0

            Tc_prev = Tc
            Twg_prev = Twg
            Twc_prev = Twc
            Taw_prev = Taw

            if Tc_L1 <= data_in['tol'] and Twg_L1 <= data_in['tol'] and Twc_L1 <= data_in['tol'] and Taw_L1 <= data_in['tol']:
                break
    print('Total Iteration Temperature: ' + str(i+1))


def optimize_channel2(data_in, data_out):
    flag1 = False
    flag2 = False
    if data_in['dim_constant'] == 'FT':
        dim_const = 'FT'
        dim_var = 'CCW'
    else:
        dim_const = 'CCW'
        dim_var = 'FT'

    geometry(data_in, data_out)
    m = 0
    for i in range(0, data_out['size']):
        if data_out['r2'][i] < data_out['r2'][m]:
            m = i
    dim_max = (2*np.pi*data_out['r2'][m])/data_out['N'][m] - data_in[dim_var + '_min']

    if dim_max-data_in[dim_const + '_min'] <= 0:
        print('Maior dimensão geométrica é menor que dimensão mínima.')
        return

    dim = (dim_max+data_in[dim_const + '_min'])/2
    x = np.array([data_in['CCH'] , dim])
    data_in[dim_const] = dim
    iteration(data_in, data_out)
    Q = max(data_out['Q'])
    Q_prev = Q
    Q0 = Q
    
    w = data_in['w']

    for opt in range(0,data_in['max_iterations_opt']):
        grad = np.gradient(x)
        xn = x - w*grad

        if xn[0] <= data_in['CCH_min'] and flag1 == False:
            flag1 = True
            print('CCH_min')
        if (xn[1] <= data_in[dim_const+'_min'] or xn[1] >= dim_max) and flag2 == False:
            flag2 = True
            print(dim_const+' min or max')

        if flag1 == True:
            xn[0] = x[0]
        
        if flag2 == True:
            xn[1] = x[1]

        data_in['CCH'] = xn[0]
        data_in[dim_const] = xn[1]
        iteration(data_in, data_out)
        Q = max(data_out['Q'])
        if Q-Q_prev < 0:
            w=w*-1
            print(w)
            continue
        
        x = xn
        print('Opt #{} Q:{} CCH:{} {}:{}'.format(opt, Q, x[0], dim_const, x[1]))

        Q_diff = abs(Q-Q_prev)/Q0
        Q_prev = Q

        if Q_diff <= data_in['tol_opt']:
            break

def plot(data_out):
    data = []
    for i in range(0,data_out['size']):
        data_row = [ data_out['z'][i], data_out['Q'][i], data_out['Taw'][i], data_out['Twg'][i], data_out['Twc'][i], data_out['Tc'][i]]
        data.append(data_row)

    with open('rcc_plot_data.csv', mode='w', encoding='utf-8') as data_file:
        csv_writer = csv.writer(data_file, delimiter=',')
        csv_writer.writerows(data)
    
    p = subprocess.Popen("gnuplot \'rcc_plot_config.gnu\'", shell = True)
    os.waitpid(p.pid, 0)
    '''p2 = subprocess.Popen("ristretto \'temps.png\'", shell = True)
    os.waitpid(p2.pid, 1)'''

def calc_wall_thickness(p, path, sigmae, n=1):

    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        r = float(list(csv_reader)[0][1])

    def f_t(t):
        sigma1 = (p*r)/t
        sigma2 = (p*r)/(2*t)
        return math.sqrt(sigma1**2-sigma1*sigma2+sigma2**2)-sigmae

    return (optimize.bisect(f_t, 1e-8, 1, rtol=8.881784197001252e-16))*n