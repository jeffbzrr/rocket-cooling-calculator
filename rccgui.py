'''
Interface gŕafica da biblioteca librcc.py
Jefferson Bezerra
'''

import PySimpleGUI as sg
import sys
import json
import librcc as rcc
import csv
import pyCEA

def output_files(data_in, data_out):
    data = []
    data_row = ['k_c', 'CCH', 'R_c', 'L', 'eta_o', 'r2', 'Taw', 'm', 'Tc', 'Ae/At', 'size', 'V_c', 'Twc', 'Aa', 'z', 'mi_c', 'M', 'zt', 'R_g', 'eta_f', 'Pr_c', 'Rt', 'At', 'gama', 'p_drop', 'N', 'cp_c', 'Q', 'error_code', 'Twg', 'q', 'Re_c', 'Ae', 'mi_s', 'h_g', 'deltap', 'D_h', 'T_static', 'r3', 'p_static', 'r1', 'R', 'R_w', 'f', 'hl', 'ro', 'cp', 'h_c', 'FT', 'CCW', 'Atotal', 't']
    data.append(data_row)
    for i in range(0,data_out['size']):
        data_row = [data_out['k_c'][i], data_out['CCH'][i], data_out['R_c'][i], data_out['L'][i], data_out['eta_o'][i], data_out['r2'][i], data_out['Taw'][i], data_out['m'][i], data_out['Tc'][i], data_out['Ae/At'][i], data_out['size'], data_out['V_c'][i], data_out['Twc'][i], data_out['Aa'][i], data_out['z'][i], data_out['mi_c'][i], data_out['M'][i], data_out['zt'], data_out['R_g'][i], data_out['eta_f'][i], data_out['Pr_c'][i], data_out['Rt'], data_out['At'], data_out['gama'][i], data_out['p_drop'], data_out['N'][i], data_out['cp_c'][i], data_out['Q'][i], data_out['error_code'], data_out['Twg'][i], data_out['q'][i], data_out['Re_c'][i], data_out['Ae'][i], data_out['mi_s'][i], data_out['h_g'][i], data_out['deltap'][i], data_out['D_h'][i], data_out['T_static'][i], data_out['r3'][i], data_out['p_static'][i], data_out['r1'][i], data_out['R'][i], data_out['R_w'][i], data_out['f'][i], data_out['hl'][i], data_out['ro'][i], data_out['cp'][i], data_out['h_c'][i], data_out['FT'][i], data_out['CCW'][i], data_out['Atotal'][i], data_out['t']]
        data.append(data_row)

    with open('{}_data_output.csv'.format(data_in['motor_name']), mode='w', encoding='utf-8') as data_file:
        csv_writer = csv.writer(data_file, delimiter=',')
        csv_writer.writerows(data)

menu_def = [['Arquivo', ['Abrir', 'Salvar como...', 'Sair']],     
['Ajuda', 'Sobre...'], ]   
 
frame_calcular = [
[sg.Text('Tolerânica', size=(18, 1)), sg.InputText('1e-6', size=(8, 1), key='tol')],
[sg.Text('Máximo de iterações', size=(18, 1)), sg.InputText('100', size=(8, 1), key='max_iterations')],
[sg.Text(' ')],
[sg.Submit('Calcular')]
]

frame_otimizar = [
[sg.Text('Tolerânica', size=(23, 1)), sg.InputText('1e-6' ,size=(8, 1), key='tol_opt'),
sg.Text('Máximo de iterações', size=(19, 1)), sg.InputText('4000', size=(8, 1), key='max_iterations_opt')],
[sg.Text('Espessura mínima de aleta', size=(23, 1)), sg.InputText(size=(8, 1), key='FT_min'),
sg.Text('Altura mínima do canal', size=(19, 1)), sg.InputText(size=(8, 1), key='CCH_min')],
[sg.Text('Largura mínima do canal', size=(23, 1)), sg.InputText(size=(8, 1), key='CCW_min'),
sg.Text('Passo das iterações', size=(19, 1)), sg.InputText('0.2' ,size=(8, 1), key='w')],
[sg.Submit('Otimizar')]
]

layout = [
    [sg.Menu(menu_def, tearoff=True)],
    [sg.Text('Nome do motor', size=(22, 1)), sg.InputText(size=(22, 1), key='motor_name'),
    sg.Text('Geometria', size=(15, 1)), sg.InputText(size=(22, 1), key='geometry_path'), sg.FileBrowse(file_types=(('Arquivos RCC', '*.csv'),),)],
    [sg.Text('Posição dos canais (m)', size=(22, 1)), sg.InputText(size=(22, 1), key='channel_number_dim'),
    sg.Text('Número de canais', size=(15, 1)), sg.InputText(size=(22, 1), key='channel_number_qt')],
    [sg.Text('Espessura da parede interna (m)', size=(28, 1)), sg.InputText(size=(8, 1), key='IWT'),
    sg.Text('Dimensão constante', size=(17, 1)), sg.InputCombo(('Largura da aleta', 'Largura do canal'), size=(21, 1), key='constant_dim')],   
    [sg.Text('Altura do canal (m)', size=(16, 1)), sg.InputText(size=(8, 1), key='CCH'),
    sg.Text('Espessura da aleta (m)', size=(19, 1)), sg.InputText(size=(8, 1), key='FT'), sg.Text('Largura do canal (m)', size=(18, 1)), sg.InputText(size=(8, 1), key='CCW')],
    [sg.Text('Temperaturas iniciais:')],
    [sg.Text('Refrigerante (K)', size=(25, 1)), sg.InputText('303', size=(8, 1), key='Tc_primary'),
    sg.Text('Base da aleta (K)', size=(23, 1)), sg.InputText('600', size=(8, 1), key='Twc_primary')],
    [sg.Text('Parede interna da câmara (K)', size=(25, 1)), sg.InputText('1000', size=(8, 1), key='Twg_primary'),
    sg.Text('Gases da combustão (K)', size=(23, 1)), sg.InputText('3000', size=(8, 1), key='Taw_primary')],
    [sg.Text('Refrigerante'), sg.InputCombo(('RP-1', 'C2H5OH(L)'), size=(10, 1), key='coolant'),
    sg.Text('Combustível'), sg.InputCombo(('RP-1', 'C2H5OH(L)'), size=(10, 1), key='fuel'),
    sg.Text('Oxidante'), sg.InputCombo(('O2(L)', 'N2O'), size=(10, 1), key='oxidizer')],
    [sg.Text('Razão de mistura', size=(19, 1)), sg.InputText(size=(8, 1), key='of'),
    sg.Text('Vazão mássica do combustível (kg/s)', size=(31, 1)), sg.InputText(size=(8, 1), key='m._c')],
    [sg.Text('Condutividade térmica da câmara e aletas (W m^-1 K^-1)', size=(47, 1)), sg.InputText(size=(8, 1), key='k_w'),
    sg.Text('Pressão de estagnação (MPa)', size=(25, 1)), sg.InputText(size=(8, 1), key='p0')],
    [sg.Text('Rugosidade dos canais', size=(19, 1)), sg.InputText(size=(8, 1), key='e'),
    sg.Text('Fator de segurança da parede', size=(25, 1)), sg.InputText('1', size=(8, 1), key='fs_wall'),
    sg.Text('Limite de escoamento (MPa)', size=(24, 1)), sg.InputText(size=(8, 1), key='sigmae')],
    [sg.Frame('Dados para cálculo', frame_calcular), sg.Frame('Dados para otimização', frame_otimizar)]]

window = sg.Window('RCC GUI', layout)

while True:
    event, values = window.Read()
    if event is None or event == 'Sair':
        break
    if event == 'Salvar como...':
        save_path = sg.PopupGetFile('Salvar como...', save_as=True, file_types=(('Arquivos RCC', '*.rcc'),),)
        with open(save_path, 'w') as outfile:
            json.dump(values, outfile)
    if event == 'Abrir':
        load_path = sg.PopupGetFile('Abrir', file_types=(('Arquivos RCC', '*.rcc'),),)
        if load_path != None and load_path != '':
            json_file = open(load_path)
            values = json.load(json_file)
            del values['0']
            del values['Browse']
            for v in values:
                window.Element(v).Update(values[v])
    if event == 'Sobre...':
        sg.Popup('Versão 1.0', 'https://github.com/jeffersonmsb')
    if event == 'Calcular' or event == 'Otimizar':
        sg.PopupNoWait('Acompanhe a execução através da janela terminal.', 'Feche a janela terminal para interromper a execução.',title='Calculando...', auto_close=True)
        entrada = values.copy()
        entrada.pop(0, None)
        entrada.pop('Browse', None)
        entrada['e'] = float(entrada['e'])
        entrada['sigmae'] = float(entrada['sigmae'])*1e6
        entrada['max_iterations'] = int(entrada['max_iterations'])
        entrada['max_iterations_opt'] = int(entrada['max_iterations_opt'])
        entrada['Twc_primary'] = float(entrada['Twc_primary'])
        entrada['Tc_primary'] = float(entrada['Tc_primary'])
        entrada['Twg_primary'] = float(entrada['Twg_primary'])
        entrada['Taw_primary'] = float(entrada['Taw_primary'])
        entrada['m._c'] = float(entrada['m._c'])
        entrada['FT_min'] = float(entrada['FT_min'])
        entrada['CCW_min'] = float(entrada['CCW_min'])
        entrada['CCW'] = float(entrada['CCW'])
        entrada['CCH_min'] = float(entrada['CCH_min'])
        entrada['CCH'] = float(entrada['CCH'])
        entrada['IWT'] = float(entrada['IWT'])
        entrada['fs_wall'] = float(entrada['fs_wall'])
        entrada['k_w'] = float(entrada['k_w'])
        entrada['tol_opt'] = float(entrada['tol_opt'])
        entrada['tol'] = float(entrada['tol'])
        entrada['w'] = float(entrada['w'])
        entrada['of'] = float(entrada['of'])
        entrada['p0'] = float(entrada['p0'])*1e6
        if entrada['constant_dim'] == 'Largura da aleta':
            entrada['dim_constant'] = 'FT'
            entrada['FT'] = float(entrada['FT'])
        else:
            entrada['dim_constant'] = 'CCW'
            entrada['CCW'] = float(entrada['CCW'])
        entrada['channel_number_dim'] = entrada['channel_number_dim'].split(',')
        entrada['channel_number_dim'] = list(map(float, entrada['channel_number_dim']))
        entrada['channel_number_qt'] = entrada['channel_number_qt'].split(',')
        entrada['channel_number_qt'] = list(map(int, entrada['channel_number_qt']))
        
        entrada['channel_number'] = []
        for i in range(0,len(entrada['channel_number_qt'])):
            entrada['channel_number'].append([entrada['channel_number_dim'][i],entrada['channel_number_qt'][i]])

        saida = {}
        if event == 'Calcular':
            rcc.iteration(entrada, saida)
        else:
            rcc.optimize_channel2(entrada, saida)

        saida['t'] = rcc.calc_wall_thickness(entrada['p0'], entrada['geometry_path'], entrada['sigmae'], entrada['fs_wall'])
        output_files(entrada, saida)
        sg.PopupNoWait('Finalizado', 'Arquivo {}_data_output.csv criado.'.format(entrada['motor_name']),title='Finalizado')

window.Close()