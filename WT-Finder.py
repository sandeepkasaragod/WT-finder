from pyteomics import mgf, pepxml, mass
from itertools import islice
import re
import os
#import b_y_with_mod
import numpy as np
import matplotlib.pyplot as plt
from os import path

def reverse_peptide_for_decoy(peptide, modifications):
    rev_pep = peptide[::-1][1:] + peptide[-1]
    reverse_mod = ""
    if len(modifications) > 0:
        for mods in modifications.split(';'):
            full_name = re.search('\((.*?)\)', mods).group(0)
            location = re.findall('\d+',re.search('(.*?)\(', mods).group(1))[0]
            alphabet = re.findall('[A-Z]', re.search('(.*?)\(', mods).group(1))[0]
            #print (rev_pep, alphabet + str(len(peptide)-int(location)) + full_name)
            if len(reverse_mod) == 0:
                reverse_mod = reverse_mod + alphabet + str(len(peptide)-int(location)) + full_name
            else:
                reverse_mod = reverse_mod + "; " + alphabet + str(len(peptide)-int(location)) + full_name
    else:
        reverse_mod = " "
    return rev_pep, reverse_mod

input_fl = 'Input_data/BT20_LysC.txt'
infile_path = '/homes/homedirs8/bunit/NGS/Sandeep/mydata/raw_files/BT20/LysC/'
output_file = open('WTF-BT20_test.txt', 'w')

modification_list = {}

for mods in open('static_files/PTM_Modification.txt'):
    split_mods = mods.split('\t')
    modification_list[split_mods[0]] = split_mods[1]

dicts_mod_pos = {}
dicts_mz = {} # Store MGF files
dicts_pepmass = {} #Stores precursor ion
dicts_scan = {} #Store scan number
dicts_charge = {}
dicts_rt = {}
pepmass = []

def generate_score_theo(list_mz, pep_len): #list of m/z, lenght of peptide
    calc = 0.0
    for i in range(len(list_mz)):
        calc = calc + (((i + 1) / pep_len) * list_mz[i])
    return calc

def generate_socre_matched(list_mz, list_idx, pep_len):
    calc = 0.0
    for i in range(len(list_mz)):
        calc = calc + (((list_idx[i] + 1) / pep_len) * list_mz[i])
    return calc

def read_mgf(infile):
    file_name = infile.split('/')[-1].split('.')[0]
    print (file_name)
    if file_name  not in dicts_mz:
        dicts_mz[file_name] = {}
        dicts_pepmass[file_name] = {}
        dicts_scan[file_name] = {}
        dicts_charge[file_name] = {}
        dicts_rt[file_name] = {}
    for i in open(infile):
        if "BEGIN IONS" in i.rstrip():
            pepmass = ""
            rt = ""
            mz = ""
            scan = ""
            charge = ""
        if "scan=" in i.rstrip():
            scan =  i.rstrip().split("scan=")[1].replace('"', '')
        if "PEPMASS" in i.rstrip():
            pepmass = i.rstrip().split("=")[1]
        if "CHARGE" in i.rstrip():
            charge = re.findall('\d+', i.rstrip())[0]
        if "RTINSECONDS" in i.rstrip():
            rt = int(float(i.split('=')[1].rstrip()))
            rts = str(float(i.split('=')[1].rstrip()))
            lst = []
        if i.rstrip()[0].isdigit():
            lst.append(i.rstrip())
        if "END IONS" in i.rstrip():
            if rt not in dicts_mz[file_name]:
                dicts_mz[file_name][rt] = [lst]
                
            else:
                dicts_mz[file_name][rt].append(lst)
                
            if rt not in dicts_pepmass[file_name]:
                dicts_pepmass[file_name][rt] = [pepmass]
                dicts_scan[file_name][rt] = [scan]
                dicts_charge[file_name][rt] = [charge]
                dicts_rt[file_name][rt] = [rts]                
            else:
                dicts_pepmass[file_name][rt].append(pepmass)
                dicts_scan[file_name][rt].append(scan)
                dicts_charge[file_name][rt].append(charge)
                dicts_rt[file_name][rt].append(rts)
                

def generate_label_raw(charge):
    ap = [('b' + str(x+1), 'y' + str(x+1)) for x in range(charge)]
    return sorted([x for y in ap for x in y])

def generate_label(charge):
    ap = [('b' + str(x+1), 'y' + str(x+1)) for x in range(charge)]
    b = [('b' + str(x+1)) for x in range(charge)]
    y = [('y' + str(x+1)) for x in range(charge)]
    reverse_b = (b[::-1])
    combine = reverse_b + y
    return [x for x in combine]

def fragments(peptide, types=('b', 'y'), maxcharge=1):
    """
    The function generates all possible m/z for fragments of types 
    `types` and of charges from 1 to `maxharge`.
    """
    for i in range(1, len(peptide)):
        for ion_type in types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(
                            peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    yield mass.fast_mass(peptide[i:], ion_type=ion_type, charge=charge)

def get_theo_spectra_with_mod(peptide, charge, ptm_mods):
    global dicts_mod_pos
    
    frags = fragments(peptide, ('b', 'y'), charge)
    
    ls = [x for x in frags] # Adding the fragments to list 
    keys = generate_label_raw(charge)  #genares the list like ['b1', 'y1', 'b2' ...................]
    keys_count = len(keys)
    b_and_y = {keys[index]: ls[index:len(ls):keys_count] for index in range(0,keys_count)}
    mod_pos = {}
    mod_pos_rev = {}
    if re.search('\((.*?)\)', ptm_mods):
        for string in ptm_mods.split(';'):
            if re.search('\((.*?)\)', string).group(1).strip() in modification_list:
                # storing actual modification location -1 to match the list items of dictioanry keys
                mod_pos[int(re.search('\d+', string).group(0).strip()) - 1] = [float(modification_list[re.search('\((.*?)\)', string).group(1)])] # for b ions
                mod_pos_rev[(len(peptide)-1) - (int(re.search('\d+', string).group(0).strip()) - 1)] = [float(modification_list[re.search('\((.*?)\)', string).group(1)])] #for y ions

    for m, n in b_and_y.items():
        a = 0.0
        if 'b' in m:
            for i in range(len(n)):
                if '1' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))
                        
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + mod_pos[i][0] + float(a)]
                            a+=mod_pos[i][0]
                        else:
                            dicts_mod_pos[m].append(n[i] + mod_pos[i][0] + float(a))
                            a+=mod_pos[i][0]
                        
                elif '2' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))                        
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + (mod_pos[i][0]/2) + float(a)]
                            a+=(mod_pos[i][0]/2)
                        else:
                            dicts_mod_pos[m].append(n[i] + (mod_pos[i][0]/2) + float(a))
                            a+=(mod_pos[i][0]/2)
                        
                elif '3' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))
                        
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + (mod_pos[i][0]/3) + float(a)]
                            a+=(mod_pos[i][0]/3)
                        else:
                            dicts_mod_pos[m].append(n[i] + (mod_pos[i][0]/3) + float(a))
                            a+=(mod_pos[i][0]/3)
                        
                elif '4' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))

                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + (mod_pos[i][0]/4) + float(a)]
                            a+=(mod_pos[i][0]/4)
                        else:
                            dicts_mod_pos[m].append(n[i] + (mod_pos[i][0]/4) + float(a))
                            a+=(mod_pos[i][0]/4)
      
                elif '5' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + (mod_pos[i][0]/5) + float(a)]
                            a+=(mod_pos[i][0]/5)
                        else:
                            dicts_mod_pos[m].append(n[i] + (mod_pos[i][0]/5) + float(a))
                            a+=(mod_pos[i][0]/5)
                            
                elif '6' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + (mod_pos[i][0]/6) + float(a)]
                            a+=(mod_pos[i][0]/6)
                        else:
                            dicts_mod_pos[m].append(n[i] + (mod_pos[i][0]/6) + float(a))
                            a+=(mod_pos[i][0]/6)
                            
                elif '7' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + (mod_pos[i][0]/7) + float(a)]
                            a+=(mod_pos[i][0]/7)
                        else:
                            dicts_mod_pos[m].append(n[i] + (mod_pos[i][0]/7) + float(a))
                            a+=(mod_pos[i][0]/7)
                            
                elif '8' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))
                            
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + (mod_pos[i][0]/8) + float(a)]
                            a+=(mod_pos[i][0]/8)
                        else:
                            dicts_mod_pos[m].append(n[i] + (mod_pos[i][0]/8) + float(a))
                            a+=(mod_pos[i][0]/8)
                            
                elif '9' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))
                            
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + (mod_pos[i][0]/9) + float(a)]
                            a+=(mod_pos[i][0]/9)
                        else:
                            dicts_mod_pos[m].append(n[i] + (mod_pos[i][0]/9) + float(a))
                            a+=(mod_pos[i][0]/9)
                            
                elif '10' in m:
                    if i not in mod_pos:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n[i] + (mod_pos[i][0]/10) + float(a)]
                            a+=(mod_pos[i][0]/10)
                        else:
                            dicts_mod_pos[m].append(n[i] + (mod_pos[i][0]/10) + float(a))
                            a+=(mod_pos[i][0]/10)
                

        #For y ion reverse the list and add the items 
        if 'y' in m:
            n_rev = [x for x in n[::-1]] # Reverse the m/z for y ions
            for i in range(len(n_rev)):
                if '1' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))
                            
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + mod_pos_rev[i][0] + float(a)]
                            a+=mod_pos_rev[i][0]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + mod_pos_rev[i][0] + float(a))
                            a+=mod_pos_rev[i][0]
                    
                elif '2' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + (mod_pos_rev[i][0]/2) + float(a)]
                            a+=(mod_pos_rev[i][0]/2)
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + (mod_pos_rev[i][0]/2) + float(a))
                            a+=(mod_pos_rev[i][0]/2)
                        
                elif '3' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))

                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + (mod_pos_rev[i][0]/3) + float(a)]
                            a+=(mod_pos_rev[i][0]/3)
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + (mod_pos_rev[i][0]/3) + float(a))
                            a+=(mod_pos_rev[i][0]/3)
                        
                elif '4' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + (mod_pos_rev[i][0]/4) + float(a)]
                            a+=(mod_pos_rev[i][0]/4)
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + (mod_pos_rev[i][0]/4) + float(a))
                            a+=(mod_pos_rev[i][0]/4)
                       
                elif '5' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + (mod_pos_rev[i][0]/5) + float(a)]
                            a+=(mod_pos_rev[i][0]/5)
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + (mod_pos_rev[i][0]/5) + float(a))
                            a+=(mod_pos_rev[i][0]/5)
                            
                elif '6' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))
                            
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + (mod_pos_rev[i][0]/6) + float(a)]
                            a+=(mod_pos_rev[i][0]/6)
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + (mod_pos_rev[i][0]/6) + float(a))
                            a+=(mod_pos_rev[i][0]/6)
                        
                elif '7' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + (mod_pos_rev[i][0]/7) + float(a)]
                            a+=(mod_pos_rev[i][0]/7)
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + (mod_pos_rev[i][0]/7) + float(a))
                            a+=(mod_pos_rev[i][0]/7)
                            
                elif '8' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))
                            
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + (mod_pos_rev[i][0]/8) + float(a)]
                            a+=(mod_pos_rev[i][0]/8)
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + (mod_pos_rev[i][0]/8) + float(a))
                            a+=(mod_pos_rev[i][0]/8)
                            
                elif '9' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + (mod_pos_rev[i][0]/9) + float(a)]
                            a+=(mod_pos_rev[i][0]/9)
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + (mod_pos_rev[i][0]/9) + float(a))
                            a+=(mod_pos_rev[i][0]/9)
                            
                elif '10' in m:
                    if i not in mod_pos_rev:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + float(a)]
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + float(a))
                    else:
                        if m not in dicts_mod_pos:
                            dicts_mod_pos[m] = [n_rev[i] + (mod_pos_rev[i][0]/10) + float(a)]
                            a+=(mod_pos_rev[i][0]/10)
                        else:
                            dicts_mod_pos[m].append(n_rev[i] + (mod_pos_rev[i][0]/10) + float(a))
                            a+=(mod_pos_rev[i][0]/10)
            
    b = dicts_mod_pos
    for a, c in b.items():
        if 'y' in a:
            b[a] = c[::-1]
    
    b['b1-NH3'] = [float(x)-17.02655 for x in b['b1']]
    b['y1-NH3'] = [float(x)-17.02655 for x in b['y1']]
    
    b['b1-H2O'] = [float(x)-18.01057 for x in b['b1']]
    b['y1-H2O'] = [float(x)-18.01057 for x in b['y1']]
    return b

def reverse_y_ions(dicts_theo_mz):
    dicts_theo_mz_converted = {}
    for k, v in dicts_theo_mz.items():
        if 'y' in k:
            dicts_theo_mz_converted[k] = [x for x in v[::-1]]
        else:
            dicts_theo_mz_converted[k] = v
    return dicts_theo_mz_converted

def mh_calculator(wt_pep, ptm_mod, mh, charg, variant_aa_loc):
    print (wt_pep, ptm_mod, mh, charg, variant_aa_loc)
    mass_without_mod = mass.calculate_mass(sequence=wt_pep, ion_type='M', charge=charg)
    mass_with_mod = 0.0
    if ptm_mod != '':
        mod_val = 0.0
        for iter_mods in ptm_mod.split(';'):
            if re.search('\((.*?)\)', iter_mods).group(1).strip() in modification_list:
                mod_val = mod_val  + float(modification_list[re.search('\((.*?)\)', iter_mods).group(1).strip()]) / charg
        mass_with_mod = mass_without_mod + mod_val
    #print (mh, dicts_aa_mass[variant_aa_loc[-1]], tmp1, dicts_aa_mass[variant_aa_loc[0]], tmp2, mod_val)
    else:
        pass
    return mass_with_mod, mass_without_mod


def ppm_calc(ppm, theo_precursor_mass):
    ppm_calculate = (ppm / 1000000) * theo_precursor_mass
    return ppm_calculate
#print (ppm_calc(10, 554.55658))

raw_file_list = []

for i in os.listdir(infile_path):
    if '.mgf' in i:
        raw_file_list.append(i.split('.')[0])
        read_mgf(infile_path + "/" + i)

def generate_pos_neg_rt(rt_min):
    rt_pos = [x + rt_min + 1 for x in range(10)]
    rt_neg = [rt_min - x  - 1 for x in range(10)]
    return rt_pos, rt_neg

# Here we will be adding the - to the first row and last rows of b and y respectivelty
def draft_table(theoratical_spectra):
    dicts_theo = {}
    for k, v in theoratical_spectra.items():
        if 'y' in k:
            dicts_theo[k] = [float(str(x).split('.')[0] + "." + str(x).split('.')[1][0]  + str(x).split('.')[1][1]) for x in v]
            dicts_theo[k].insert(0, '-')
        else:
            dicts_theo[k] = [float(str(x).split('.')[0] + "." + str(x).split('.')[1][0]  + str(x).split('.')[1][1]) for x in v]
            dicts_theo[k].insert(len(dicts_theo[k]), '-')
    return dicts_theo

# Creating the table look like Proteome Discoverer
def create_table(peptide, charge, theo_spectra):
    #print (peptide)
    dict_len = (len(theo_spectra))
    trans_table = ""
    for i in range(len(peptide)):
        if dict_len == 6: # Charge 1
            trans_table = trans_table +  str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b1'][i]) + '\t' + \
                          peptide[i] + '\t' + str(theo_spectra['y1'][i]) + '\t' + str(theo_spectra['y1-H2O'][i]) + '\t' + str(theo_spectra['y1-NH3'][i]) + '\n'
        elif dict_len == 8: #charge 2
            trans_table = trans_table + str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b2'][i]) + '\t' + \
                          str(theo_spectra['b1'][i]) + '\t' + peptide[i] + '\t' + str(theo_spectra['y1'][i]) + '\t' + str(theo_spectra['y2'][i]) + \
                          '\t' + str(theo_spectra['y1-NH3'][i]) + '\t' + str(theo_spectra['y1-H2O'][i]) + '\n'
        elif dict_len == 10: #charge 3
            trans_table = trans_table + str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b3'][i]) + '\t' + \
                          str(theo_spectra['b2'][i]) + '\t' + str(theo_spectra['b1'][i]) + '\t' + peptide[i] + '\t' + str(theo_spectra['y1'][i]) + \
                          '\t' + str(theo_spectra['y2'][i]) + '\t' + str(theo_spectra['y3'][i]) +  '\t' + str(theo_spectra['y1-NH3'][i]) + '\t' + \
                          str(theo_spectra['y1-H2O'][i]) + '\n'            
        elif dict_len == 12: #charge 4
            trans_table = trans_table + str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b4'][i]) + '\t' + \
                          str(theo_spectra['b3'][i]) + '\t' + str(theo_spectra['b2'][i]) + '\t' + str(theo_spectra['b1'][i]) + '\t' + \
                          peptide[i] + '\t' + str(theo_spectra['y1'][i]) + '\t' + str(theo_spectra['y2'][i]) + '\t' + str(theo_spectra['y3'][i]) + \
                          '\t' + str(theo_spectra['y4'][i])+ '\t' + str(theo_spectra['y1-NH3'][i]) + '\t' + str(theo_spectra['y1-H2O'][i]) + '\n'
        elif dict_len == 14: #charge 5
            trans_table = trans_table + str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b5'][i]) + '\t' + \
                          str(theo_spectra['b4'][i]) + '\t' + str(theo_spectra['b3'][i]) + '\t' + str(theo_spectra['b2'][i]) + '\t' \
                          + str(theo_spectra['b1'][i]) + '\t'  + peptide[i] + '\t' + str(theo_spectra['y1'][i]) + '\t' + str(theo_spectra['y2'][i]) + \
                          '\t' + str(theo_spectra['y3'][i]) + '\t' + str(theo_spectra['y4'][i]) + '\t' + str(theo_spectra['y5'][i]) + \
                          '\t' + str(theo_spectra['y1-NH3'][i]) + '\t' + str(theo_spectra['y1-H2O'][i]) + '\n'
        elif dict_len == 16: #charge 6
            trans_table = trans_table + str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b6'][i]) + '\t' + \
                          str(theo_spectra['b5'][i]) + '\t' + str(theo_spectra['b4'][i]) + '\t' + str(theo_spectra['b3'][i]) + '\t' + str(theo_spectra['b2'][i]) +  \
                          '\t' + str(theo_spectra['b1'][i]) + '\t'  + peptide[i] + '\t' + str(theo_spectra['y1'][i]) + '\t' + str(theo_spectra['y2'][i]) + \
                          '\t' + str(theo_spectra['y3'][i]) + '\t' + str(theo_spectra['y4'][i]) + '\t' + str(theo_spectra['y5'][i])  + '\t' + \
                          str(theo_spectra['y6'][i]) + '\t' + str(theo_spectra['y1-NH3'][i]) + '\t' + str(theo_spectra['y1-H2O'][i]) + '\n'            
        elif dict_len == 18: #charge 7
            trans_table = trans_table + str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b7'][i]) + '\t' \
                          + str(theo_spectra['b6'][i]) + '\t' + str(theo_spectra['b5'][i]) + '\t' + str(theo_spectra['b4'][i]) + '\t' \
                          + str(theo_spectra['b3'][i]) + '\t' + str(theo_spectra['b2'][i]) + '\t' + str(theo_spectra['b1'][i]) + '\t'  + peptide[i] + '\t' \
                          + str(theo_spectra['y1'][i]) + '\t' + str(theo_spectra['y2'][i]) + '\t' + str(theo_spectra['y3'][i]) + '\t' + str(theo_spectra['y4'][i]) \
                          + '\t' + str(theo_spectra['y5'][i])  + '\t' + str(theo_spectra['y6'][i]) + '\t' + str(theo_spectra['y7'][i]) + '\t' + str(theo_spectra['y1-NH3'][i]) \
                          + '\t' + str(theo_spectra['y1-H2O'][i]) + '\n'
        elif dict_len == 20: #charge 8
            trans_table = trans_table + str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b8'][i]) + '\t' + \
                          str(theo_spectra['b7'][i]) + '\t' + str(theo_spectra['b6'][i]) + '\t' + str(theo_spectra['b5'][i]) + '\t' \
                          + str(theo_spectra['b4'][i]) + '\t' + str(theo_spectra['b3'][i]) + '\t' + str(theo_spectra['b2'][i]) + '\t' + \
                          str(theo_spectra['b1'][i]) + '\t'  + peptide[i] + '\t' + str(theo_spectra['y1'][i]) + '\t' + str(theo_spectra['y2'][i]) \
                          + '\t' + str(theo_spectra['y3'][i]) + '\t' + str(theo_spectra['y4'][i]) + '\t' + str(theo_spectra['y5'][i])  + '\t' + \
                          str(theo_spectra['y6'][i]) + '\t' + str(theo_spectra['y7'][i]) + '\t' + str(theo_spectra['y8'][i]) + '\t' + str(theo_spectra['y1-NH3'][i]) \
                          + '\t' + str(theo_spectra['y1-H2O'][i])+ '\n'
        elif dict_len == 22: #charge 9
            trans_table = trans_table + str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b9'][i]) + \
                          '\t' + str(theo_spectra['b8'][i]) + '\t' + str(theo_spectra['b7'][i]) + '\t' + str(theo_spectra['b6'][i]) + '\t' \
                          + str(theo_spectra['b5'][i]) + '\t' + str(theo_spectra['b4'][i]) + '\t' + str(theo_spectra['b3'][i]) + '\t' +\
                          str(theo_spectra['b2'][i]) + '\t' + str(theo_spectra['b1'][i]) +'\t'  + peptide[i] + '\t' + str(theo_spectra['y1'][i]) + \
                          '\t' + str(theo_spectra['y2'][i]) + '\t' + str(theo_spectra['y3'][i]) + '\t' + str(theo_spectra['y4'][i]) + '\t' + \
                          str(theo_spectra['y5'][i])  + '\t' + str(theo_spectra['y6'][i]) + '\t' + str(theo_spectra['y7'][i]) + '\t' + \
                          str(theo_spectra['y8'][i]) + '\t' + str(theo_spectra['y9'][i]) + '\t' + str(theo_spectra['y1-NH3'][i]) + '\t' + \
                          str(theo_spectra['y1-H2O'][i])+ '\n'
        elif dict_len == 24: #charge 10
            trans_table = trans_table + str(theo_spectra['b1-H2O'][i]) + '\t' + str(theo_spectra['b1-NH3'][i]) + '\t' + str(theo_spectra['b10'][i]) + \
                          '\t' + str(theo_spectra['b9'][i]) + '\t' + str(theo_spectra['b8'][i]) + '\t' + str(theo_spectra['b7'][i]) + '\t' \
                          + str(theo_spectra['b6'][i]) + '\t' + str(theo_spectra['b5'][i]) + '\t' + str(theo_spectra['b4'][i]) + '\t' + \
                          str(theo_spectra['b3'][i]) + '\t' + str(theo_spectra['b2'][i]) + '\t' + str(theo_spectra['b1'][i]) +'\t' \
                          + peptide[i] + '\t' + str(theo_spectra['y1'][i]) + '\t' + str(theo_spectra['y2'][i]) + '\t' + str(theo_spectra['y3'][i]) \
                          + '\t' + str(theo_spectra['y4'][i]) + '\t' + str(theo_spectra['y5'][i])  + '\t' + str(theo_spectra['y6'][i]) + \
                          '\t' + str(theo_spectra['y7'][i]) + '\t' + str(theo_spectra['y8'][i]) + '\t' + str(theo_spectra['y9'][i]) + '\t' \
                          + str(theo_spectra['y10'][i]) + '\t' + str(theo_spectra['y1-NH3'][i]) + '\t' + str(theo_spectra['y1-H2O'][i]) + '\n'

    return trans_table

def label_to_index(charge):
    ''' conver the labels to value for plot the table '''
    dict_label_idx = {}
    if charge == 1:
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b1':2, 'y1':4, 'y1-NH3':5, 'y1-H2O':6}
    elif charge == 2: 
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b2':2, 'b1':3, 'y1':5, 'y2':6, 'y1-NH3':7, 'y1-H2O':8}
    elif charge == 3: 
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b3':2, 'b2':3, 'b1':4, 'y1':6, 'y2':7, 'y3':8, 'y1-NH3':9, 'y1-H2O':10}
    elif charge == 4: 
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b4':2, 'b3':3, 'b2':4, 'b1':5, 'y1':7, 'y2':8, 'y3':9, 'y4':10, 'y1-NH3':11, 'y1-H2O':12}
    elif charge == 5: 
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b5':2, 'b4':3, 'b3':4, 'b2':5, 'b1':6, 'y1':8, 'y2':9, 'y3':10, 'y4':11, 'y5':12, 'y1-NH3':13, 'y1-H2O':14}
    elif charge == 6: 
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b6':2, 'b5':3, 'b4':4, 'b3':5, 'b2':6, 'b1':7, 'y1':9, 'y2':10, 'y3':11, 'y4':12, 'y5':13, 'y6':14, 'y1-NH3':15, \
                          'y1-H2O':16}
    elif charge == 7: 
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b7':2, 'b6':3, 'b5':4, 'b4':5, 'b3':6, 'b2':7, 'b1':8, 'y1':10, 'y2':11, 'y3':12, 'y4':13, 'y5':14, 'y6':15, \
                          'y7':16, 'y1-NH3':17, 'y1-H2O':18}
    elif charge == 8: 
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b8':2, 'b7':3, 'b6':4, 'b5':5, 'b4':6, 'b3':7, 'b2':8, 'b1':9, 'y1':11, 'y2':12, 'y3':13, 'y4':14, 'y5':15, 'y6':16, \
                          'y7':17, 'y8':18, 'y1-NH3':19, 'y1-H2O':20}
    elif charge == 9: 
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b9':2, 'b8':3, 'b7':4, 'b6':5, 'b5':6, 'b4':7, 'b3':8, 'b2':9, 'b1':10, 'y1':12, 'y2':13, 'y3':14, 'y4':15, 'y5':16, \
                          'y6':17, 'y7':18, 'y8':19, 'y9':20, 'y1-NH3':21, 'y1-H2O':22}
    elif charge == 10: 
        dict_label_idx = {'b1-H2O':0, 'b1-NH3':1, 'b10':2, 'b9':3, 'b8':4, 'b7':5, 'b6':6, 'b5':7, 'b4':8, 'b3':9, 'b2':10, 'b1':11, 'y1':13, 'y2':14, 'y3':15, 'y4':16, \
                          'y5':17, 'y6':18, 'y7':19, 'y8':20, 'y9':21, 'y10':22, 'y1-NH3':23, 'y1-H2O':24}
        
    return dict_label_idx
        
def table_plot(lst, peptide, charge, theo_spectra, raw_file_name_for_store):
    fig, ax = plt.subplots()
    ax.axis('off')
    ax.axis('tight')
    peps = list(peptide)
    collabel = generate_label(charge)
    collabel.insert(charge, "seq")
    collabel = (list((x) for x in collabel))
    collabel.insert(0, "b1-H2O")
    collabel.insert(1, "b1-NH3")
    collabel.insert(len(collabel) + 1 , "y1-NH3")
    collabel.insert(len(collabel) + 2 , "y1-H2O")
    ###print (collabel)
##    for x, y in theo_spectra.items():
##        print (x, len(y), len(peps))
    if charge == 1:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], peps , theo_spectra['b1'], \
                                                 peps, theo_spectra['y1'], theo_spectra['y1-NH3'], theo_spectra['y1-H2O'])) \
                                                ,colLabels=collabel,loc='center')
    elif charge == 2:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], theo_spectra['b2'],theo_spectra['b1'], peps, \
                                                 theo_spectra['y1'], theo_spectra['y2'], theo_spectra['y1-NH3'], \
                                                 theo_spectra['y1-H2O'])), colLabels=collabel,loc='center')
    elif charge == 3:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], theo_spectra['b3'],theo_spectra['b2'], theo_spectra['b1'], \
                                                 peps, theo_spectra['y1'], theo_spectra['y2'], theo_spectra['y3'], theo_spectra['y1-NH3'], \
                                                 theo_spectra['y1-H2O'])),colLabels=collabel,loc='center')
    elif charge == 4:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], theo_spectra['b4'],theo_spectra['b3'], \
                                                 theo_spectra['b2'],theo_spectra['b1'], peps, theo_spectra['y1'], theo_spectra['y2'], \
                                                 theo_spectra['y3'], theo_spectra['y4'], theo_spectra['y1-NH3'], \
                                                 theo_spectra['y1-H2O'])),colLabels=collabel,loc='center')
    elif charge == 5:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], theo_spectra['b5'],theo_spectra['b4'], theo_spectra['b3'], \
                                                 theo_spectra['b2'],theo_spectra['b1'], peps, theo_spectra['y1'], theo_spectra['y2'], \
                                                 theo_spectra['y3'], theo_spectra['y4'], theo_spectra['y5'], theo_spectra['y1-NH3'], \
                                                 theo_spectra['y1-H2O'])),colLabels=collabel,loc='center')
    elif charge == 6:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], theo_spectra['b6'],theo_spectra['b5'], \
                                                 theo_spectra['b4'],theo_spectra['b3'],theo_spectra['b2'],theo_spectra['b1'], peps, theo_spectra['y1'], \
                                                 theo_spectra['y2'], theo_spectra['y3'], theo_spectra['y4'], theo_spectra['y5'], \
                                                 theo_spectra['y6'], theo_spectra['y1-NH3'], theo_spectra['y1-H2O'])),colLabels=collabel,loc='center')
    elif charge == 7:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], theo_spectra['b7'],theo_spectra['b6'], theo_spectra['b5'], \
                                                 theo_spectra['b4'],theo_spectra['b3'],theo_spectra['b2'], theo_spectra['b1'], \
                                                 peps, theo_spectra['y1'],theo_spectra['y2'], theo_spectra['y3'], theo_spectra['y4'], theo_spectra['y5'], \
                                                 theo_spectra['y6'], theo_spectra['y7'], theo_spectra['y1-NH3'], \
                                                 theo_spectra['y1-H2O'])),colLabels=collabel,loc='center')
    elif charge == 8:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], theo_spectra['b8'],theo_spectra['b7'], \
                                                 theo_spectra['b6'],theo_spectra['b5'],theo_spectra['b4'],theo_spectra['b3'], \
                                                 theo_spectra['b2'], theo_spectra['b1'], peps, theo_spectra['y1'],theo_spectra['y2'], theo_spectra['y3'], \
                                                 theo_spectra['y4'], theo_spectra['y5'], theo_spectra['y6'], theo_spectra['y7'], theo_spectra['y8'], \
                                                 theo_spectra['y1-NH3'], theo_spectra['y1-H2O'])),colLabels=collabel,loc='center')
    elif charge == 9:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], theo_spectra['b9'],theo_spectra['b8'], theo_spectra['b7'], \
                                                 theo_spectra['b6'], theo_spectra['b5'],theo_spectra['b4'], \
                                                 theo_spectra['b3'], theo_spectra['b2'], theo_spectra['b1'], peps, theo_spectra['y1'], \
                                                 theo_spectra['y2'], theo_spectra['y3'], theo_spectra['y4'], \
                                                 theo_spectra['y5'], theo_spectra['y6'], theo_spectra['y7'], theo_spectra['y8'], theo_spectra['y9'], \
                                                 theo_spectra['y1-NH3'], theo_spectra['y1-H2O'])),colLabels=collabel,loc='center')
    elif charge == 10:
        tbl = ax.table(cellText=np.column_stack((theo_spectra['b1-H2O'], theo_spectra['b1-NH3'], theo_spectra['b10'],theo_spectra['b9'], \
                                                 theo_spectra['b8'],theo_spectra['b7'],theo_spectra['b6'],theo_spectra['b5'], \
                                                 theo_spectra['b4'], theo_spectra['b3'] , theo_spectra['b2'], theo_spectra['b1'], peps, \
                                                 theo_spectra['y1'],theo_spectra['y2'], theo_spectra['y3'], theo_spectra['y4'], theo_spectra['y5'], \
                                                 theo_spectra['y6'], theo_spectra['y7'], theo_spectra['y8'],theo_spectra['y9'], \
                                                 theo_spectra['y10'], theo_spectra['y1-NH3'], theo_spectra['y1-H2O'])),colLabels=collabel,loc='center')

##    tbl._cells[9, 1].set_facecolor("#33b5e5")
    ###print (lst, "All list")
    
    for color_mz in range(len(lst)):
        split_color_mz = lst[color_mz].split('_')
        if charge == 1:
            if  int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1 , int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
        elif charge == 2:
            if  int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1 , int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
        elif charge == 3:
            if  int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1 , int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
        elif charge == 4:
            if int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1, int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
        elif charge == 5:
            if  int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1 , int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
        elif charge == 6:
            if  int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1 , int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
        elif charge == 7:
            if  int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1 , int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
        elif charge == 8:
            if  int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1 , int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
        elif charge == 9:
            if  int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1 , int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
        elif charge == 10:
            if  int(split_color_mz[1]) <= charge + 1:
                tbl._cells[(int(split_color_mz[0]) + 1 , int(split_color_mz[1]))].set_facecolor("#f06292")
            else:
                tbl._cells[(int(split_color_mz[0]) + 2, int(split_color_mz[1]))].set_facecolor("#33b5e5")
    plt.savefig(peptide + "_" + raw_file_name_for_store +".png", bbox_inches='tight', dpi=600)
    plt.close()
    
def match_precursor_mass(peptide, charge, modifications, exp_mz, raw_files_for_store):
    dicts_matched_mz = {}
    dicts_matched_mz_idx = {}
    dicts_matched_mz_idx_plot = []
    dicts_all_mz = {}
    theo_spec = get_theo_spectra_with_mod(peptide, charge, modifications) #Peptide, charge
    draft_tabl = draft_table(theo_spec) # adding - to the table
    create_tabl = create_table(peptide, charge, draft_tabl) # A table like PD
    get_table_idx = label_to_index(charge) # get the index of the table
    for k in exp_mz:#dicts_mz[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass]:
            split_k = k.split(' ')
            for x, y in theo_spec.items():
                #print (x, y, len(split_i[8]))
                score_theo = generate_score_theo(y, len(peptide))
                #print (y, peptide)
                if x not in dicts_all_mz:
                    dicts_all_mz[x] = score_theo #Score for theo mz
                #Calculate santosh formula here for theoratical m/z
                for z in range(len(y)):
                    if abs(float(split_k[0]) - float(y[z])) <= 0.05:
                        if x not in dicts_matched_mz:
                            #print (str(z) + "_" + str(get_table_idx[x]), y[z])
                            dicts_matched_mz[x] = [y[z]]
                            dicts_matched_mz_idx[x] = [z]
                            #print (str(z) + "_" + str(get_table_idx[x]), x, y[z])
                            dicts_matched_mz_idx_plot.append(str(z) + "_" + str(get_table_idx[x]))
                            #print (str(z) + "_" + str(get_table_idx[x]), x, y[z])
                        else:
                            #print (str(z) + "_" + str(get_table_idx[x]), x, y[z])
                            dicts_matched_mz_idx_plot.append(str(z) + "_" + str(get_table_idx[x]))
                            dicts_matched_mz[x].append(y[z])
                            dicts_matched_mz_idx[x].append(z)
    cnt_maches = 0
    table_plot(dicts_matched_mz_idx_plot, peptide, int(charge), draft_tabl, raw_files_for_store)

    for o, l in dicts_matched_mz.items():
        cnt_maches+=len(set(l))
    return cnt_maches

##    #print (dicts_matched_mz_idx_plot)
##    dicts_final_score = {}
##    table_plot(dicts_matched_mz_idx_plot, peptide, int(charge), draft_tabl, raw_files_for_store)
##    for a, b in dicts_matched_mz.items():
##        ###print (a, b, dicts_matched_mz_idx[a],peptide)
##        matched_score = generate_socre_matched(b, dicts_matched_mz_idx[a], len(peptide))
##        dicts_final_score[a] = matched_score / dicts_all_mz[a]
##    for m, n in dicts_final_score.items():
##        ##print (m, (n*100))
##        pass
##    print("*********************************************************************")

def decoy_search(rev_peptide, charge, modification, exp_mz, raw_files_to_store):
    decoy_dicts_matched_mz = {}
    decoy_dicts_matched_mz_idx = {}
    decoy_dicts_matched_mz_idx_plot = []
    decoy_dicts_all_mz = {}
    rev_pep_and_mod = reverse_peptide_for_decoy(rev_peptide, modification)
    rev_theo_spec = get_theo_spectra_with_mod(rev_pep_and_mod[0], charge, rev_pep_and_mod[1]) #Peptide, charge
    rev_draft_tabl = draft_table(rev_theo_spec) # adding - to the table
    rev_create_tabl = create_table(rev_pep_and_mod[0], charge, rev_draft_tabl) # A table like PD
    get_table_idx = label_to_index(charge)
    for k in exp_mz:#dicts_mz[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass]:
            split_k = k.split(' ')
            for x, y in rev_theo_spec.items():
                #print (x, y, len(split_i[8]))
                score_theo = generate_score_theo(y, len(rev_peptide))
                #print (y, peptide)
                if x not in decoy_dicts_all_mz:
                    decoy_dicts_all_mz[x] = score_theo #Score for theo mz
                #Calculate santosh formula here for theoratical m/z
                for z in range(len(y)):
                    if abs(float(split_k[0]) - float(y[z])) <= 0.05:
                        if x not in decoy_dicts_matched_mz:
                            #print (str(z) + "_" + str(get_table_idx[x]), y[z])
                            decoy_dicts_matched_mz[x] = [y[z]]
                            decoy_dicts_matched_mz_idx[x] = [z]
                            #print (str(z) + "_" + str(get_table_idx[x]), x, y[z])
                            decoy_dicts_matched_mz_idx_plot.append(str(z) + "_" + str(get_table_idx[x]))
                            #print (str(z) + "_" + str(get_table_idx[x]), x, y[z])
                        else:
                            #print (str(z) + "_" + str(get_table_idx[x]), x, y[z])
                            decoy_dicts_matched_mz_idx_plot.append(str(z) + "_" + str(get_table_idx[x]))
                            decoy_dicts_matched_mz[x].append(y[z])
                            decoy_dicts_matched_mz_idx[x].append(z)

    cnt_matches = 0
    for o, l in decoy_dicts_matched_mz.items():
        cnt_matches+=len(set(l))
    return cnt_matches
    

output_file.write('id' + '\t' + 'TD' + '\t' + 'ScanNr' + '\t' + 'numberOfMatchingPeaks' + '\t' + 'charge' + '\t' + 'Exp.Mz' + '\t' + 'RTinSec' + '\t' +  'Theo.Mz' + '\t' + 'Delta_Mass_Error(PPM)' + '\t' + 'length' + '\t' + 'peptide' + '\t' + 'protein' + '\n')                     
with open(input_fl) as f:
    for i in islice(f, 1, None):
        split_i = i.split('\t')
        a = mh_calculator(split_i[8], split_i[1], float(split_i[3]), int(split_i[2]), split_i[9].rstrip())
        raw_file_lists = raw_file_list# #Without file extension example:- Q0145
        for calc_precursor_mass in a: #mh with mod and without mod
            if calc_precursor_mass != 0.0:
                precursor_mass_calc = float("%.5f" % calc_precursor_mass)
                rt_calc = generate_pos_neg_rt(int(float(split_i[5]) * 60))  # add the for loop after this line for chekc in previous and next raw files
                ppm_window = ppm_calc(10, precursor_mass_calc)
                for raw_file in raw_file_lists:
                    for iter_rt in range(len(rt_calc[0])): # rt + increment = positive values
                        if rt_calc[0][iter_rt] in dicts_pepmass[raw_file]:
                            for iter_precursor_mass in range(len(dicts_pepmass[raw_file][rt_calc[0][iter_rt]])):
                                #print (precursor_mass_calc, float(dicts_pepmass[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass].split(' ')[0]))
                                if abs(precursor_mass_calc - float(dicts_pepmass[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass].split(' ')[0])) <= ppm_window:
                                    #print (dicts_pepmass[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass], iter_precursor_mass, precursor_mass_calc, "Pos Yes", split_i[8], calc_precursor_mass, precursor_mass_calc, split_i[2] + "Charge" + raw_file)
                                    aa = match_precursor_mass(split_i[8], int(split_i[2]), split_i[1], dicts_mz[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass], raw_file)
                                    id_col = raw_file + ":" + str(dicts_scan[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass]) + ":" + split_i[8] + ":" + dicts_charge[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass] + ":" + split_i[1]
                                    output_file.write(id_col + '\t' +  "1" + '\t' + str(dicts_scan[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass]) + '\t' + str(aa) + '\t' + dicts_charge[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass] + '\t' + dicts_pepmass[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass].split(' ')[0] + '\t' + dicts_rt[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass] + '\t' + str(precursor_mass_calc) + '\t' +  str(float(dicts_pepmass[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass].split(' ')[0]) - float(precursor_mass_calc)) + '\t' + str(len(split_i[8])) + '\t' + split_i[8] + '\t' + split_i[11].rstrip() + '\n')
                                    
                                    dicts_mod_pos = {}
                                    rev_aa = decoy_search(split_i[8], int(split_i[2]), split_i[1], dicts_mz[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass], raw_file)
                                    pep_and_mods = reverse_peptide_for_decoy(split_i[8], split_i[1])
                                    id_col = raw_file + ":" + str(dicts_scan[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass]) + ":" + "decoy_" + pep_and_mods[0] + ":" + dicts_charge[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass] + ":" + pep_and_mods[1]
                                    output_file.write(id_col + '\t' + "-1" + '\t' + str(dicts_scan[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass]) + '\t' + str(rev_aa) + '\t' + dicts_charge[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass] + '\t' + dicts_pepmass[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass].split(' ')[0] + '\t' + dicts_rt[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass] + '\t' + str(precursor_mass_calc) + '\t' + str(float(dicts_pepmass[raw_file][rt_calc[0][iter_rt]][iter_precursor_mass].split(' ')[0]) - float(precursor_mass_calc)) + '\t' + str(len(pep_and_mods[0])) + '\t' + pep_and_mods[0] + '\t' + "decoy_" + split_i[11].rstrip() + '\n')
                                    dicts_mod_pos = {}

                    for iter_rt in range(len(rt_calc[1])): # rt - increment = negative values
                        if rt_calc[1][iter_rt] in dicts_pepmass[raw_file]:
                            #print (rt_calc[1][iter_rt], "!!!!!!!!!!!!!!! Negative", split_i[5], iter_rt, split_i[8])
                            for iter_precursor_mass in range(len(dicts_pepmass[raw_file][rt_calc[1][iter_rt]])):
                                if abs(precursor_mass_calc - float(dicts_pepmass[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass].split(' ')[0])) <= ppm_window:
                                    #print (dicts_pepmass[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass], iter_precursor_mass, precursor_mass_calc, "Neg Yes", split_i[8], calc_precursor_mass, precursor_mass_calc, split_i[2] + "Charge" + raw_file)
                                    # variant_peptide wildtype_peptide, m/z:intensity, charge, scan_number, RT, raw_file
                                    aa = match_precursor_mass(split_i[8], int(split_i[2]), split_i[1], dicts_mz[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass], raw_file)
                                    id_col = raw_file + ":" + str(dicts_scan[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass]) + ":" + split_i[8] + ":" + dicts_charge[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass] + ":" + split_i[1]
                                    output_file.write(id_col + '\t' +  "1" + '\t' + str(dicts_scan[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass]) + "\t" + str(aa) + '\t' + dicts_charge[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass] + '\t' + dicts_pepmass[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass].split(' ')[0] + '\t' + dicts_rt[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass] + '\t' + str(precursor_mass_calc) + '\t' + str(float(dicts_pepmass[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass].split(' ')[0]) - float(precursor_mass_calc)) + '\t'  + str(len(split_i[8])) + '\t' + split_i[8] + '\t' + split_i[11].rstrip() + '\n')
                                    
                                    dicts_mod_pos = {}
                                    rev_aa = decoy_search(split_i[8], int(split_i[2]), split_i[1], dicts_mz[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass], raw_file)
                                    pep_and_mods = reverse_peptide_for_decoy(split_i[0], split_i[1])
                                    id_col = raw_file + ":" + str(dicts_scan[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass]) + ":" + "decoy_" + pep_and_mods[0] + ":" + dicts_charge[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass] + ":" + pep_and_mods[1]
                                    output_file.write(id_col + '\t' + "-1" + '\t' + str(dicts_scan[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass]) + "\t" + str(rev_aa) + '\t' + dicts_charge[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass] + '\t' + dicts_pepmass[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass].split(' ')[0] + '\t' + dicts_rt[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass] + '\t' + str(precursor_mass_calc) + '\t' + str(float(dicts_pepmass[raw_file][rt_calc[1][iter_rt]][iter_precursor_mass].split(' ')[0]) - float(precursor_mass_calc)) + '\t' + str(len(pep_and_mods[0])) + '\t' + pep_and_mods[0] + '\t' + "decoy_" + split_i[11].rstrip() + '\n')
                                    dicts_mod_pos = {}
output_file.close()
                                
