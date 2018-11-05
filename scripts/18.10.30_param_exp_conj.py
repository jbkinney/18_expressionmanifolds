import os
import pylab as py
import pandas as pd
import matplotlib
import numpy as np
import scipy.optimize as opt
from matplotlib import pyplot as plt
import sys
from random import randint
from random import uniform
import pickle

py.close('all')

#specify folder for data files
datapath='../data/'

#plate_dat is a pandas panel, indexed by plate name, that stores well, construct name, 
#miller measurement, OD600, reaction rate, cAMP concentration, volume, glycerol location, 
#and whether or not the seq is good
plate_dat=pd.read_pickle(datapath+'plate_panel.pkl')

#glycerol stocks contains info about each construct
cons_names=pd.read_excel(datapath+'glycerol_stocks.xlsx', index_col=0)

#seq spacing has edited sequence and spacing info about a limited number of constructs
seq_sum=pd.read_excel(datapath+'seq_spacing.xlsx', index_col=0)

#library clusters has dictionaries of all constructs belonging to each library
exec(open(datapath+'library_clusters.py').read(), globals())

#this script is only interested in the occlusion library
curves={'c60':'c60r18.10', 'c61':'c61r18', 'c62':'c62r18.10', 'c63':'c63r18.10', 'c64':'c64r18.10', 
        'c65':'c65r18.10', 'c66':'c66r18.10', 'c71':'c71r18', 'c72':'c72r18.10', 'c76':'c76r18.10', 
        'c81':'c81r18.10', 'c82':'c82r18.10'}

raw_vals={}

for plate in plate_dat:
    for val in range(len(plate_dat[plate])):
        nam=plate_dat[plate]["glycerol loc"][val]
        
        if nam=='RDM' or nam=='blank':
            continue
            
        if plate_dat[plate]['OK seq?'][val]=='FALSE' or float(plate_dat[plate]['OD600'][val])>0.5:
            continue
        
        if nam not in raw_vals:
            raw_vals[nam]={'0.0':[], '2.5':[], '5.0':[], '10.0':[],
                           '25.0':[], '50.0':[], '125.0':[], '250.0':[]}
        
        if plate_dat[plate]['cAMP'][val] not in raw_vals[nam]:
                raw_vals[nam][plate_dat[plate]['cAMP'][val]]=[]
        
        if np.isnan(float(plate_dat[plate]['Miller meas'][val])): continue
        if float(plate_dat[plate]['Miller meas'][val])<0: continue
        
        raw_vals[nam][plate_dat[plate]['cAMP'][val]].append(float(plate_dat[plate]['Miller meas'][val]))
        
mins=2  

corrected_vals={}

for cons in raw_vals:
    if not cons_names['valid_clone'][cons]: continue
    corrected_vals[cons]={}
    for conc in raw_vals[cons]:
        if conc=='0.0':
            corrected_vals[cons][conc]=[0.855*x for x in raw_vals[cons][conc]]
        else:
            corrected_vals[cons][conc]=raw_vals[cons][conc]
            
corrected_medians={}

for sample in corrected_vals:
    if not cons_names['valid_clone'][sample]: continue
    corrected_medians[sample]={}
    
    for conc in corrected_vals[sample]:
        bas_temp=np.nan
        raw_bas_screened=[]
        
        for val in corrected_vals[sample][conc]:
            if val>0:
                raw_bas_screened.append(val)
                
        if len(raw_bas_screened)>mins:
            bas_temp=np.median(raw_bas_screened)
            
        corrected_medians[sample][conc]=bas_temp


all_vects={}
med_vals={}
for count in range(len(seq_sum['location'])):
    cons=seq_sum['location'][count]
    med_vals[cons]={'bas':corrected_medians[cons]['0.0'], 'ind':corrected_medians[cons]['250.0']}

for spac in curves:
    lib=curves[spac]
    all_vects[lib]={'lt+':[], 'lt-':[], 'locs':[]}
    for count in range(len(seq_sum['location'])):
        cons=seq_sum['location'][count]
        if cons not in library_groups[lib]['all']: continue
        if cons in library_groups[lib]['outliers']: continue
        if np.isnan(med_vals[cons]['bas']): continue
        if np.isnan(med_vals[cons]['ind']): continue
            
        all_vects[lib]['lt+'].append(np.log(med_vals[cons]['ind']))
        all_vects[lib]['lt-'].append(np.log(med_vals[cons]['bas']))
        all_vects[lib]['locs'].append(cons)

def comp_ltaui(lrn,lthetas): 
    """
    Calculates and returns the log of the induced transcription.
    
    Parameters
    ----------
    lrn : float
        The log of the P for this promoter.
        
    lthetas : float array
        lthetas[0] is the log of the saturated transcription
        lthetas[1] is the log of the background transcription
        lthetas[2] is the log of the eta (eta = (1+Falpha)/(1+F))
    
    Returns
    -------
    float
        The log of the calculated induced transcription,
        induced transc = ts*(eta*P/(1+eta*P))
    """

    P=np.exp(lrn)
    ts=np.exp(lthetas[0])
    tbg=np.exp(lthetas[1])
    et=np.exp(lthetas[2])

    
    #partition fn
    Z=1+P*et
    
    #total transc
    ti=(ts*P*et/Z)+tbg
    
    return np.log(ti)
    
#calcs log of tau basal from same givens as above
def comp_ltaub(lrn,lthetas): 
    """
    Calculates and returns the log of the basal transcription.
    
    Parameters
    ----------
    lrn : float
        The log of the P for this promoter.
        
    lthetas : float array
        lthetas[0] is the log of the saturated transcription
        lthetas[1] is the log of the background transcription
        lthetas[2] is the log of the eta (eta = (1+Falpha)/(1+F)) (not used)
    
    Returns
    -------
    float
        The log of the calculated basal transcription,
        basal transc = ts*(P/(1+P))
    """

    P=np.exp(lrn)
    ts=np.exp(lthetas[0])
    tbg=np.exp(lthetas[1])
    et=np.exp(lthetas[2])
    
    Z=1+P
    
    #total transc
    tb=(ts*P/Z)+tbg
    
    return np.log(tb)
    
#objective function, given the measured values for log basal transc and log induced transc,
#and a vector for the log of all variables fitted
#returns the squared difference in calculated log induced transc from each measured log induc transc we've got
#plus calculated log basal transc from each measured log basal transc
def comp_error(x, lbasal, linduced, weights):
    
    err=0
    lib_ct=len(weights)
    counter=0
    
    for n, lib in enumerate(weights):
        tsat=x[-(2*lib_ct+1)]
        tbackg=x[-(2*lib_ct-n)]
        eta=x[-(lib_ct-n)]
        
        errbas=[]
        errind=[]
        
        lt_minus=lbasal[n]
        lt_plus=linduced[n]
        
        mocker=[]
        
        for count in range(len(lib)):
            if lib[count]==0:
                mocker.append(0)
            else:
                mocker.append(x[counter])
                counter+=1
                
        for y in range(len(mocker)):
            
            errb_temp=lt_minus[y]-comp_ltaub(mocker[y], [tsat, tbackg, eta])
            erri_temp=lt_plus[y]-comp_ltaui(mocker[y], [tsat, tbackg, eta])
            
            errbas.append(errb_temp)
            errind.append(erri_temp)
            
        err_temp=np.transpose(lib)@(np.square(errbas)+np.square(errind))
        
        err+=err_temp
        
    return err

tsat_set=[np.log(1),np.log(10),np.log(100)]
tbg_set=[np.log(0.0001),np.log(0.001),np.log(0.01)]
alphp_set=[np.log(50),np.log(500),np.log(5000)]
p_set=[np.log(0.01),np.log(0.1)]
    
fit_set=[]
for lib in curves:
    if 'c-' in lib: continue
    if 'oc' in lib: continue
    if 'r17' in lib: continue
    
    fit_set.append(curves[lib])

for cons in fit_set:
    temp=[]
    for x in all_vects[cons]['lt-']:
        rtemp=p_set[0]
        temp.append(rtemp)
        
    all_vects[cons]['initial']=temp
    all_vects[cons]['full_init']=temp+[tsat_set[0],tbg_set[0],alphp_set[0]]
    
conjoined={'lt+':[], 'lt-':[], 'locs':[], 'initial':[], 'curve_list':[]}

for lib in fit_set:
    if 'c-' in lib: continue
    conjoined['curve_list'].append(lib)
    conjoined['lt+'].append(all_vects[lib]['lt+'])
    conjoined['lt-'].append(all_vects[lib]['lt-'])
    
    for category in ['locs', 'initial']:
        for cons in all_vects[lib][category]:
            conjoined[category].append(cons)
            
all_vects['conj']=conjoined
    
fits={}
mfts={}

weight_test=[]

for lib in all_vects['conj']['curve_list']:
    weight_test.append(np.ones_like(all_vects[lib]['lt-']))

def objective_func_creator_conj(lib):
    bas=all_vects[lib]['lt-']
    ind=all_vects[lib]['lt+']
    lammy=lambda m: comp_error(m,bas,ind,weight_test)
    return lammy

all_vects['conj']['objective_func']=objective_func_creator_conj('conj')

inits_set={
    'aa':[tsat_set[0],tbg_set[0],alphp_set[0]],
    'ab':[tsat_set[0],tbg_set[0],alphp_set[1]],
    'ac':[tsat_set[0],tbg_set[0],alphp_set[2]],
    'ad':[tsat_set[0],tbg_set[1],alphp_set[0]],
    'ae':[tsat_set[0],tbg_set[1],alphp_set[1]],
    'af':[tsat_set[0],tbg_set[1],alphp_set[2]],
    'ag':[tsat_set[0],tbg_set[2],alphp_set[0]],
    'ah':[tsat_set[0],tbg_set[2],alphp_set[1]],
    'ai':[tsat_set[0],tbg_set[2],alphp_set[2]],
    'ba':[tsat_set[1],tbg_set[0],alphp_set[0]],
    'bb':[tsat_set[1],tbg_set[0],alphp_set[1]],
    'bc':[tsat_set[1],tbg_set[0],alphp_set[2]],
    'bd':[tsat_set[1],tbg_set[1],alphp_set[0]],
    'be':[tsat_set[1],tbg_set[1],alphp_set[1]],
    'bf':[tsat_set[1],tbg_set[1],alphp_set[2]],
    'bg':[tsat_set[1],tbg_set[2],alphp_set[0]],
    'bh':[tsat_set[1],tbg_set[2],alphp_set[1]],
    'bi':[tsat_set[1],tbg_set[2],alphp_set[2]],
    'ca':[tsat_set[2],tbg_set[0],alphp_set[0]],
    'cb':[tsat_set[2],tbg_set[0],alphp_set[1]],
    'cc':[tsat_set[2],tbg_set[0],alphp_set[2]],
    'cd':[tsat_set[2],tbg_set[1],alphp_set[0]],
    'ce':[tsat_set[2],tbg_set[1],alphp_set[1]],
    'cf':[tsat_set[2],tbg_set[1],alphp_set[2]],
    'cg':[tsat_set[2],tbg_set[2],alphp_set[0]],
    'ch':[tsat_set[2],tbg_set[2],alphp_set[1]],
    'ci':[tsat_set[2],tbg_set[2],alphp_set[2]],
}

for inits in inits_set:
    print(inits)
    temp=[]
    for x in all_vects['conj']['lt-']:
        for y in x:
            temp.append(p_set[0])
        
    all_vects['conj']['initial']=temp
    
    init_conj=[inits_set[inits][0]]
    for cons in conjoined['curve_list']:
        init_conj.append(inits_set[inits][1])
    for cons in conjoined['curve_list']:
        init_conj.append(inits_set[inits][2])
    
    all_vects['conj']['full_init']=temp+init_conj
    
    fits[inits]=opt.minimize(all_vects['conj']['objective_func'], all_vects['conj']['full_init'], method='L-BFGS-B', options={'maxfun':1000000})

for inits in inits_set:
    print(inits)
    temp=[]
    for x in all_vects['conj']['lt-']:
        for y in x:
            temp.append(p_set[1])
        
    all_vects['conj']['initial']=temp
    
    init_conj=[inits_set[inits][0]]
    for cons in conjoined['curve_list']:
        init_conj.append(inits_set[inits][1])
    for cons in conjoined['curve_list']:
        init_conj.append(inits_set[inits][2])
    
    all_vects['conj']['full_init']=temp+init_conj
    
    fits['d'+inits]=opt.minimize(all_vects['conj']['objective_func'], all_vects['conj']['full_init'], method='L-BFGS-B', options={'maxfun':1000000})

with open('../intermediates/conj_init_param_exp.pkl', 'wb') as handle:
    pickle.dump(fits, handle, protocol=pickle.HIGHEST_PROTOCOL)