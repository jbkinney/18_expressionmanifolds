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
curves={'occlusion':'occlusion'}

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
        
#returns calculated log induced transc given logs of:
#rn=RNAP binding weight to this sequence, tm=max WT lac transc, tbg=background transc, 
#ec=CRP binding weight
def comp_ltaui_oc(lrn,lthetas): 
    rn=np.exp(lrn)
    tm=np.exp(lthetas[0])
    tbg=np.exp(lthetas[1])
    ec=np.exp(lthetas[2])

    
    #partition fn
    Z=1+rn+ec
    
    #total transc
    ti=(tm*rn/Z)+tbg
    
    return np.log(ti)

#calcs log of tau basal from same givens as above
def comp_ltaub_oc(lrn,lthetas): 
    rn=np.exp(lrn)
    tm=np.exp(lthetas[0])
    tbg=np.exp(lthetas[1])
    
    Z=1+rn
    
    #total transc
    tb=(tm*rn/Z)+tbg
    
    return np.log(tb)
    
def comp_error_oc(x, lbasal, linduced, weights):
    err=[]
    mocker=[]
    counter=0
    
    for count in range(len(weights)):
        if weights[count]==0:
            mocker.append(0)
        else:
            mocker.append(x[counter])
            counter+=1
    
    for y in range(len(lbasal)):
        
        if mocker[y]==0: continue
        
        #subtract the calculated log tau basal and log tau induced from measured
        errb=lbasal[y]-comp_ltaub_oc(mocker[y],x[-3:]) 
        erri=linduced[y]-comp_ltaui_oc(mocker[y],x[-3:])
        
        #add squares of those differences
        err.append(((errb**2)+(erri**2))*weights[y])
        
    return sum(err)
    
all_vects['occlusion']['objective_func']=lambda m: comp_error_oc(m,all_vects['occlusion']['lt-'],all_vects['occlusion']['lt+'],np.ones_like(all_vects['occlusion']['lt-']))

tsat_set=[np.log(1),np.log(10),np.log(100)]
tbg_set=[np.log(0.0001),np.log(0.001),np.log(0.01)]
f_set=[np.log(1),np.log(10),np.log(100)]
p_set=[np.log(0.01),np.log(0.1)]

tm=tsat_set[0]
tbg=tbg_set[0]
ec=f_set[0]

init_oc=[tm,tbg,ec]
fit_oc={}
mfts_oc={}

temp=[]
for x in all_vects['occlusion']['lt-']:
    temp.append(p_set[0])
        
all_vects['occlusion']['initial']=temp
all_vects['occlusion']['full_init']=temp+init_oc

inits_set={
    'aa':[tsat_set[0],tbg_set[0],f_set[0]],
    'ab':[tsat_set[0],tbg_set[0],f_set[1]],
    'ac':[tsat_set[0],tbg_set[0],f_set[2]],
    'ad':[tsat_set[0],tbg_set[1],f_set[0]],
    'ae':[tsat_set[0],tbg_set[1],f_set[1]],
    'af':[tsat_set[0],tbg_set[1],f_set[2]],
    'ag':[tsat_set[0],tbg_set[2],f_set[0]],
    'ah':[tsat_set[0],tbg_set[2],f_set[1]],
    'ai':[tsat_set[0],tbg_set[2],f_set[2]],
    'ba':[tsat_set[1],tbg_set[0],f_set[0]],
    'bb':[tsat_set[1],tbg_set[0],f_set[1]],
    'bc':[tsat_set[1],tbg_set[0],f_set[2]],
    'bd':[tsat_set[1],tbg_set[1],f_set[0]],
    'be':[tsat_set[1],tbg_set[1],f_set[1]],
    'bf':[tsat_set[1],tbg_set[1],f_set[2]],
    'bg':[tsat_set[1],tbg_set[2],f_set[0]],
    'bh':[tsat_set[1],tbg_set[2],f_set[1]],
    'bi':[tsat_set[1],tbg_set[2],f_set[2]],
    'ca':[tsat_set[2],tbg_set[0],f_set[0]],
    'cb':[tsat_set[2],tbg_set[0],f_set[1]],
    'cc':[tsat_set[2],tbg_set[0],f_set[2]],
    'cd':[tsat_set[2],tbg_set[1],f_set[0]],
    'ce':[tsat_set[2],tbg_set[1],f_set[1]],
    'cf':[tsat_set[2],tbg_set[1],f_set[2]],
    'cg':[tsat_set[2],tbg_set[2],f_set[0]],
    'ch':[tsat_set[2],tbg_set[2],f_set[1]],
    'ci':[tsat_set[2],tbg_set[2],f_set[2]],
}

for inits in inits_set:
    print(inits)
    temp=[]
    for x in all_vects['occlusion']['lt-']:
        temp.append(p_set[0])
        
    all_vects['occlusion']['initial']=temp
    all_vects['occlusion']['full_init']=temp+inits_set[inits]
    
    fit_oc[inits]=opt.minimize(all_vects['occlusion']['objective_func'], all_vects['occlusion']['full_init'], method='L-BFGS-B', options={'maxfun':1000000})
    mfts_oc[inits]=fit_oc[inits]['x'][-3:]
    
for inits in inits_set:
    print(inits)
    temp=[]
    for x in all_vects['occlusion']['lt-']:
        temp.append(p_set[1])
        
    all_vects['occlusion']['initial']=temp
    all_vects['occlusion']['full_init']=temp+inits_set[inits]
    
    fit_oc['d'+inits]=opt.minimize(all_vects['occlusion']['objective_func'], all_vects['occlusion']['full_init'], method='L-BFGS-B', options={'maxfun':1000000})
    mfts_oc['d'+inits]=fit_oc['d'+inits]['x'][-3:]
    


with open('../intermediates/occlusion_init_param_exp.pkl', 'wb') as handle:
    pickle.dump(fit_oc, handle, protocol=pickle.HIGHEST_PROTOCOL)