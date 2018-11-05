import os
import pylab as py
import pandas as pd
import matplotlib
import numpy as np
import scipy.optimize as opt
import glob
from matplotlib import pyplot as plt
import sys
from random import randint
from random import uniform
import pickle

import pdb

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

curves={'c60':'c60r18.10', 'c61':'c61r18', 'c62':'c62r18.10', 'c63':'c63r18.10', 'c64':'c64r18.10', 
        'c65':'c65r18.10', 'c66':'c66r18.10', 'c71':'c71r18', 'c72':'c72r18.10', 'c76':'c76r18.10', 
        'c81':'c81r18.10', 'c82':'c82r18.10', 'c-':'c-', 'occlusion':'occlusion', 
        'c61r17':'c61r17.cons', 'c41':'gal.c41'}

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
        
fits={}
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
        if np.isnan(corrected_medians[cons]['250.0']): continue
        if np.isnan(corrected_medians[cons]['0.0']): continue
            
        all_vects[lib]['lt+'].append(np.log(corrected_medians[cons]['250.0']))
        all_vects[lib]['lt-'].append(np.log(corrected_medians[cons]['0.0']))
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
    
    
init_param_exp=pickle.load(open('../intermediates/occlusion_init_param_exp.pkl', 'rb'))

fit_oc={}
mfts_oc={}

for inits in init_param_exp:
    mfts_oc[inits]=init_param_exp[inits]['x'][-3:]

tms=[]
tbgs=[]
fs=[]
ps=[]

for cond in mfts_oc:
    tms.append(np.exp(mfts_oc[cond][0]))
    tbgs.append(np.exp(mfts_oc[cond][1]))
    fs.append(np.exp(mfts_oc[cond][2]))
    ps.append(init_param_exp[cond]['x'][:-3])
    
ps_split={}
init_ps=[]

for n,cons in enumerate(all_vects['occlusion']['locs']):
    tmp=[]
    for fit in init_param_exp:
        tmp.append(np.exp(init_param_exp[fit]['x'][n]))
        
    ps_split[cons]=tmp
    
    init_ps.append(np.log(np.percentile(tmp,50)))
    
tsat_oc=np.log(np.percentile(tms,50))
tbg_oc=np.log(np.percentile(tbgs,50))
F_oc=np.log(np.percentile(fs,50))

all_vects['occlusion']['initial']=init_ps
all_vects['occlusion']['full_init']=init_ps+[tsat_oc,tbg_oc,F_oc]

init=[tsat_oc,tbg_oc,F_oc]

all_vects['occlusion']['objective_func']=lambda m: comp_error_oc(m,all_vects['occlusion']['lt-'],all_vects['occlusion']['lt+'],np.ones_like(all_vects['occlusion']['lt-']))

fit_oc['occlusion']=opt.minimize(all_vects['occlusion']['objective_func'], all_vects['occlusion']['full_init'], method='L-BFGS-B', options={'maxfun':1000000})
mfts_oc['occlusion']=fit_oc['occlusion']['x'][-3:]
fits['occlusion']=fit_oc['occlusion']

        
#necessary functions for c61 fitting
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
        lthetas[2] is the log of the alpha_prime (a' = (1 + F*alpha) / (1 + F)) (not used)
        lthetas[3] is the log of the beta_prime (b' = (1 + alpha*beta*F) / (1 + alpha*F)) (not used)
        
    
    Returns
    -------
    float
        The log of the calculated basal transcription,
        basal transc = ts*(P/(1+P))
    """

    P=np.exp(lrn)
    ts=np.exp(lthetas[0])
    tbg=np.exp(lthetas[1])
    
    Z=1+P
    
    #total transc
    tb=(ts*P/Z)+tbg
    
    return np.log(tb)
    
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
        lthetas[2] is the log of the alpha_prime (a' = (1 + F*alpha) / (1 + F))
    
    Returns
    -------
    float
        The log of the calculated induced transcription,
        induced transc = ts*(ap*P/(1+ap*P))+tbg
    """

    P=np.exp(lrn)
    ts=np.exp(lthetas[0])
    tbg=np.exp(lthetas[1])
    ap=np.exp(lthetas[2])

    
    #partition fn
    Z = 1 + P*ap
    
    #total transc
    ti = ts*P*ap/Z + tbg
    
    return np.log(ti)
    
#objective function, given the measured values for log basal transc and log induced transc,
#and a vector for the log of all variables fitted
#returns the squared difference in calculated log induced transc from each measured log induc transc we've got
#plus calculated log basal transc from each measured log basal transc
def comp_error(x, lbasal, linduced, weights):
    
    mocker=[]
    counter=0
    err=[]
    
    for count in range(len(weights)):
        if weights[count]==0:
            mocker.append(0)
        else:
            mocker.append(x[counter])
            counter+=1
    
    for y in range(len(lbasal)):
        
        if mocker[y]==0: continue
        
        #subtract the calculated log tau basal and log tau induced from measured
        errb=lbasal[y]-comp_ltaub(mocker[y],x[-3:]) 
        erri=linduced[y]-comp_ltaui(mocker[y],x[-3:])
        
        #add squares of those differences
        err.append(((errb**2)+(erri**2))*weights[y])
        
    return sum(err)

init_param_exp=pickle.load(open('../intermediates/c61_init_param_exp.pkl', 'rb'))


mfts={}

for inits in init_param_exp:
    mfts[inits]=init_param_exp[inits]['x'][-3:]

c61_tms=[]
c61_tbgs=[]
c61_aps=[]
c61_ps=[]

for cond in mfts:
    c61_tms.append(np.exp(mfts[cond][0]))
    c61_tbgs.append(np.exp(mfts[cond][1]))
    c61_aps.append(np.exp(mfts[cond][2]))
    c61_ps.append(init_param_exp[cond]['x'][:-3])
    
c61_tsat_fit=np.log(np.percentile(c61_tms,50))
c61_tbg_fit=np.log(np.percentile(c61_tbgs,50))
c61_ap_fit=np.log(np.percentile(c61_aps,50))

c61_init=[c61_tsat_fit, c61_tbg_fit, c61_ap_fit]

c61_ps_split={}
c61_init_ps=[]

for n,cons in enumerate(all_vects['c61r18']['locs']):
    tmp=[]
    for fit in init_param_exp:
        tmp.append(np.exp(init_param_exp[fit]['x'][n]))
        
    c61_ps_split[cons]=tmp
    
    c61_init_ps.append(np.log(np.percentile(tmp,50)))
    
all_vects['c61r18']['initial']=c61_init_ps
all_vects['c61r18']['full_init']=c61_init_ps+c61_init

all_vects['c61r18']['objective_func']=lambda m: comp_error(m,all_vects['c61r18']['lt-'],all_vects['c61r18']['lt+'],np.ones_like(all_vects['c61r18']['lt-']))

fits['c61r18']=opt.minimize(all_vects['c61r18']['objective_func'], all_vects['c61r18']['full_init'], method='L-BFGS-B', options={'maxfun':1000000})
mfts['c61r18']=fits['c61r18']['x'][-3:]

        

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
    
curves={'c60':'c60r18.10', 'c61':'c61r18', 'c62':'c62r18.10', 'c63':'c63r18.10', 'c64':'c64r18.10', 
        'c65':'c65r18.10', 'c66':'c66r18.10', 'c71':'c71r18', 'c72':'c72r18.10', 'c76':'c76r18.10', 
        'c81':'c81r18.10', 'c82':'c82r18.10'}
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


init_param_exp=pickle.load(open('../intermediates/conj_init_param_exp.pkl', 'rb'))

for inits in init_param_exp:
    mfts[inits]=init_param_exp[inits]['x']

conj_terms=[]

for cond in init_param_exp:
    conj_terms.append(init_param_exp[cond]['x'])

conj_terms_split={}
conj_inits=[]

columns=['log_tmax', 
         'log_tbg_'+all_vects['conj']['curve_list'][0],
         'log_tbg_'+all_vects['conj']['curve_list'][1],
         'log_tbg_'+all_vects['conj']['curve_list'][2],
         'log_tbg_'+all_vects['conj']['curve_list'][3],
         'log_tbg_'+all_vects['conj']['curve_list'][4],
         'log_tbg_'+all_vects['conj']['curve_list'][5],
         'log_tbg_'+all_vects['conj']['curve_list'][6],
         'log_tbg_'+all_vects['conj']['curve_list'][7],
         'log_tbg_'+all_vects['conj']['curve_list'][8],
         'log_tbg_'+all_vects['conj']['curve_list'][9],
         'log_tbg_'+all_vects['conj']['curve_list'][10],
         'log_tbg_'+all_vects['conj']['curve_list'][11],
         'log_alphap_'+all_vects['conj']['curve_list'][0],
         'log_alphap_'+all_vects['conj']['curve_list'][1],
         'log_alphap_'+all_vects['conj']['curve_list'][2],
         'log_alphap_'+all_vects['conj']['curve_list'][3],
         'log_alphap_'+all_vects['conj']['curve_list'][4],
         'log_alphap_'+all_vects['conj']['curve_list'][5],
         'log_alphap_'+all_vects['conj']['curve_list'][6],
         'log_alphap_'+all_vects['conj']['curve_list'][7],
         'log_alphap_'+all_vects['conj']['curve_list'][8],
         'log_alphap_'+all_vects['conj']['curve_list'][9],
         'log_alphap_'+all_vects['conj']['curve_list'][10],
         'log_alphap_'+all_vects['conj']['curve_list'][11]]

for n,cons in enumerate(all_vects['conj']['locs']+columns):
    tmp=[]
    for fit in init_param_exp:
        tmp.append(np.exp(init_param_exp[fit]['x'][n]))
        
    conj_terms_split[cons]=tmp
    
    conj_inits.append(np.log(np.percentile(tmp,50)))
    
all_vects['conj']['initial']=[np.log(np.percentile(conj_terms_split[x],50)) for x in all_vects['conj']['locs']]
inits=[np.log(np.percentile(conj_terms_split[x],50)) for x in columns]
init_conj=inits

all_vects['conj']['full_init']=all_vects['conj']['initial']+inits


fits['conj']=opt.minimize(all_vects['conj']['objective_func'], all_vects['conj']['full_init'], method='L-BFGS-B', options={'maxfun':1000000})



with open('../intermediates/weighted_all_thetas_fit.pkl', 'wb') as handle:
    pickle.dump(fits, handle, protocol=pickle.HIGHEST_PROTOCOL)


