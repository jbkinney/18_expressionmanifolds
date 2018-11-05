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

raw_vals={}

for plate in plate_dat:
    for val in range(len(plate_dat[plate])):
        nam=plate_dat[plate]["glycerol loc"][val]
        
        if nam=='RDM' or nam=='blank':
            continue
            
        if plate_dat[plate]['OK seq?'][val]=='no' or float(plate_dat[plate]['OD600'][val])>0.5:
            continue
        
        if nam not in raw_vals:
            raw_vals[nam]={'0.0':[], '2.5':[], '5.0':[], '10.0':[],
                           '25.0':[], '50.0':[], '125.0':[], '250.0':[]}
        
        if plate_dat[plate]['cAMP'][val] not in raw_vals[nam]:
                raw_vals[nam][plate_dat[plate]['cAMP'][val]]=[]
        
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

curves={'c60':'c60r18.10', 'c61':'c61r18', 'c62':'c62r18.10', 'c63':'c63r18.10', 'c64':'c64r18.10', 
        'c65':'c65r18.10', 'c66':'c66r18.10', 'c71':'c71r18', 'c72':'c72r18.10', 'c76':'c76r18.10', 
        'c81':'c81r18.10', 'c82':'c82r18.10', 'c-':'c-', 'occlusion':'occlusion', 
        'c61r17':'c61r17.cons', 'c41':'gal.c41'}

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
        

resamplings=pd.read_pickle('../intermediates/weighted_conjoined_thetas_100.pkl')

ltsat=np.percentile(resamplings['log_tmax'], 50)

#calcs log of tau basal from same givens as above
def comp_ltaub(lrn,lthetas): 
    """
    Calculates and returns the log of the basal transcription.
    
    Parameters
    ----------
    lrn : float
        The log of the P for this promoter.
        
    lthetas : float array
        lthetas[0] is the log of the background transcription
        lthetas[1] is the log of the alpha_prime (a' = (1 + F*alpha) / (1 + F)) (not used)
        lthetas[2] is the log of the beta_prime (b' = (1 + alpha*beta*F) / (1 + alpha*F)) (not used)
        
    
    Returns
    -------
    float
        The log of the calculated basal transcription,
        basal transc = ts*(P/(1+P))
    """

    P=np.exp(lrn)
    ts=np.exp(ltsat)
    tbg=np.exp(lthetas[0])
    
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
        lthetas[0] is the log of the background transcription
        lthetas[1] is the log of the alpha_prime (a' = (1 + F*alpha) / (1 + F))
        lthetas[2] is the log of the beta_prime (b' = (1 + alpha*beta*F) / (1 + alpha*F))
    
    Returns
    -------
    float
        The log of the calculated induced transcription,
        induced transc = bp*ts*(ap*P/(1+ap*P))+tbg
    """

    P=np.exp(lrn)
    ts=np.exp(ltsat)
    tbg=np.exp(lthetas[0])
    ap=np.exp(lthetas[1])
    bp=np.exp(lthetas[2])

    
    #partition fn
    Z = 1 + P*ap
    
    #total transc
    ti = bp*ts*P*ap/Z + tbg
    
    return np.log(ti)
    
#objective function, given the measured values for log basal transc and log induced transc,
#and a vector for the log of all variables fitted
#returns the squared difference in calculated log induced transc from each measured log induc transc we've got
#plus calculated log basal transc from each measured log basal transc
def comp_error(x, lbasal, linduced, weights):
    
    err=[]
    counter=0
    
    for y in range(len(lbasal)):
        
        if weights[y]==0: continue
        
        #subtract the calculated log tau basal and log tau induced from measured
        errb=lbasal[y]-comp_ltaub(x[counter],x[-3:]) 
        erri=linduced[y]-comp_ltaui(x[counter],x[-3:])
        
        #add squares of those differences
        err.append(((errb**2)+(erri**2))*weights[y])
        
        counter+=1
        
    return sum(err)
    
def boot_resampler(construct, cycles):

    """
	Using all_vects, resamples the data medians and fits a curve to the resampled dataset.
	
	Parameters
	----------
	construct : str
		the construct library to be resampled from all_vects
	
	cycles : int
		number of cycles of successfully fitted resampled data to require
	
	Returns
	-------
	array of fitted terms
	array of arrays of the indices of the datapoints resampled to be fit
    """
    
    print("\nresampling "+construct)
    lbasal=all_vects[construct]['lt-']
    linduced=all_vects[construct]['lt+']
    rns=all_vects[construct]['initial']
    init=all_vects[construct]['full_init'][-3:]
    
    lthetas=[]
    sets_used=[]
    lpvals=[]
    
    cntr=0
    
    while cntr<cycles:
        np.random.seed(cntr)
    
        x_temp=[]
        
        #choose weights for these data points
        weighter=np.zeros_like(lbasal)
    
        for j in range(len(lbasal)):
            a=randint(0,len(lbasal)-1)
            weighter[a]+=1
        
        #for every nonzero weight, add it to the set of parameters fitted
        for n in range(len(weighter)):
            if weighter[n]==0: continue
            x_temp.append(rns[n])
        
        comp_init=x_temp+init
        
        # Specify data on which to compute objective function
        objective_func = lambda m: comp_error(m,lbasal,linduced,weighter)
        
        fits=opt.minimize(objective_func,comp_init, method='L-BFGS-B')
        
        if fits['success']:
            lthetas.append(fits['x'][-3:])
            sets_used.append(weighter)
            lpvals.append(fits['x'][:-3])
            
            cntr+=1
            
            # a counter to see what run it's on
            sys.stdout.write('\r'+str(cntr))
        
    lthetas=np.array(lthetas)
    lpvals=np.array(lpvals)
    
    return lthetas, sets_used, lpvals

#objective function for each library is defined as the comp_error run on whatever is handed to the objective function,
#plus the measured basal and induced median values
#we then optimize the error
def objective_func_creator(lib):
    bas=all_vects[lib]['lt-']
    ind=all_vects[lib]['lt+']
    dumweights=np.ones_like(all_vects[lib]['lt-'])
    lammy=lambda m: comp_error(m,bas,ind,dumweights)
    return lammy
    
#make a new lib of all c61 constructs, r18 and r17 combined
conjoined={'lt+':[], 'lt-':[], 'locs':[]}

for lib in ['c61r18', 'c61r17.cons']:
    for category in ['lt+', 'lt-', 'locs']:
        for cons in all_vects[lib][category]:
            conjoined[category].append(cons)
            
all_vects['c61']=conjoined

#of the libs, make a list to fit, excluding the c- and occlusion ones
fit_set=[]

for lib in all_vects.keys():
    if 'oc' in lib: continue
    if 'c-' in lib: continue
    else:
        fit_set.append(lib)
        
#make objective functions to minimize for all libs to fit
for lib in fit_set:
    all_vects[lib]['objective_func']=objective_func_creator(lib)

#initializing values are set here
tbg=np.log(0.001)
et=np.log(50)
bet=np.log(2)

init=[tbg,et,bet]

for cons in all_vects:
    temp=[]
    for x in all_vects[cons]['lt-']:
        rtemp=np.exp(x)/np.exp(ltsat)
        temp.append(np.log(rtemp))
        
    all_vects[cons]['initial']=temp
    all_vects[cons]['full_init']=temp+init
    
fits={}

for lib in fit_set:
    print(lib)
    fits[lib]=opt.minimize(all_vects[lib]['objective_func'], all_vects[lib]['full_init'], method='L-BFGS-B', options={'maxfun':1000000})

with open('../intermediates/beta_weighted_all_thetas_fit.pkl', 'wb') as handle:
    pickle.dump(fits, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
recalc=True
cycles=100
constructs=['c60r18.10', 'c61r18', 'c62r18.10', 'c63r18.10', 'c64r18.10', 'c65r18.10', 
            'c71r18', 'c72r18.10', 'c81r18.10', 'c82r18.10', 'c61', 'gal.c41']
resamplings={}

if recalc:
    
    for cons in constructs:
        lthetas, sets_used, lpvals = boot_resampler(cons, cycles)
    
        resamp=pd.DataFrame(lthetas, columns=['log_tbackground','log_alphap', 'log_betap'])
        resamp['P_vals']=[x for x in lpvals]
        resamp['sets_used']=[x for x in sets_used]
    
        resamp.to_pickle("../intermediates/beta_weighted_"+cons+"_thetas_"+str(cycles)+".pkl")
    
        resamplings[cons]=resamp
    
    
else:
    
    for cons in constructs:
        resamplings[cons]=pd.read_pickle("../intermediates/beta_weighted_"+cons+"_thetas_"+str(cycles)+".pkl")

