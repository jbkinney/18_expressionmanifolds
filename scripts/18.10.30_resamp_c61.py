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
curves={'c61':'c61r18'}

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
    
def boot_resampler(construct, initializing, cycles):

    """
	Using all_vects, resamples the data medians and fits a curve to the resampled dataset.
	
	Parameters
	----------
	construct : str
		the construct library to be resampled
	
	initializing : float array
		the terms about which randomized values are selected to initialize the fit
	
	cycles : int
		number of cycles of successfully fitted resampled data to require
	
	Returns
	-------
	array of fitted terms
	array of arrays of the indices of the datapoints resampled to be fit
    """
    
    print("resampling "+construct)
    lbasal=all_vects[construct]['lt-']
    linduced=all_vects[construct]['lt+']
    rns=all_vects[construct]['initial']
    
    lres_fits=[]
    lsets_used=[]
    lsets_used_temp=[]
    lpvals=[]
    
    
    cntr=0
    
    while cntr<cycles:
        np.random.seed(cntr)
    
        x_temp=[]
        
        weighter=np.zeros_like(lbasal)
        
        tm=initializing[0]
        tbg=initializing[1]
        et=initializing[2]
        init=[tm,tbg,et]
    
        for j in range(len(lbasal)):
            a=randint(0,len(lbasal)-1)
            weighter[a]+=1
        
        for n in range(len(weighter)):
            if weighter[n]==0: continue
            x_temp.append(rns[n])
        
        
        comp_init=x_temp+init
        lsets_used_temp.append(weighter)
        
        
    
        # Specify data on which to compute objective function
        objective_func = lambda m: comp_error(m,lbasal,linduced,weighter)
        
        fits=opt.minimize(objective_func,comp_init, method='L-BFGS-B')
        lres_fits.append([fits['success'],fits['x'][-3:],fits['x'][:-3]])
        
        if fits['success']: 
            cntr+=1
            
            # a counter to see what run it's on
            sys.stdout.write('\r'+str(cntr))
        
    lthetas=[]
    for ind, res in enumerate(lres_fits):
        #if the fit ended successfully, add the theta values
        if res[0]:
            lthetas.append(res[1])
            lsets_used.append(lsets_used_temp[ind])
            lpvals.append(res[2])
        
    lthetas=np.array(lthetas)
    lpvals=np.array(lpvals)
    
    return lthetas, lsets_used, lpvals

init_param_exp=pickle.load(open('../intermediates/c61_init_param_exp.pkl', 'rb'))

fits={}
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

recalc=True
cycles=100
construct='c61r18'
resamplings={}

if recalc:
    
    lthetas, lsets_used, lpvals = boot_resampler(construct, c61_init, cycles)
    
    resamp=pd.DataFrame(lthetas, columns=['log_tmax','log_tbackground','log_alphap'])
    resamp['P_vals']=[x for x in lpvals]
    resamp['sets_used']=[x for x in lsets_used]
    
    resamp.to_pickle("../intermediates/weighted_"+construct+"_thetas_"+str(cycles)+".pkl")
    
    resamplings[construct]=resamp
