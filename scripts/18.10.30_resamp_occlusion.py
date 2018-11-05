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
    
def boot_resampler_oc(construct, initializing, cycles):

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

    cntr=0
    
    lbasal=all_vects[construct]['lt-']
    linduced=all_vects[construct]['lt+']
    rns=all_vects[construct]['initial']
    
    lres_fits=[]
    lsets_used=[]
    lsets_used_temp=[]
    lpvals=[]
    
    while cntr<cycles:
        np.random.seed(cntr)
        
        x_temp=[]
        
        weighter=np.zeros_like(lbasal)
        
        tm=initializing[0]
        tbg=initializing[1]
        ec=initializing[2]
        init=[tm,tbg,ec]
    
        for j in range(len(lbasal)):
            a=randint(0,len(lbasal)-1)
            weighter[a]+=1
        
        for n in range(len(weighter)):
            if weighter[n]==0: continue
            x_temp.append(rns[n])
            
        comp_init=x_temp+init
        lsets_used_temp.append(weighter)

    
        # Specify data on which to compute objective function
        objective_func = lambda m: comp_error_oc(m,lbasal,linduced, weighter)
        
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
    
all_vects['occlusion']['objective_func']=lambda m: comp_error_oc(m,all_vects['occlusion']['lt-'],all_vects['occlusion']['lt+'],np.ones_like(all_vects['occlusion']['lt-']))

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

fit_oc['occlusion']=opt.minimize(all_vects['occlusion']['objective_func'], all_vects['occlusion']['full_init'], method='L-BFGS-B', options={'maxfun':1000000})
mfts_oc['occlusion']=fit_oc['occlusion']['x'][-3:]

recalc=True
cycles=100
construct='occlusion'
resamplings={}

if recalc:
    
    lthetas, lsets_used, lpvals = boot_resampler_oc(construct, init, cycles)
    
    resamp=pd.DataFrame(lthetas, columns=['log_tmax','log_tbackground','log_F'])
    resamp['P_vals']=[x for x in lpvals]
    resamp['sets_used']=[x for x in lsets_used]
    
    resamp.to_pickle("../intermediates/weighted_"+construct+"_thetas_"+str(cycles)+".pkl")
    
    resamplings[construct]=resamp