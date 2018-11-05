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
        
curves={'c60':'c60r18.10', 'c61':'c61r18', 'c62':'c62r18.10', 'c63':'c63r18.10', 'c64':'c64r18.10', 
        'c65':'c65r18.10', 'c66':'c66r18.10', 'c71':'c71r18', 'c72':'c72r18.10', 'c76':'c76r18.10', 
        'c81':'c81r18.10', 'c82':'c82r18.10', 'c-':'c-', 'occlusion':'occlusion', 'c61r17':'c61r17.cons'}

all_vects={}
cAMP_concs=['250.0','125.0','50.0','25.0','10.0','5.0','2.5']

for spac in ['c71','occlusion']:
    lib=curves[spac]
    all_vects[lib]={'ltplus':{'2.5':[], '5.0':[], '10.0':[],
                           '25.0':[], '50.0':[], '125.0':[], '250.0':[]},
                    'ltminus':{'2.5':[], '5.0':[], '10.0':[],
                           '25.0':[], '50.0':[], '125.0':[], '250.0':[]},
                    'locs':{'2.5':[], '5.0':[], '10.0':[],
                           '25.0':[], '50.0':[], '125.0':[], '250.0':[]}}
    for cons in corrected_medians:

        if cons not in library_groups[lib]['all']: continue
        if cons in library_groups[lib]['outliers']: continue
        if np.isnan(corrected_medians[cons]['0.0']): continue
            
        for conc in cAMP_concs:
            if np.isnan(corrected_medians[cons][conc]): continue
            all_vects[lib]['ltplus'][conc].append(np.log(corrected_medians[cons][conc]))
            all_vects[lib]['ltminus'][conc].append(np.log(corrected_medians[cons]['0.0']))
            all_vects[lib]['locs'][conc].append(cons)
            

resamplings={}
resamplings['occlusion']=pd.read_pickle('../intermediates/weighted_occlusion_thetas_100.pkl')
resamplings['conjoined']=pd.read_pickle('../intermediates/weighted_conjoined_thetas_100.pkl')

ltsat=np.log(np.percentile(np.exp(resamplings['conjoined']['log_tmax']),50))
ltbg_c71=np.log(np.percentile(np.exp(resamplings['conjoined']['log_tbg_c71r18']),50))
lalphap_c71=np.log(np.percentile(np.exp(resamplings['conjoined']['log_alphap_c71r18']),50))
ltbg_ocl=np.log(np.percentile(np.exp(resamplings['occlusion']['log_tbackground']),50))
lF_ocl=np.log(np.percentile(np.exp(resamplings['occlusion']['log_F']),50))

lalph_c71=np.log((np.exp(lalphap_c71)*(1+np.exp(lF_ocl))-1)/np.exp(lF_ocl))

#calcs log of tau basal from same givens as above
def comp_ltaub(lP,lF): 
    """
    Calculates and returns the log of the basal transcription.
    
    Parameters
    ----------
    lP : float
        The log of the P for this promoter.
        
    lF : float
        The log of F for this concentration.
        
    
    Returns
    -------
    float
        The log of the calculated basal transcription,
        basal transc = ts*(P/(1+P))+tbg
    """

    P=np.exp(lP)
    F=np.exp(lF)
    ts=np.exp(ltsat)
    tbg=np.exp(ltbg_c71)
    alph=np.exp(lalph_c71)
    
    Z=1+P
    
    #total transc
    tb=(ts*P/Z)+tbg
    
    return np.log(tb)
    
def comp_ltaui(lP,lF): 
    """
    Calculates and returns the log of the induced transcription.
    
    Parameters
    ----------
    lP : float
        The log of the P for this promoter.
        
    lF : float
        The log of F for this concentration.
    
    Returns
    -------
    float
        The log of the calculated induced transcription,
        induced transc = ts*((P+alph*P)/(1+F+P+alph*P))+tbg
    """

    P=np.exp(lP)
    F=np.exp(lF)
    ts=np.exp(ltsat)
    tbg=np.exp(ltbg_c71)
    alph=np.exp(lalph_c71)

    
    #partition fn
    Z = 1+F+P+alph*F*P
    
    #total transc
    ti = ts*(P+alph*F*P)/Z + tbg
    
    return np.log(ti)
    
#returns calculated log induced transc given logs of:
#rn=RNAP binding weight to this sequence, tm=max WT lac transc, tbg=background transc, 
#ec=CRP binding weight
def comp_ltaui_oc(lP,lF): 
    """
    Calculates and returns the log of the induced transcription for occluded constructs.
    
    Parameters
    ----------
    lP : float
        The log of the P for this promoter.
        
    lF : float
        The log of F for this concentration.
    
    Returns
    -------
    float
        The log of the calculated induced transcription,
        induced transc = ts*((P)/(1+F+P))+tbg
    """
    
    P=np.exp(lP)
    F=np.exp(lF)
    ts=np.exp(ltsat)
    tbg=np.exp(ltbg_ocl)

    
    #partition fn
    Z=1+F+P
    
    #total transc
    ti=(ts*P/Z)+tbg
    
    return np.log(ti)

#calcs log of tau basal from same givens as above
def comp_ltaub_oc(lP,lF): 
    """
    Calculates and returns the log of the basal transcription for occluded constructs.
    
    Parameters
    ----------
    lP : float
        The log of the P for this promoter.
        
    lF : float
        The log of F for this concentration.
    
    Returns
    -------
    float
        The log of the calculated basal transcription,
        basal transc = ts*((P)/(1+P))+tbg
    """
    
    P=np.exp(lP)
    F=np.exp(lF)
    ts=np.exp(ltsat)
    tbg=np.exp(ltbg_ocl)

    
    #partition fn
    Z=1+P
    
    #total transc
    tb=(ts*P/Z)+tbg
    
    return np.log(tb)
    
def comp_error_oc(x, conc, weights):
    lbasal=all_vects['occlusion']['ltminus'][conc]
    linduced=all_vects['occlusion']['ltplus'][conc]
    
    #ltsat and ltbg_ocl are taken from earlier fitting
    
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
        errb=lbasal[y]-comp_ltaub_oc(mocker[y],x[-1]) 
        erri=linduced[y]-comp_ltaui_oc(mocker[y],x[-1])
        
        #add squares of those differences
        err.append(((errb**2)+(erri**2))*weights[y])
        
    return sum(err)
    
def comp_error_c71(x, conc, weights):
    lbasal=all_vects['c71r18']['ltminus'][conc]
    linduced=all_vects['c71r18']['ltplus'][conc]
    
    #ltsat and ltbg_ocl are taken from earlier fitting
    
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
        errb=lbasal[y]-comp_ltaub(mocker[y],x[-1]) 
        erri=linduced[y]-comp_ltaui(mocker[y],x[-1])
        
        #add squares of those differences
        err.append(((errb**2)+(erri**2))*weights[y])
        
    return sum(err)
    
dum_weights={}
for construct in all_vects:
    dum_weights[construct]={}
    for conc in all_vects[construct]['ltplus']:
        dum_weights[construct][conc]=np.ones_like(all_vects[construct]['ltplus'][conc])
        
for construct in all_vects:
    all_vects[construct]['inits']={}
    
    for conc in all_vects[construct]['locs']:
        all_vects[construct]['inits'][conc]=[]
        
        for transc in all_vects[construct]['ltminus'][conc]:
            tbas=np.exp(transc)
            all_vects[construct]['inits'][conc].append(np.log(tbas/np.exp(ltsat)))
            
        all_vects[construct]['inits'][conc].append(np.log(30))
        
def obj_func_oc(conc,weights):
    """
    Creates an objective function for occlusion libraries 
    given a concentration in the all_vects dictionary.
    
    Parameters
    ----------
    conc : string 
        concentration of cAMP, in uM
        
    Returns
    -------
    lambda function
        function that calculates the squared error for parameters from the measured vals
    """
    
    lammy=lambda m: comp_error_oc(m,conc,weights)
    
    return lammy

def obj_func_c71(conc,weights):
    """
    Creates an objective function for c71r18 libraries 
    given a concentration in the all_vects dictionary.
    
    Parameters
    ----------
    conc : string 
        concentration of cAMP, in uM
        
    Returns
    -------
    lambda function
        function that calculates the squared error for parameters from the measured vals
    """
    
    lammy=lambda m: comp_error_c71(m,conc,weights)
    
    return lammy

all_vects['occlusion']['objective_func']={}
all_vects['c71r18']['objective_func']={}

#generate the objective functions of the known cAMP concentrations
for conc in all_vects['occlusion']['ltminus']:
    all_vects['occlusion']['objective_func'][conc]=obj_func_oc(conc,dum_weights['occlusion'][conc])
    
for conc in all_vects['c71r18']['ltminus']:
    all_vects['c71r18']['objective_func'][conc]=obj_func_c71(conc,dum_weights['c71r18'][conc])
    
fits={}
fit_terms={}

for construct in all_vects:
    fits[construct]={}
    fit_terms[construct]={}
    for conc in all_vects[construct]['objective_func']:
        fits[construct][conc]=opt.minimize(all_vects[construct]['objective_func'][conc], all_vects[construct]['inits'][conc], method='L-BFGS-B')
        fit_terms[construct][conc]=fits[construct][conc]['x']
        
with open('../intermediates/oc_and_c71_weighted_cAMPdil_thetas_fit.pkl', 'wb') as handle:
    pickle.dump(fits, handle, protocol=pickle.HIGHEST_PROTOCOL)

def boot_resampler_oc(cycles, conc):
    print("\nresampling occlusion at "+conc+'uM')
    cntr=0
    sets_used=[]
    lthetas=[]
    Fs_fit=[]
    
    #continue attempting fits until there's cycles number of successful fits
    while cntr<cycles:
        np.random.seed(cntr)
    
        #start with a vector of weights with each equal to zero
        weights=np.zeros_like(all_vects['occlusion']['ltminus'][conc])
        
        #choose the weight indices randomly and add one to the weight at that index
        for j in range(len(all_vects['occlusion']['locs'][conc])):
            a=randint(0,len(all_vects['occlusion']['locs'][conc])-1)
            weights[a]+=1
            
        #the initial values for unchosen datapoints need to be removed
        comp_init=[]
        
        for ind, weight in enumerate(weights):
            if weight==0: continue
            
            comp_init.append(all_vects['occlusion']['inits'][conc][ind])
            
        #initial value for F set to 30
        comp_init.append(np.log(30))
        
        #objective function using the occlusion error function
        obj_function=lambda m: comp_error_oc(m,conc,weights)
        
        #minimize the objective function using the selected initial values
        fits=opt.minimize(obj_function, comp_init, method='L-BFGS-B')
        
        #if this fit ends successfully: 
        #add the fitted terms to lthetas and the weights to sets_used, then increment the counter
        if fits['success']:
            
            lthetas.append(fits['x'][:-1])
            Fs_fit.append(fits['x'][-1])
            sets_used.append(weights)
            
            cntr+=1
            
            # a counter to see what run it's on
            sys.stdout.write('\r'+str(cntr))
            
    lthetas=np.array(lthetas)
    
    return Fs_fit, lthetas, sets_used
    
recalc=True
cycles=100
concentrations=['2.5', '5.0','10.0', '25.0', '50.0', '125.0', '250.0']
resamplings_oc={}

if recalc:
    
    for conc in concentrations:
        Fs_fit, lthetas, lsets_used = boot_resampler_oc(cycles, conc)
    
        resamp=pd.DataFrame(Fs_fit, columns=['log_F'])
        resamp['lp_vals']=[x for x in lthetas]
        resamp['sets_used']=[x for x in lsets_used]
        
        resamplings_oc[conc]=resamp
    
    with open('../intermediates/occlusion_weighted_cAMPdil_thetas_100.pkl', 'wb') as handle:
        pickle.dump(resamplings_oc, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
else:
    
    resamplings_oc=pickle.load(open('../intermediates/occlusion_weighted_cAMPdil_thetas_100.pkl', 'rb'))

def boot_resampler_c71(cycles, conc):
    print("\nresampling c71r18 at "+conc+'uM')
    cntr=0
    sets_used=[]
    lthetas=[]
    Fs_fit=[]
    
    #continue attempting fits until there's cycles number of successful fits
    while cntr<cycles:
        np.random.seed(cntr)
    
        #start with a vector of weights with each equal to zero
        weights=np.zeros_like(all_vects['c71r18']['ltminus'][conc])
        
        #choose the weight indices randomly and add one to the weight at that index
        for j in range(len(all_vects['c71r18']['locs'][conc])):
            a=randint(0,len(all_vects['c71r18']['locs'][conc])-1)
            weights[a]+=1
            
        #the initial values for unchosen datapoints need to be removed
        comp_init=[]
        
        for ind, weight in enumerate(weights):
            if weight==0: continue
            
            comp_init.append(all_vects['c71r18']['inits'][conc][ind])
            
        #initial value for F set to 30
        comp_init.append(np.log(30))
        
        #objective function using the occlusion error function
        obj_function=lambda m: comp_error_c71(m,conc,weights)
        
        #minimize the objective function using the selected initial values
        fits=opt.minimize(obj_function, comp_init, method='L-BFGS-B')
        
        #if this fit ends successfully: 
        #add the fitted terms to lthetas and the weights to sets_used, then increment the counter
        if fits['success']:
            
            lthetas.append(fits['x'][:-1])
            Fs_fit.append(fits['x'][-1])
            sets_used.append(weights)
            
            cntr+=1
            
            # a counter to see what run it's on
            sys.stdout.write('\r'+str(cntr))
            
    lthetas=np.array(lthetas)
    
    return Fs_fit, lthetas, sets_used
    
recalc=True
cycles=100
concentrations=['2.5', '5.0','10.0', '25.0', '50.0', '125.0', '250.0']
resamplings_c71={}

if recalc:
    
    for conc in concentrations:
        Fs_fit, lthetas, lsets_used = boot_resampler_c71(cycles, conc)
    
        resamp=pd.DataFrame(Fs_fit, columns=['log_F'])
        resamp['lp_vals']=[x for x in lthetas]
        resamp['sets_used']=[x for x in lsets_used]
        
        resamplings_c71[conc]=resamp
    
    with open('../intermediates/c71r18_weighted_cAMPdil_thetas_100.pkl', 'wb') as handle:
        pickle.dump(resamplings_c71, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
else:
    
    resamplings_c71=pickle.load(open('../intermediates/occlusion_weighted_cAMPdil_thetas_100.pkl', 'rb'))

