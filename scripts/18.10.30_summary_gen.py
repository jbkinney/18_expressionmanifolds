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
from pandas import ExcelWriter

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

writer = ExcelWriter('../results.xlsx')

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
        
outliers={}

for lib in library_groups:
    for cons in library_groups[lib]['all']:
        outliers[cons]=cons in library_groups[lib]['outliers']
        
#generate the overall summary for all constructs used, from the list of constructs in seq_sum
corrected_meas_sum={'name':[],'location':[], 'log_t+':[], 'dlog_t+':[], 'log_t-':[], 'dlog_t-':[], 
                    'num_t+':[], 'num_t-':[], 'outlier':[], 'spacing':[], 'sequence':[]}

for n in range(len(seq_sum['location'])):
    loc=seq_sum['location'][n]
    corrected_meas_sum['location'].append(loc)
    corrected_meas_sum['name'].append(cons_names['short_name'][loc])
    corrected_meas_sum['log_t+'].append(np.log(corrected_medians[loc]['250.0']))
    if np.isnan(corrected_medians[loc]['250.0']):
        corrected_meas_sum['dlog_t+'].append(np.NaN)
    else:
        corrected_meas_sum['dlog_t+'].append((np.std(np.log(corrected_vals[loc]['250.0']), ddof=1)))
    corrected_meas_sum['log_t-'].append(np.log(corrected_medians[loc]['0.0']))
    if np.isnan(corrected_medians[loc]['0.0']):
        corrected_meas_sum['dlog_t-'].append(np.NaN)
    else:
        corrected_meas_sum['dlog_t-'].append((np.std(np.log(corrected_vals[loc]['0.0']), ddof=1)))
    corrected_meas_sum['num_t+'].append(len(corrected_vals[loc]['250.0']))
    corrected_meas_sum['num_t-'].append(len(corrected_vals[loc]['0.0']))
    corrected_meas_sum['spacing'].append(seq_sum['spacing'][n])
    corrected_meas_sum['sequence'].append(seq_sum['sequence'][n])
    if loc in outliers:
        corrected_meas_sum['outlier'].append(outliers[loc])
    else:
        corrected_meas_sum['outlier'].append(np.NaN)
        
#save the summary of all the constructs
corrected_sum_df = pd.DataFrame(corrected_meas_sum, index=corrected_meas_sum['name'], columns=['location', 'log_t+', 'dlog_t+', 'log_t-', 'dlog_t-', 'num_t+', 'num_t-', 'outlier', 'spacing', 'sequence'])
corrected_sum_df.to_csv('../summaries/measurements_summary.txt', sep='\t', index_label='name')
corrected_sum_df.to_excel(writer, 'measurements_summary', index_label='name')


#cAMP summary uses all concentrations for occlusion and c71r18
cAMP_sum={'name':[], 'location':[], 'log_t_250.0':[], 'dlog_t_250.0':[], 'log_t_125.0':[], 'dlog_t_125.0':[],
          'log_t_50.0':[], 'dlog_t_50.0':[], 'log_t_25.0':[], 'dlog_t_25.0':[],
          'log_t_10.0':[], 'dlog_t_10.0':[], 'log_t_5.0':[], 'dlog_t_5.0':[],
          'log_t_2.5':[], 'dlog_t_2.5':[], 'log_t_0.0':[], 'dlog_t_0.0':[]}

for cons in library_groups['c71r18']['all']:
    if cons not in corrected_medians: continue
    if cons in library_groups['c71r18']['outliers']: continue
    cAMP_sum['name'].append(cons_names['short_name'][cons])
    cAMP_sum['location'].append(cons)
    
    for conc in ['250.0', '125.0', '50.0', '25.0', '10.0', '5.0', '2.5', '0.0']:
        cAMP_sum['log_t_'+conc].append(np.log(corrected_medians[cons][conc]))
        if np.isnan(corrected_medians[cons][conc]):
            cAMP_sum['dlog_t_'+conc].append(np.NaN)
        else:
            cAMP_sum['dlog_t_'+conc].append(np.std(np.log(corrected_vals[cons][conc]), ddof=1))
            
for cons in library_groups['occlusion']['all']:
    if cons in library_groups['occlusion']['outliers']: continue
    if cons not in corrected_medians: continue
    cAMP_sum['name'].append(cons_names['short_name'][cons])
    cAMP_sum['location'].append(cons)
    
    for conc in ['250.0', '125.0', '50.0', '25.0', '10.0', '5.0', '2.5', '0.0']:
        cAMP_sum['log_t_'+conc].append(np.log(corrected_medians[cons][conc]))
        if np.isnan(corrected_medians[cons][conc]):
            cAMP_sum['dlog_t_'+conc].append(np.NaN)
        else:
            cAMP_sum['dlog_t_'+conc].append(np.std(np.log(corrected_vals[cons][conc]), ddof=1))
            
#generate the summary of just c71 and occlusion cAMP dilution data
cAMP_sum_df = pd.DataFrame(cAMP_sum, index=cAMP_sum['name'], columns=['location', 'log_t_250.0', 'dlog_t_250.0', 'log_t_125.0', 'dlog_t_125.0', 'log_t_50.0', 'dlog_t_50.0', 'log_t_25.0', 'dlog_t_25.0', 'log_t_10.0', 'dlog_t_10.0', 'log_t_5.0', 'dlog_t_5.0', 'log_t_2.5', 'dlog_t_2.5', 'log_t_0.0', 'dlog_t_0.0'])
cAMP_sum_df.to_csv('../summaries/cAMP_summary.txt', sep='\t', index_label='name')
cAMP_sum_df.to_excel(writer, 'cAMP_summary', index_label='name')

#building up all_vects for the occlusion, c61, and conj summary
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
        

#weighted_all_thetas_fit.pkl has the fit using all data points for occlusion, c61, and conjoined
fits=pickle.load(open('../intermediates/weighted_all_thetas_fit.pkl', 'rb'))
fit_oc={}
mfts_oc={}
fit_oc['occlusion']=fits['occlusion']
mfts_oc['occlusion']=fit_oc['occlusion']['x'][-3:]

resamplings={}
resamplings['occlusion']=pd.read_pickle("../intermediates/weighted_occlusion_thetas_100.pkl")
resamplings['c61r18']=pd.read_pickle("../intermediates/weighted_c61r18_thetas_100.pkl")

    
oc_sum={'run':[], 'log_tsat':[], 'log_tbg':[], 'log_F':[], 'lPs':[], 'weights':[]}

oc_sum['run'].append('fit')
oc_sum['log_tsat'].append(mfts_oc['occlusion'][0])
oc_sum['log_tbg'].append(mfts_oc['occlusion'][1])
oc_sum['log_F'].append(mfts_oc['occlusion'][2])
oc_sum['lPs'].append(fit_oc['occlusion']['x'][:-3])
oc_sum['weights'].append(np.ones_like(all_vects['occlusion']['lt-']))


#
for run, n in enumerate(resamplings['occlusion']['log_tmax']):
    oc_sum['run'].append('samp_%02d'%(run))
    oc_sum['log_tsat'].append(resamplings['occlusion']['log_tmax'][run])
    oc_sum['log_tbg'].append(resamplings['occlusion']['log_tbackground'][run])
    oc_sum['log_F'].append(resamplings['occlusion']['log_F'][run])
    oc_sum['lPs'].append(resamplings['occlusion']['P_vals'][run])
    oc_sum['weights'].append(resamplings['occlusion']['sets_used'][run])
    
for n, cons in enumerate(all_vects['occlusion']['locs']):
    oc_sum[cons]=[]
    oc_sum[cons+'_weight']=[]
    
for run, m in enumerate(oc_sum['weights']):
    count=0
        
    for n, cons in enumerate(all_vects['occlusion']['locs']):
        weight=oc_sum['weights'][run][n]
        oc_sum[cons+'_weight'].append(weight)
        if weight==0:
            oc_sum[cons].append(np.nan)
            continue
        oc_sum[cons].append(oc_sum['lPs'][run][count])
        count+=1
        
cols={x for x in oc_sum.keys()}
cols.remove('run')
cols.remove('log_tsat')
cols.remove('log_tbg')
cols.remove('log_F')
cols.remove('lPs')
cols.remove('weights')

oc_sum_df = pd.DataFrame(oc_sum, index=oc_sum['run'], columns=['log_tsat', 'log_tbg', 'log_F']+[x for x in sorted(cols)])
oc_sum_df.to_csv('../summaries/occlusion_resamp.txt', sep='\t', index_label='run')
oc_sum_df.to_excel(writer, 'occlusion_resamp', index_label='run')

        
mfts={}
mfts['c61r18']=fits['c61r18']['x'][-3:]



c61_sum={'run':[], 'log_tsat':[], 'log_tbg':[], 'log_alphap':[], 'lPs':[], 'weights':[]}

c61_sum['run'].append('fit')
c61_sum['log_tsat'].append(mfts['c61r18'][0])
c61_sum['log_tbg'].append(mfts['c61r18'][1])
c61_sum['log_alphap'].append(mfts['c61r18'][2])
c61_sum['lPs'].append(fits['c61r18']['x'][:-3])
c61_sum['weights'].append(np.ones_like(all_vects['c61r18']['lt-']))


#
for run, n in enumerate(resamplings['c61r18']['log_tmax']):
    c61_sum['run'].append('samp_%02d'%(run))
    c61_sum['log_tsat'].append(resamplings['c61r18']['log_tmax'][run])
    c61_sum['log_tbg'].append(resamplings['c61r18']['log_tbackground'][run])
    c61_sum['log_alphap'].append(resamplings['c61r18']['log_alphap'][run])
    c61_sum['lPs'].append(resamplings['c61r18']['P_vals'][run])
    c61_sum['weights'].append(resamplings['c61r18']['sets_used'][run])
    
for n, cons in enumerate(all_vects['c61r18']['locs']):
    c61_sum[cons]=[]
    c61_sum[cons+'_weight']=[]
    
for run, m in enumerate(c61_sum['weights']):
    count=0
        
    for n, cons in enumerate(all_vects['c61r18']['locs']):
        weight=c61_sum['weights'][run][n]
        c61_sum[cons+'_weight'].append(weight)
        if weight==0:
            c61_sum[cons].append(np.nan)
            continue
        c61_sum[cons].append(c61_sum['lPs'][run][count])
        count+=1
        
cols={x for x in c61_sum.keys()}
cols.remove('run')
cols.remove('log_tsat')
cols.remove('log_tbg')
cols.remove('log_alphap')
cols.remove('lPs')
cols.remove('weights')

c61_sum_df = pd.DataFrame(c61_sum, index=c61_sum['run'], columns=['log_tsat', 'log_tbg', 'log_alphap']+[x for x in sorted(cols)])
c61_sum_df.to_csv('../summaries/c61_resamp.txt', sep='\t', index_label='run')
c61_sum_df.to_excel(writer, 'c61_resamp', index_label='run')

        
#resampled set for conjoined fit
resamplings=pd.read_pickle('../intermediates/weighted_conjoined_thetas_100.pkl')

conj_curves=['c60', 'c61', 'c62', 'c63', 'c64', 'c65', 'c66', 'c71', 'c72', 'c76', 'c81', 'c82']
tsat_set=[np.log(1),np.log(10),np.log(100)]
tbg_set=[np.log(0.0001),np.log(0.001),np.log(0.01)]
alphp_set=[np.log(50),np.log(500),np.log(5000)]
p_set=[np.log(0.01),np.log(0.1)]
    
fit_set=[]
for lib in conj_curves:
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

weight_test=[]

for lib in all_vects['conj']['curve_list']:
    weight_test.append(np.ones_like(all_vects[lib]['lt-']))

conj_sum={'run':[]}
for val in resamplings.keys():
    conj_sum[val]=[]
    
fit_vals=fits['conj']['x']
fit_terms=['log_tmax',
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

#initial fit is added to all columns
conj_sum['run'].append('fit')
for n in range(len(fit_terms)):
    conj_sum[fit_terms[len(fit_terms)-(n+1)]].append(fit_vals[-(n+1)])

conj_sum['P_vals'].append(fit_vals[:-len(fit_terms)])
conj_sum['sets_used'].append(weight_test)


for run, n in enumerate(resamplings['log_tmax']):
    conj_sum['run'].append('samp_%02d'%(run))
    
    for term in fit_terms:
        conj_sum[term].append(resamplings[term][run])
        
    conj_sum['P_vals'].append(resamplings['P_vals'][run])
    conj_sum['sets_used'].append(resamplings['sets_used'][run])

    
cons_cols_list=[]
for cons in all_vects['conj']['locs']:
    conj_sum[cons+'_log_P']=[]
    conj_sum[cons+'_weight']=[]
    cons_cols_list.append(cons+'_log_P')
    cons_cols_list.append(cons+'_weight')

for run, m in enumerate(conj_sum['sets_used']):
    count=0
    
    weight_unpk=[]
    for weighter in conj_sum['sets_used'][run]:
        for n in weighter:
            weight_unpk.append(n)
    
    for n, cons in enumerate(all_vects['conj']['locs']):
        weight=weight_unpk[n]
        conj_sum[cons+'_weight'].append(weight)
        if weight==0:
            conj_sum[cons+'_log_P'].append(np.nan)
            continue
        conj_sum[cons+'_log_P'].append(conj_sum['P_vals'][run][count])
        count+=1
        
conj_sum_df=pd.DataFrame(conj_sum, index=conj_sum['run'], columns=fit_terms+cons_cols_list)
conj_sum_df.to_csv('../summaries/conjoined_resamp.txt', sep='\t', index_label='run')
conj_sum_df.to_excel(writer, 'conjoined_resamp', index_label='run')


#the cAMP summary requires a different structure for all_vects
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
     
#constructing a set of weights for each concentration & construct
dum_weights={}
for construct in all_vects:
    dum_weights[construct]={}
    for conc in all_vects[construct]['ltplus']:
        dum_weights[construct][conc]=np.ones_like(all_vects[construct]['ltplus'][conc])
        
#fit_terms loads in the fit using all data points for each conc for c71 and occ
fit_terms=pickle.load(open('../intermediates/oc_and_c71_weighted_cAMPdil_thetas_fit.pkl', 'rb'))

#load in the resampling summaries for all concentrations for c71 and occlusion
resamplings_oc=pickle.load(open('../intermediates/occlusion_weighted_cAMPdil_thetas_100.pkl', 'rb'))
resamplings_c71=pickle.load(open('../intermediates/c71r18_weighted_cAMPdil_thetas_100.pkl', 'rb'))

for conc in resamplings_c71:
    temp_sum={'run':[], 'log_F':[], 'lPs':[], 'weights':[]}
    
    temp_sum['run'].append('fit')
    temp_sum['log_F'].append(fit_terms['c71r18'][conc]['x'][-1])
    temp_sum['lPs'].append(fit_terms['c71r18'][conc]['x'][:-1])
    temp_sum['weights'].append(np.ones_like(all_vects['c71r18']['ltminus'][conc]))
    
    for run, n in enumerate(resamplings_c71[conc]['log_F']):
        temp_sum['run'].append('samp_%02d'%(run))
        temp_sum['log_F'].append(resamplings_c71[conc]['log_F'][run])
        temp_sum['lPs'].append(resamplings_c71[conc]['lp_vals'][run])
        temp_sum['weights'].append(resamplings_c71[conc]['sets_used'][run])
        
    plist=[]
    for n, cons in enumerate(all_vects['c71r18']['locs'][conc]):
        temp_sum[cons+'_log_P']=[]
        temp_sum[cons+'_weight']=[]
        plist.append(cons+'_log_P')
        plist.append(cons+'_weight')
        
    for run, m in enumerate(temp_sum['weights']):
        count=0
        for n, cons in enumerate(all_vects['c71r18']['locs'][conc]):

            weight=temp_sum['weights'][run][n]
            temp_sum[cons+'_weight'].append(weight)
            if weight==0:
                temp_sum[cons+'_log_P'].append(np.nan)
                continue
            temp_sum[cons+'_log_P'].append(temp_sum['lPs'][run][count])
            count+=1
       
    #generate a summary of each concentration for c71     
    temp_sum_df=pd.DataFrame(temp_sum, index=temp_sum['run'], columns=['log_F']+plist)
    temp_sum_df.to_csv('../summaries/c71_resamp_'+conc+'uM.txt', sep='\t', index_label='run')
    temp_sum_df.to_excel(writer, 'c71_resamp_'+conc+'uM', index_label='run')
    
    
            
for conc in resamplings_oc:
    temp_sum={'run':[], 'log_F':[], 'lPs':[], 'weights':[]}
    
    temp_sum['run'].append('fit')
    temp_sum['log_F'].append(fit_terms['occlusion'][conc]['x'][-1])
    temp_sum['lPs'].append(fit_terms['occlusion'][conc]['x'][:-1])
    temp_sum['weights'].append(np.ones_like(all_vects['occlusion']['ltminus'][conc]))
    
    for run, n in enumerate(resamplings_oc[conc]['log_F']):
        temp_sum['run'].append('samp_%02d'%(run))
        temp_sum['log_F'].append(resamplings_oc[conc]['log_F'][run])
        temp_sum['lPs'].append(resamplings_oc[conc]['lp_vals'][run])
        temp_sum['weights'].append(resamplings_oc[conc]['sets_used'][run])
        
    plist=[]
    for n, cons in enumerate(all_vects['occlusion']['locs'][conc]):
        temp_sum[cons+'_log_P']=[]
        temp_sum[cons+'_weight']=[]
        plist.append(cons+'_log_P')
        plist.append(cons+'_weight')
        
    for run, m in enumerate(temp_sum['weights']):
        count=0
        
        for n, cons in enumerate(all_vects['occlusion']['locs'][conc]):
            weight=temp_sum['weights'][run][n]
            temp_sum[cons+'_weight'].append(weight)
            if weight==0:
                temp_sum[cons+'_log_P'].append(np.nan)
                continue
            temp_sum[cons+'_log_P'].append(temp_sum['lPs'][run][count])
            count+=1
            
    #generate a summary of each concentration for occlusion
    temp_sum_df=pd.DataFrame(temp_sum, index=temp_sum['run'], columns=['log_F']+plist)
    temp_sum_df.to_csv('../summaries/occlusion_resamp_'+conc+'uM.txt', sep='\t', index_label='run')
    temp_sum_df.to_excel(writer, 'occlusion_resamp_'+conc+'uM', index_label='run')
    


#beta fitting

curves={'c60':'c60r18.10', 'c61':'c61r18', 'c62':'c62r18.10', 'c63':'c63r18.10', 'c64':'c64r18.10', 
        'c65':'c65r18.10', 'c71':'c71r18', 'c72':'c72r18.10',
        'c81':'c81r18.10', 'c82':'c82r18.10', 'c-':'c-', 'occlusion':'occlusion', 
        'c61r17':'c61r17.cons', 'c41':'gal.c41'}

#starting all_vects over fresh, as cAMP dilution data had a different structure
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
        
    
#make a new lib of all c61 constructs, r18 and r17 combined
conjoined={'lt+':[], 'lt-':[], 'locs':[]}

for lib in ['c61r18', 'c61r17.cons']:
    for category in ['lt+', 'lt-', 'locs']:
        for cons in all_vects[lib][category]:
            conjoined[category].append(cons)
            
all_vects['c61']=conjoined

#of the libs, make a list to fit, excluding the c- and occlusion ones
fit_set=['c60r18.10', 'c61r18', 'c62r18.10', 'c63r18.10', 'c64r18.10', 'c65r18.10', 'c71r18', 'c72r18.10', 'c81r18.10', 'c82r18.10', 'c61', 'gal.c41']
        
#load in the beta fit for each construct done with all data points
fits=pickle.load(open('../intermediates/beta_weighted_all_thetas_fit.pkl', 'rb'))

constructs=['c60r18.10', 'c61r18', 'c62r18.10', 'c63r18.10', 'c64r18.10', 'c65r18.10', 
            'c71r18', 'c72r18.10', 'c81r18.10', 'c82r18.10', 'c61', 'gal.c41']

spac={'c60r18.10':'c60', 'c61r18':'c61r18', 'c62r18.10':'c62', 'c63r18.10':'c63', 'c64r18.10':'c64', 'c65r18.10':'c65',
      'c71r18':'c71', 'c72r18.10':'c72', 'c81r18.10':'c81', 'c82r18.10':'c82', 'c61':'c61', 'gal.c41':'c41'}

resamplings={}

for cons in constructs:
    resamplings[cons]=pd.read_pickle('../intermediates/beta_weighted_'+cons+'_thetas_100.pkl')

for lib in resamplings:
    temp_sum={'run':[], 'log_tbg':[], 'log_alphap':[], 'log_betap':[], 'lPs':[], 'weights':[]}
    
    temp_sum['run'].append('fit')
    temp_sum['log_tbg'].append(fits[lib]['x'][-3])
    temp_sum['log_alphap'].append(fits[lib]['x'][-2])
    temp_sum['log_betap'].append(fits[lib]['x'][-1])
    temp_sum['lPs'].append(fits[lib]['x'][:-3])
    temp_sum['weights'].append(np.ones_like(all_vects[lib]['lt-']))
    
    for run, n in enumerate(resamplings[lib]['log_tbackground']):
        temp_sum['run'].append('samp_%02d'%(run))
        temp_sum['log_tbg'].append(resamplings[lib]['log_tbackground'][run])
        temp_sum['log_alphap'].append(resamplings[lib]['log_alphap'][run])
        temp_sum['log_betap'].append(resamplings[lib]['log_betap'][run])
        temp_sum['lPs'].append(resamplings[lib]['P_vals'][run])
        temp_sum['weights'].append(resamplings[lib]['sets_used'][run])
        
    plist=[]
    for n, cons in enumerate(all_vects[lib]['locs']):
        temp_sum[cons+'_log_P']=[]
        temp_sum[cons+'_weight']=[]
        plist.append(cons+'_log_P')
        plist.append(cons+'_weight')
        
    for run, m in enumerate(temp_sum['weights']):
        count=0
        
        for n, cons in enumerate(all_vects[lib]['locs']):
            weight=temp_sum['weights'][run][n]
            temp_sum[cons+'_weight'].append(weight)
            if weight==0:
                temp_sum[cons+'_log_P'].append(np.nan)
                continue
            temp_sum[cons+'_log_P'].append(temp_sum['lPs'][run][count])
            count+=1
            
    temp_sum_df=pd.DataFrame(temp_sum, index=temp_sum['run'], columns=['log_tbg', 'log_alphap', 'log_betap']+plist)
    temp_sum_df.to_csv('../summaries/'+spac[lib]+'_beta_resamp.txt', sep='\t', index_label='run')
    temp_sum_df.to_excel(writer, spac[lib]+'_beta_resamp', index_label='run')
    
    
            

writer.save()