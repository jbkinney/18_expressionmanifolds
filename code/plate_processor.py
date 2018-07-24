#!/usr/bin/env python3
import os
import pylab as py
import pandas as pd
import matplotlib
import numpy as np
import scipy.optimize as opt
import glob

import pdb

class Reaction: pass;


exec(open(os.pardir+'/code/library_clusters.py').read(), globals())

#code is in a folder at the same level as data
datapath=os.pardir+'/data/'

#platelist is a list of all the spec sheets in the folder
platelist = glob.glob('../data/plate_reader/metadata/specs_*.py')

shared_folder = os.pardir+'/intermediate/'

#raw_data is the path for the folder containing the raw data for each plate
raw_data=datapath+'plate_reader/raw/'

py.close('all')

cons_names=pd.read_excel(datapath+'glycerol_stocks.xlsx', index_col=0)

#if a timepoint for the well has higher A600 than threshold, it will get excluded
od_threshold = 2.0

row_names=['A','B','C','D','E','F','G','H']
col_names=['1','2','3','4','5','6','7','8','9','10','11','12']

#plate_processor takes in the spec file for a plate and returns a dictionary of processed reaction data for the wells
def plate_processor(specs):
    plate = specs.split('/')[-1]
    #specs=plate_metadata+plate
    plate_rx=[]
    plate_rxn_dic={}
    
    #execfile(specs, globals())
    print('Processing plate %s...'%specs)
    exec(open(specs).read(), globals())
    
    #parse culture locations
    cultures_mat = py.zeros([num_rows, num_cols])
    lines=[l for l in culture_locations_string.split('\n') if len(l.split())==num_cols]
    for r in range(num_rows):
        for c in range(num_cols):
            cultures_mat[r,c] = lines[r].split()[c]

    # Parse culture volumes
    culture_volumes = py.zeros([num_rows, num_cols])
    lines = [l for l in culture_volumes_string.split('\n') if len(l.split())==num_cols]
    for r in range(num_rows):
        for c in range(num_cols):
            culture_volumes[r,c] = lines[r].split()[c]

    # Parse cAMP concentrations
    camp_concs = py.zeros([num_rows, num_cols])
    lines = [l for l in camp_concentrations_string.split('\n') if len(l.split())==num_cols]
    for r in range(num_rows):
        for c in range(num_cols):
                camp_concs[r,c] = lines[r].split()[c]
                
                
    #read in the A420 data
    rxn_dataframe=pd.read_csv(raw_data+miller_csv, header=1)

    #read in the A550 data
    optics_dataframe = pd.read_csv(raw_data+optics_csv, header=1)
    
    rtim=[]
    otim=[]

    #processes the hh:mm:ss time into float minutes and appends them to rtim or otim
    for n in range(len(rxn_dataframe['Time'])):
        rtemp=rxn_dataframe['Time'][n].split(':')
        rtemp=list(map(float,rtemp))
    
        otemp=optics_dataframe['Time'][n].split(':')
        otemp=list(map(float,otemp))
    
        rtim.append(rtemp[0]*60.0+rtemp[1]+rtemp[2]/60.0)
        otim.append(otemp[0]*60.0+otemp[1]+otemp[2]/60.0)
    
    
    rxn_dataframe=rxn_dataframe.assign(Times = rtim)
    optics_dataframe=optics_dataframe.assign(Times=otim)

    ts = rxn_dataframe['Times']


    for r in row_names:
        for c in col_names:
            #this does the float mapping for all the OD values
            rxn_dataframe[r+c]=list(map(float,rxn_dataframe[r+c]))
            optics_dataframe[r+c]=list(map(float,optics_dataframe[r+c]))
        

    with open(raw_data+growth_file, 'U') as f:
        sdata = [l.split()[1:(num_cols+1)] for l in f.readlines()[2:(2+num_rows)]]
        growth_data = np.array([[float(c) for c in row] for row in sdata])
        growth_data=growth_data.astype(float)
    
    #grabs all wells with just culture media to determine bg A600 value
    mask = (cultures_mat == 0)
    growth_base = py.median(growth_data[mask])
    #subtracts bg A600 value from the plate's A600 data
    gods = growth_data - growth_base
    #sets any well with negative corrected A600 to epsilon
    gods[gods <= 0] = np.finfo(float).eps
    
    
    
    #builds all the Reactions, gives their well location and all data and conditions for the well

    for r in range(num_rows):
        for c in range(num_cols):
        
            rxn = Reaction()
            
            rxn.r = r
            rxn.c = c

            #the Reaction is given the corrected normed values in .culture_od600
            rxn.culture_od600 = gods[r, c]
            
            #the Reaction is given the culture ID in .culture, the name in .name
            rxn.culture = cultures_mat[r, c]
            rxn.loc = cultures_dict[rxn.culture]
            
            #if the well was a blank or contained RDM, set the name to blank and good sequence status to FALSE
            if cultures_dict[rxn.culture]=='RDM' or cultures_dict[rxn.culture]=='blank':
                rxn.name='blank'
                rxn.seq_stat='FALSE'
                
            else:
                rxn.name = cons_names['short_name'][rxn.loc]
                rxn.seq_stat = cons_names['valid_clone'][rxn.loc]
        
            rxn.well = row_names[r]+col_names[c]
            
            #the Reaction is given the cAMP concentration in .camp
            rxn.camp = camp_concs[r, c]

            #the Reaction is given the culture volume in .culture_vol
            rxn.culture_vol = culture_volumes[r, c]
            
            #puts all the Miller measurements for the timecourse in this well in .od420s
            rxn.od420s = rxn_dataframe[row_names[r]+col_names[c]]
            rxn.od550s = optics_dataframe[row_names[r]+col_names[c]]
            
            
            #this should put a vector of bools for whether or not a well is within the timepoints of interest
            #and also both the miller measurements and the condensation measurements are below threshold into .ok_indices
            rxn.ok_indices = (ts <= stop_time)&(ts >= start_time)&(rxn.od420s < od_threshold)&(rxn.od550s < od_threshold)&(rxn.od420s - rxn.od550s<0.5)
            
            
            #the signal is set to the difference between the miller and growth values minus the smallest difference
            rxn.signal = rxn.od420s - rxn.od550s - min(rxn.od420s - rxn.od550s)
            
            #fits first order polynomial to signal vs timeplot for the ok points
            ##checks that there is a minimum number of good values
            if sum(rxn.ok_indices)<2:
                rxn.coeffs = []
                rxn.fit = []
                rxn.rate = np.nan
                
            else:
                rxn.coeffs = np.polyfit(ts[rxn.ok_indices], rxn.signal[rxn.ok_indices], 1)
            
                #the .fit vector is the signal values predicted from the polyfit for every timepoint
                rxn.fit = py.polyval(rxn.coeffs,ts)
            
                #.rate gets the slope of the polyfit
                rxn.rate = rxn.coeffs[-2]
                
            #this puts the plate name in .extension
            rxn.extension = plate
            
            rxn.millers = pow(10,3.6)*(rxn.rate) / (rxn.culture_od600 * rxn.culture_vol)
            
            #adds this Reaction to the array of rxns initialized earlier
            plate_rx.append(rxn)
            plate_rxn_dic[rxn.well]=rxn
    
    return plate_rx, plate_rxn_dic, cultures_dict
    
plate_rxns={}
strain_locs={}
rxn_dics={}
strain_names={}

#for all spec sheets in the metadata folder, process and pull the reaction data
for val in platelist:
    plate_rxns[val], rxn_dics[val], strain_locs[val] = plate_processor(val)

       
#dictionary of reactions by plate
temp_df={}

for plate in platelist:
    temp=[]
    
    #for each reaction, grab the pertinent info and store it in temp
    for rxn in plate_rxns[plate]:
        temp.append([rxn.well, rxn.name, rxn.millers, rxn.culture_od600, rxn.rate, rxn.camp, rxn.culture_vol, rxn.loc, rxn.seq_stat])
    
    temp=np.array(temp)

    #store the plate's reaction data in temp_df
    temp_df[plate]=temp
    
#make a panel of the reaction data for all plates
temp_panel=pd.Panel(temp_df, minor_axis=["well","name","Miller meas","OD600","rxn rate","cAMP","volume","glycerol loc","OK seq?"])

#This takes the panel of data and summarizes it
raw_vals={}
med_vals={}
raw_screened={}

#First, go through the panel of data and grab all the raw values, by glycerol stock location
for plate in temp_panel:
    for val in range(len(temp_panel[plate])):
        nam=temp_panel[plate]["glycerol loc"][val]
        
        if nam=='RDM':
            continue
            
        if nam=='blank': continue
        
        if cons_names['valid_clone'][nam]=='FALSE': continue
        if cons_names['valid_clone'][nam]=='': continue
            
        if nam not in raw_vals:
            raw_vals[nam]={'bas':[], 'ind':[]}
        
        if temp_panel[plate]['cAMP'][val]=='0.0':
            raw_vals[nam]['bas'].append(float(temp_panel[plate]['Miller meas'][val]))
            
        elif temp_panel[plate]['cAMP'][val]=='250.0':
            raw_vals[nam]['ind'].append(float(temp_panel[plate]['Miller meas'][val]))
            
#goes through raw data for each glycerol stock location and takes the median
for sample in raw_vals:
    
    bas_temp=np.nan
    ind_temp=np.nan
    
    raw_bas_screened=[]
    raw_ind_screened=[]
    
    #grab all positive raw vals for this construct
    for val in raw_vals[sample]['bas']:
        if val>0:
            raw_bas_screened.append(val)
            
    for val in raw_vals[sample]['ind']:
        if val>0:
            raw_ind_screened.append(val)
    
    #if you have at least 3 positive raw vals, find the median
    if len(raw_bas_screened)>2:
        bas_temp=np.median(raw_bas_screened)
        
    if len(raw_ind_screened)>2:
        ind_temp=np.median(raw_ind_screened)
    
    med_vals[sample]={'bas':bas_temp, 'ind':ind_temp}
    raw_screened[sample]={'bas':raw_bas_screened, 'ind':raw_ind_screened}
    

summary={'name':[], 'location':[], 'median_u':[], 'median_i':[], 'num_u':[], 'num_i':[], 'all_u':[], 'all_i':[]}

for cons in med_vals:
    if cons=='blank' or cons=='RDM': continue
    summary['location'].append(cons)
    summary['name'].append(cons_names['short_name'][cons])
    summary['median_u'].append(med_vals[cons]['bas'])
    summary['median_i'].append(med_vals[cons]['ind'])
    summary['all_u'].append(raw_screened[cons]['bas'])
    summary['all_i'].append(raw_screened[cons]['ind'])
    summary['num_u'].append(len(raw_screened[cons]['bas']))
    summary['num_i'].append(len(raw_screened[cons]['ind']))
    
#go one folder up then down into the intermediate folder to save the measurements as measurements.txt
shared_folder = os.pardir+'/intermediate/'
namfile='measurements.txt'
all_sum = pd.DataFrame(summary, index=summary['name'], columns=[ 'location', 'median_u', 'median_i', 'num_u', 'num_i', 'all_u', 'all_i'])
all_sum.to_csv(shared_folder+namfile, sep='\t', na_rep=np.nan)

#this puts the row labels back to rights, as 'name' gets clipped out
with open(shared_folder+namfile, 'r+') as named_summary:
    temp=named_summary.read()

    with open(shared_folder+namfile, 'w+') as named_summary:
        named_summary.write('name'+temp)
        
        



lib_groups=['c61r18','c71r18','occlusion',
            '61c-r18.35','wtc60r18','wtc61r18.10',
            'wtc71r18.10','c60r18.10','c62r18.10',
            'c81r18.10','c63r18.10','c64r18.10',
            'c65r18.10','c66r18.10','c76r18.10',
            'c82r18.10','c72r18.10', 'gal', 'c40', 'c41'
           ]

group_summary={'name':[], 'loc':[], 'median_basal':[], 'median_induced':[], 'inliers':[]}
    
for lib in lib_groups:
    for construct in library_groups[lib]['all']:
        if construct not in med_vals: continue
        
        group_summary['loc'].append(construct)
    
        if '41.gal' not in cons_names['short_name'][construct]:
            group_summary['name'].append(cons_names['short_name'][construct])
        else:
            group_summary['name'].append(cons_names['short_name'][construct][4:])
    
        group_summary['median_basal'].append(med_vals[construct]['bas'])
    
        group_summary['median_induced'].append(med_vals[construct]['ind'])
    
        group_summary['inliers'].append(construct not in library_groups[lib]['outliers'])
        
for lib in rgap_cons:
    for construct in rgap_cons[lib]:
        if construct not in med_vals: continue
        
        group_summary['loc'].append(construct)
    
        group_summary['name'].append(cons_names['short_name'][construct])
    
        group_summary['median_basal'].append(med_vals[construct]['bas'])
    
        group_summary['median_induced'].append(med_vals[construct]['ind'])
    
        group_summary['inliers'].append(True)
        
group_nam='clonal_measurements.txt'
group_sum = pd.DataFrame(group_summary, index=group_summary['name'], columns=[ 'loc', 'median_basal', 'median_induced', 'inliers'])
group_sum.to_csv(shared_folder+group_nam, sep='\t', na_rep=np.nan)

with open(shared_folder+group_nam, 'r+') as named_summary:
    temp=named_summary.read()

    with open(shared_folder+group_nam, 'w+') as named_summary:
        named_summary.write('name'+temp)
        
def comp_ltaui(lrn,lthetas): 
    """
    Calculates and returns the log of the induced transcription.
    
    Parameters
    ----------
    lrn : float
        The log of the e^-r for this promoter.
        
    lthetas : float array
        lthetas[0] is the log of the maximum transcription
        lthetas[1] is the log of the background transcription
        lthetas[2] is the log of the eta (eta = (1+e^-c-i)/(1+e^-c))
    
    Returns
    -------
    float
        The log of the calculated induced transcription,
        induced transc = tm*(eta*e^-r/(1+eta*e^-r))
    """

    rn=np.exp(lrn)
    tm=np.exp(lthetas[0])
    tbg=np.exp(lthetas[1])
    et=np.exp(lthetas[2])

    
    #partition fn
    Z=1+rn*et
    
    #total transc
    ti=(tm*rn*et/(1+rn*et))+tbg
    
    return np.log(ti)

#calcs log of tau basal from same givens as above
def comp_ltaub(lrn,lthetas): 
    """
    Calculates and returns the log of the basal transcription.
    
    Parameters
    ----------
    lrn : float
        The log of the e^-r for this promoter.
        
    lthetas : float array
        lthetas[0] is the log of the maximum transcription
        lthetas[1] is the log of the background transcription
        lthetas[2] is the log of the eta (eta = (1+e^-c-i)/(1+e^-c)) (optional, not used)
    
    Returns
    -------
    float
        The log of the calculated basal transcription,
        basal transc = tm*(e^-r/(1+e^-r))
    """

    rn=np.exp(lrn)
    tm=np.exp(lthetas[0])
    tbg=np.exp(lthetas[1])
    et=np.exp(lthetas[2])
    
    Z=1+rn
    
    #total transc
    tb=(tm*rn/(1+rn))+tbg
    
    return np.log(tb)
    
    
#objective function, given the measured values for log basal transc and log induced transc,
#and a vector for the log of all variables fitted
#returns the squared difference in calculated log induced transc from each measured log induc transc we've got
#plus calculated log basal transc from each measured log basal transc
def comp_error(x, lbasal, linduced):
    
    err=[]
    
    for y in range(len(lbasal)):
        
        #subtract the calculated log tau basal and log tau induced from measured
        errb=lbasal[y]-comp_ltaub(x[y],x[-3:]) 
        erri=linduced[y]-comp_ltaui(x[y],x[-3:])
        
        #add squares of those differences
        err.append((errb**2)+(erri**2))
        
    return sum(err)

    
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
    
    
def comp_error_oc(x, lbasal, linduced):
    
    err=[]
    
    for y in range(len(lbasal)):
        
        #subtract the calculated log tau basal and log tau induced from measured
        errb=lbasal[y]-comp_ltaub_oc(x[y],x[-3:]) 
        erri=linduced[y]-comp_ltaui_oc(x[y],x[-3:])
        
        #add squares of those differences
        err.append((errb**2)+(erri**2))
        
    return sum(err)
    
#objective function for each library is defined as the comp_error run on whatever is handed to the objective function,
#plus the measured basal and induced median values
#we then optimize the error
def objective_func_creator(lib):
    bas=all_vects[lib]['bas']
    ind=all_vects[lib]['ind']
    lammy=lambda m: comp_error(m,bas,ind)
    return lammy


all_vects={}

for val in lib_groups:
    bas_temp=[]
    ind_temp=[]
    
    for cons in library_groups[val]['all']:
        if cons in library_groups[val]['outliers']: continue 
        if cons not in med_vals: continue
        
        cur_bas=np.log(med_vals[cons]['bas'])
        cur_ind=np.log(med_vals[cons]['ind'])
        
        if np.isnan(cur_bas) or np.isnan(cur_ind):
            continue
        
        bas_temp.append(cur_bas)
        ind_temp.append(cur_ind)
        
            
    all_vects[val]={'bas':bas_temp,'ind':ind_temp}
    
    
standard_libs=[]
occlusion_libs=[]

for lib in all_vects.keys():
    if 'oc' in lib:
        occlusion_libs.append(lib)
    else:
        standard_libs.append(lib)
        
for lib in standard_libs:
    all_vects[lib]['objective_func']=objective_func_creator(lib)
    
    
all_vects['occlusion']['objective_func']=lambda m: comp_error_oc(m,all_vects['occlusion']['bas'],all_vects['occlusion']['ind'])
        

fit_set=standard_libs

#initializing values are set here
tm=np.log(200)
tbg=np.log(0.0001)
et=5

init=[tm,tbg,et]

fits={}
mfts={}

for cons in fit_set:
    temp=[]
    for x in all_vects[cons]['bas']:
        rtemp=(np.exp(x)-np.exp(tbg))/(np.exp(tm)-np.exp(x)+np.exp(tbg))
        temp.append(np.log(rtemp))
        
    all_vects[cons]['initial']=temp
    all_vects[cons]['full_init']=temp+init
    
for cons in fit_set:
    fits[cons]=opt.minimize(all_vects[cons]['objective_func'], all_vects[cons]['full_init'], method='L-BFGS-B')
    mfts[cons]=fits[cons]['x'][-3:]
    



ec=1
init_oc=[tm,tbg,ec]
fit_oc={}
mfts_oc={}

for constr in occlusion_libs:
    temp=[]
    for x in all_vects[constr]['bas']:
        rtemp=(np.exp(x)-np.exp(tbg))/(np.exp(tm)-np.exp(x)+np.exp(tbg))
        temp.append(np.log(rtemp))
        
    all_vects[constr]['initial']=temp
    all_vects[constr]['full_init']=temp+init_oc
    
    fit_oc[constr]=opt.minimize(all_vects[constr]['objective_func'], all_vects[constr]['full_init'], method='L-BFGS-B')
    mfts_oc[constr]=fit_oc[constr]['x'][-3:]


resamplings={}

for cons in lib_groups:
    resamplings[cons]=pd.read_pickle(shared_folder+'resamplings/'+cons+'_thetas_100.pkl')
        

#generates resampled_params_for_occlusion.txt
oc_sum={'run':[], 'log_tsat':[], 'log_tbg':[], 'log_C':[]}

oc_sum['run'].append('best')
oc_sum['log_tsat'].append(mfts_oc['occlusion'][0])
oc_sum['log_tbg'].append(mfts_oc['occlusion'][1])
oc_sum['log_C'].append(mfts_oc['occlusion'][2])

for run, n in enumerate(resamplings['occlusion']['log_tmax']):
    oc_sum['run'].append('samp_%02d'%(run))
    oc_sum['log_tsat'].append(resamplings['occlusion']['log_tmax'][run])
    oc_sum['log_tbg'].append(resamplings['occlusion']['log_tbackground'][run])
    oc_sum['log_C'].append(resamplings['occlusion']['-c'][run])
    
oc_nam='resampled_params_for_occlusion.txt'
oc_summary = pd.DataFrame(oc_sum, index=oc_sum['run'], columns=['log_tsat', 'log_tbg', 'log_C'])
oc_summary.to_csv(shared_folder+oc_nam, sep='\t', na_rep=np.nan)

with open(shared_folder+oc_nam, 'r+') as named_summary:
    temp=named_summary.read()

    with open(shared_folder+oc_nam, 'w+') as named_summary:
        named_summary.write('run'+temp)

#generates resampled_params_for_c61.txt        
c61_sum={'run':[], 'log_tsat':[], 'log_tbg':[], 'log_eta':[]}

c61_sum['run'].append('best')
c61_sum['log_tsat'].append(mfts['c61r18'][0])
c61_sum['log_tbg'].append(mfts['c61r18'][1])
c61_sum['log_eta'].append(mfts['c61r18'][2])

for run, n in enumerate(resamplings['c61r18']['log_tmax']):
    c61_sum['run'].append('samp_%02d'%(run))
    c61_sum['log_tsat'].append(resamplings['c61r18']['log_tmax'][run])
    c61_sum['log_tbg'].append(resamplings['c61r18']['log_tbackground'][run])
    c61_sum['log_eta'].append(resamplings['c61r18']['log_interactionweight'][run])
    
c61_nam='resampled_params_for_c61.txt'
c61_summary = pd.DataFrame(c61_sum, index=c61_sum['run'], columns=['log_tsat', 'log_tbg', 'log_eta'])
c61_summary.to_csv(shared_folder+c61_nam, sep='\t', na_rep=np.nan)

with open(shared_folder+c61_nam, 'r+') as named_summary:
    temp=named_summary.read()

    with open(shared_folder+c61_nam, 'w+') as named_summary:
        named_summary.write('run'+temp)
    
#generates resampled_params_for_c71.txt    
c71_sum={'run':[], 'log_tsat':[], 'log_tbg':[], 'log_eta':[]}

c71_sum['run'].append('best')
c71_sum['log_tsat'].append(mfts['c71r18'][0])
c71_sum['log_tbg'].append(mfts['c71r18'][1])
c71_sum['log_eta'].append(mfts['c71r18'][2])

for run, n in enumerate(resamplings['c71r18']['log_tmax']):
    c71_sum['run'].append('samp_%02d'%(run))
    c71_sum['log_tsat'].append(resamplings['c71r18']['log_tmax'][run])
    c71_sum['log_tbg'].append(resamplings['c71r18']['log_tbackground'][run])
    c71_sum['log_eta'].append(resamplings['c71r18']['log_interactionweight'][run])
    
c71_nam='resampled_params_for_c71.txt'
c71_summary = pd.DataFrame(c71_sum, index=c61_sum['run'], columns=['log_tsat', 'log_tbg', 'log_eta'])
c71_summary.to_csv(shared_folder+c71_nam, sep='\t', na_rep=np.nan)

with open(shared_folder+c71_nam, 'r+') as named_summary:
    temp=named_summary.read()

    with open(shared_folder+c71_nam, 'w+') as named_summary:
        named_summary.write('run'+temp)
        
        
error_set={'c61r18':'61.5', 'c60r18.10':'60.5', 'c62r18.10':'62.5', 'c63r18.10':'63.5', 'c64r18.10':'64.5', 'c65r18.10':'65.5', 'c66r18.10':'66.5', 'c71r18':'71.5', 'c72r18.10':'72.5', 'c76r18.10':'76.5', 'c81r18.10':'81.5', 'c82r18.10':'82.5'}

fit_C=np.exp(mfts_oc['occlusion'][2])

spacing_alphas={}
    
for con in error_set:
    resamp_alphas=[]
    fit_eta = np.exp(mfts[con][2])
    fit_alpha = ((fit_eta*(1+fit_C))-1)/fit_C
    
    for val in resamplings[con]['log_interactionweight']:
        alpha=((np.exp(val)*(1+fit_C))-1)/fit_C
        resamp_alphas.append(alpha)
    
    spacing_alphas[error_set[con]]={'alpha':fit_alpha, 
                                    'alpha_16':np.percentile(resamp_alphas, 16),
                                    'alpha_50':np.percentile(resamp_alphas, 50),
                                    'alpha_84':np.percentile(resamp_alphas, 84)
                                    }
                                    
spacing_tsats={}
    
for con in error_set:
    resamp_tsats=[]
    
    for val in resamplings[con]['log_tmax']:
        resamp_tsats.append(np.exp(val))
    
    
    spacing_tsats[error_set[con]]={'tsat':np.exp(mfts[con][0]), 
                                    'tsat_16':np.percentile(resamp_tsats, 16),
                                    'tsat_50':np.percentile(resamp_tsats, 50),
                                    'tsat_84':np.percentile(resamp_tsats, 84)
                                    }
                                    
spacing_tbg={}
    
for con in error_set:
    resamp_tbg=[]
    
    for val in resamplings[con]['log_tbackground']:
        resamp_tbg.append(np.exp(val))
    
    
    spacing_tbg[error_set[con]]={'tbg':np.exp(mfts[con][1]), 
                                    'tbg_16':np.percentile(resamp_tbg, 16),
                                    'tbg_50':np.percentile(resamp_tbg, 50),
                                    'tbg_84':np.percentile(resamp_tbg, 84)
                                    }
                                    
spacings=spacing_tbg.keys()
spacings=sorted(spacings)

dist_summary={'dist':[],
              'alpha':[],'alpha_16':[],'alpha_50':[],'alpha_84':[],
              'tsat':[],'tsat_16':[],'tsat_50':[],'tsat_84':[],
              'tbg':[],'tbg_16':[],'tbg_50':[],'tbg_84':[]
             }

for dist in spacings:
    dist_summary['dist'].append(dist)
    dist_summary['alpha'].append(spacing_alphas[dist]['alpha'])
    dist_summary['alpha_16'].append(spacing_alphas[dist]['alpha_16'])
    dist_summary['alpha_50'].append(spacing_alphas[dist]['alpha_50'])
    dist_summary['alpha_84'].append(spacing_alphas[dist]['alpha_84'])
    dist_summary['tsat'].append(spacing_tsats[dist]['tsat'])
    dist_summary['tsat_16'].append(spacing_tsats[dist]['tsat_16'])
    dist_summary['tsat_50'].append(spacing_tsats[dist]['tsat_50'])
    dist_summary['tsat_84'].append(spacing_tsats[dist]['tsat_84'])
    dist_summary['tbg'].append(spacing_tbg[dist]['tbg'])
    dist_summary['tbg_16'].append(spacing_tbg[dist]['tbg_16'])
    dist_summary['tbg_50'].append(spacing_tbg[dist]['tbg_50'])
    dist_summary['tbg_84'].append(spacing_tbg[dist]['tbg_84'])
    
dist_summary_df = pd.DataFrame(dist_summary, index=dist_summary['dist'], 
                               columns=['alpha','alpha_16','alpha_50','alpha_84',
                                        'tsat','tsat_16','tsat_50','tsat_84',
                                        'tbg','tbg_16','tbg_50','tbg_84'])
dist_summary_df.to_csv(shared_folder+'params_versus_distance.txt', sep='\t')

with open(shared_folder+'params_versus_distance.txt', 'r+') as named_summary:
    temp=named_summary.read()

    with open(shared_folder+'params_versus_distance.txt', 'w+') as named_summary:
        named_summary.write('distance'+temp)