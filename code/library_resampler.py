import os
import pylab as py
import pandas as pd
import matplotlib
import numpy as np
import scipy.optimize as opt
from random import randint
from random import uniform


execfile(os.pardir+'/code/library_clusters.py')

resamp_folder = os.pardir+'/data/resamplings'

class Reaction: pass;

#code is in a folder at the same level as data
datapath=os.pardir+'/data/'

#plate_metadata is the path for the folder containing the spec sheets for each plate
plate_metadata=datapath+'plate_reader/metadata/'

#platelist is a list of all the spec sheets in the folder
platelist=os.listdir(plate_metadata)

#raw_data is the path for the folder containing the raw data for each plate
raw_data=datapath+'plate_reader/raw/'

py.close('all')

cons_names=pd.read_excel(datapath+'glycerol_stocks.xlsx', index_col=0)

#if a timepoint for the well has higher A600 than threshold, it will get excluded
od_threshold = 2.0

row_names=['A','B','C','D','E','F','G','H']
col_names=['1','2','3','4','5','6','7','8','9','10','11','12']

#plate_processor takes in the spec file for a plate and returns a dictionary of processed reaction data for the wells
def plate_processor(plate):
    specs=plate_metadata+plate
    plate_rx=[]
    plate_rxn_dic={}
    
    execfile(specs, globals())
    
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
        rtemp=map(float,rtemp)
    
        otemp=optics_dataframe['Time'][n].split(':')
        otemp=map(float,otemp)
    
        rtim.append(rtemp[0]*60.0+rtemp[1]+rtemp[2]/60.0)
        otim.append(otemp[0]*60.0+otemp[1]+otemp[2]/60.0)
    
    
    rxn_dataframe=rxn_dataframe.assign(Times = rtim)
    optics_dataframe=optics_dataframe.assign(Times=otim)

    ts = rxn_dataframe['Times']


    for r in row_names:
        for c in col_names:
            #this does the float mapping for all the OD values
            rxn_dataframe[r+c]=map(float,rxn_dataframe[r+c])
            optics_dataframe[r+c]=map(float,optics_dataframe[r+c])
        

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
    




lib_groups=['c61r18','c71r18','occlusion',
            '61c-r18.35','wtc60r18','wtc61r18.10',
            'wtc71r18.10','c60r18.10','c62r18.10',
            'c81r18.10','c63r18.10','c64r18.10',
            'c65r18.10','c66r18.10','c76r18.10',
            'c82r18.10','c72r18.10', 'gal', 'c40', 'c41'
           ]


        
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
    
    print "resampling "+construct
    lbasal=all_vects[construct]['bas']
    linduced=all_vects[construct]['ind']
    rns=all_vects[construct]['initial']
    
    lres_fits=[]
    lsets_used=[]
    lsets_used_temp=[]
    
    cntr=0
    
    while cntr<cycles:
        if cntr
    
        rtb=[]
        rti=[]
        rrn=[]
        sua=[]
        
        tm=np.log(uniform(np.exp(initializing[0]),10*np.exp(initializing[0])))
        tbg=np.log(uniform(0.1*np.exp(initializing[1]),np.exp(initializing[1])))
        et=np.log(uniform(0.5*np.exp(initializing[2]),5*np.exp(initializing[2])))
        init=[tm,tbg,et]
    
        for j in range(len(lbasal)):
            a=randint(0,len(lbasal)-1)
            rtb.append(lbasal[a])
            rti.append(linduced[a])
            rrn.append(rns[a])
            sua.append(a)
        
        comp_init=rrn+init
        lsets_used_temp.append(sua)
        
        lbasalr=rtb
        linducedr=rti
    
        # Specify data on which to compute objective function
        objective_func = lambda m: comp_error(m,lbasalr,linducedr)
        
        fits=opt.minimize(objective_func,comp_init, method='L-BFGS-B')
        lres_fits.append([fits['success'],fits['x'][-3:]])
        
        if fits['success']: cntr+=1
        
    lthetas=[]
    for ind, res in enumerate(lres_fits):
        #if the fit ended successfully, add the theta values
        if res[0]:
            lthetas.append(res[1])
            lsets_used.append(lsets_used_temp[ind])
        
    lthetas=np.array(lthetas)
    
    return lthetas, lsets_used

    
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
    
    print "resampling "+construct

    cntr=0
    
    lbasal=all_vects[construct]['bas']
    linduced=all_vects[construct]['ind']
    rns=all_vects[construct]['initial']
    
    lres_fits=[]
    lsets_used=[]
    lsets_used_temp=[]
    
    while cntr<cycles:
        rtb=[]
        rti=[]
        rrn=[]
        sua=[]
        
        tm=np.log(uniform(np.exp(initializing[0]),10*np.exp(initializing[0])))
        tbg=np.log(uniform(0.1*np.exp(initializing[1]),np.exp(initializing[1])))
        ec=np.log(uniform(0.5*np.exp(initializing[2]),5*np.exp(initializing[2])))
        init=[tm,tbg,ec]
    
        for j in range(len(lbasal)):
            a=randint(0,len(lbasal)-1)
            rtb.append(lbasal[a])
            rti.append(linduced[a])
            rrn.append(rns[a])
            sua.append(a)
        
        comp_init=rrn+init
        lsets_used_temp.append(sua)
        
        lbasalr=rtb
        linducedr=rti

    
        # Specify data on which to compute objective function
        objective_func = lambda m: comp_error_oc(m,lbasalr,linducedr)
        
        fits=opt.minimize(objective_func,comp_init, method='L-BFGS-B')
        lres_fits.append([fits['success'],fits['x'][-3:]])
        
        if fits['success']: cntr+=1
        
    lthetas=[]
    for ind, res in enumerate(lres_fits):
        #if the fit ended successfully, add the theta values
        if res[0]:
            lthetas.append(res[1])
            lsets_used.append(lsets_used_temp[ind])
        
    lthetas=np.array(lthetas)
    
    return lthetas, lsets_used
    
    

      
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
    
#objective function for each library is defined as the comp_error run on whatever is handed to the objective function,
#plus the measured basal and induced median values
#we then optimize the error
def objective_func_creator(lib):
    bas=all_vects[lib]['bas']
    ind=all_vects[lib]['ind']
    lammy=lambda m: comp_error(m,bas,ind)
    return lammy

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


resample_again=True
cycles=100
resamplings={}

if resample_again:
    for cons in fit_set:
        lthetas, lsets_used = boot_resampler(cons, init, cycles)
    
        resamp=pd.DataFrame(lthetas, columns=['log_tmax','log_tbackground','log_interactionweight'])
        resamp['sets_used']=[x for x in lsets_used]
    
        resamp.to_pickle(resamp_folder+"/"+cons+"_thetas_"+str(cycles)+".pkl")
    
        resamplings[cons]=resamp
    
    for construct in occlusion_libs:
        lthetas, lsets_used = boot_resampler_oc(construct, init_oc, cycles)
    
        resamp=pd.DataFrame(lthetas, columns=['log_tmax','log_tbackground','-c'])
        resamp['sets_used']=[x for x in lsets_used]
    
        resamp.to_pickle(resamp_folder+"/"+construct+"_thetas_"+str(cycles)+".pkl")
    
        resamplings[construct]=resamp
        
else:
    for cons in ['c61r18', 'c71r18', 'occlusion']:
        resamplings[cons]=pd.read_pickle(resamp_folder+'/'+cons+'_thetas_100.pkl')

