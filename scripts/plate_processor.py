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

temp_panel.to_pickle('../data/plate_panel.pkl')