plate_name = 'plate080u'

induced = 2

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2017_12_04_lacstar_c61r18l35lib_halfconc_plate1_uninduced_3h_OD600.txt"
raw_miller = "2017_12_04_lacstar_c61r18l35lib_halfconc_plate1_uninduced_miller_10h_420.txt"
raw_optics = "2017_12_04_lacstar_c61r18l35lib_halfconc_plate1_uninduced_miller_10h_550.txt"

num_rows = 8
num_cols = 12
date = '17.12.04'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b1B7',
	3:'b1B8',
	4:'b1B9',
	5:'b1C1',
	6:'b1C2',
	7:'b1C3',
	8:'b1C4',
	9:'b1C5',
	10:'b1C6',
	11:'b1C7',
	12:'b1C9',
	13:'b1D1',
	14:'b1D2',
	15:'b1D3',
	16:'b1D4',
	17:'b1D5',
	18:'b1D6',
	19:'b1D7',
	20:'b1D8',
	21:'b1D9',
	22:'b1E1',
	23:'b1E2',
	24:'b1F1'
	}
	
num_strains = max(cultures_dict.keys())
	
culture_locations_string = '''
1	0	5	0	9	0	13	0	17	0	21	0
1	3	5	7	9	11	13	15	17	19	21	23
1	3	5	7	9	11	13	15	17	19	21	23
0	3	0	7	0	11	0	15	0	19	0	23
2	0	6	0	10	0	14	0	18	0	22	0
2	4	6	8	10	12	14	16	18	20	22	24
2	4	6	8	10	12	14	16	18	20	22	24
0	4	0	8	0	12	0	16	0	20	0	24
'''

camp_concentrations_string_u = '''
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
'''

camp_concentrations_string_i = '''
250	250	250	250	250	250	250	250	250	250	250	250
250	250	250	250	250	250	250	250	250	250	250	250
250	250	250	250	250	250	250	250	250	250	250	250
250	250	250	250	250	250	250	250	250	250	250	250
250	250	250	250	250	250	250	250	250	250	250	250
250	250	250	250	250	250	250	250	250	250	250	250
250	250	250	250	250	250	250	250	250	250	250	250
250	250	250	250	250	250	250	250	250	250	250	250
'''

camp_concentrations_string_h = '''
125	125	125	125	125	125	125	125	125	125	125	125
125	125	125	125	125	125	125	125	125	125	125	125
125	125	125	125	125	125	125	125	125	125	125	125
125	125	125	125	125	125	125	125	125	125	125	125
125	125	125	125	125	125	125	125	125	125	125	125
125	125	125	125	125	125	125	125	125	125	125	125
125	125	125	125	125	125	125	125	125	125	125	125
125	125	125	125	125	125	125	125	125	125	125	125
'''

if induced==1: camp_concentrations_string = camp_concentrations_string_i
elif induced==2: camp_concentrations_string = camp_concentrations_string_u
elif induced==3: camp_concentrations_string = camp_concentrations_string_h

culture_volumes_string = '''
50	50	50	50	50	50	50	50	50	50	50	50
50	50	50	50	50	50	50	50	50	50	50	50
50	50	50	50	50	50	50	50	50	50	50	50
50	50	50	50	50	50	50	50	50	50	50	50
50	50	50	50	50	50	50	50	50	50	50	50
50	50	50	50	50	50	50	50	50	50	50	50
50	50	50	50	50	50	50	50	50	50	50	50
50	50	50	50	50	50	50	50	50	50	50	50
'''

