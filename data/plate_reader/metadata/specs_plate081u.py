plate_name = 'plate081u'

induced = 2

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2017_12_06_lacstar_c71r18l35lib_halfconc_plate1_uninduced_OD600_245h.txt"
raw_miller = "2017_12_06_lacstar_c71r18l35lib_halfconc_plate1_uninduced_miller_10h_420.txt"
raw_optics = "2017_12_06_lacstar_c71r18l35lib_halfconc_plate1_uninduced_miller_10h_550.txt"

num_rows = 8
num_cols = 12
date = '17.12.06'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b1C8',
	3:'b1I8',
	4:'b2A1',
	5:'b2A2',
	6:'b2A4',
	7:'b2A5',
	8:'b2A6',
	9:'b2A7',
	10:'b2A8',
	11:'b2A9',
	12:'b2B1',
	13:'b2B2',
	14:'b2B3',
	15:'b2B4',
	16:'b2B5',
	17:'b2B6',
	18:'b2B7',
	19:'b2B8',
	20:'b2B9',
	21:'b2C1',
	22:'b1I9',
	23:'b1I5',
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

