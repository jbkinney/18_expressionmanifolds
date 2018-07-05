plate_name = 'plate064i'

induced = True

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2017_10_27_barcodelibtest_2_plate2_induced_310_OD600.txt"
raw_miller = "2017_10_27_barcodelibtest_2_plate2_induced_miller_1h_420.txt"
raw_optics = "2017_10_27_barcodelibtest_2_plate2_induced_miller_1h_550.txt"

num_rows = 8
num_cols = 12
date = '17.10.27'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b8E7',
	3:'b8F2',
	4:'b8F4',
	5:'b8F6',
	6:'b8F8',
	7:'b8E6',
	8:'b8F1',
	9:'b8F3',
	10:'b8F5',
	11:'b8F7',
	12:'b8E4',
	13:'b8E5',
	14:'b8E3',
	15:'b8E1',
	16:'b8D8',
	17:'b8D6',
	18:'b8D4',
	19:'b8E2',
	20:'b8D9',
	21:'b8D7',
	22:'b8D5',
	23:'b8D3',
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

if induced: camp_concentrations_string = camp_concentrations_string_i
else: camp_concentrations_string = camp_concentrations_string_u

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

