plate_name = 'plate076i'

induced = True

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2017_11_16_lacstar_c63c64libs_plate2_induced_3h_OD600_flipped"
raw_miller = "2017_11_16_lacstar_c63c64libs_plate2_induced_miller_1h_420.txt"
raw_optics = "2017_11_16_lacstar_c63c64libs_plate2_induced_miller_1h_550.txt"

num_rows = 8
num_cols = 12
date = '17.11.16'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b9D1',
	3:'b9D2',
	4:'b9D3',
	5:'b9D4',
	6:'b9D5',
	7:'b9D7',
	8:'b9D8',
	9:'b9D9',
	10:'b9E1',
	11:'b9E2',
	12:'b9E3',
	13:'b9E4',
	14:'b9E5',
	15:'b9E6',
	16:'b9E7',
	17:'b9E8',
	18:'b9F1',
	19:'b9F2',
	20:'b9F3',
	21:'b9F4',
	22:'b9F5',
	23:'b9F6',
	24:'b1F1'
	}
	
num_strains = max(cultures_dict.keys())
	
culture_locations_string = '''
24	0	20	0	16	0	12	0	8	0	4	0
24	22	20	18	16	14	12	10	8	6	4	2
24	22	20	18	16	14	12	10	8	6	4	2
0	22	0	18	0	14	0	10	0	6	0	2
23	0	19	0	15	0	11	0	7	0	3	0
23	21	19	17	15	13	11	9	7	5	3	1
23	21	19	17	15	13	11	9	7	5	3	1
0	21	0	17	0	13	0	9	0	5	0	1
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

