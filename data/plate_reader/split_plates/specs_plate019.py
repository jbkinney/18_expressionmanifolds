growth_file = "plate019_600.txt"
miller_csv = "plate019_420.txt"
optics_csv = "plate019_550.txt"

raw_growth = "2017.3.23_lacstar_uninduced_lib_samples_plate2_OD600.txt"
raw_miller = "2017.3.24_lacstar_uninduced_lib_samples_plate2_miller_420.txt"
raw_optics = "2017.3.24_lacstar_uninduced_lib_samples_plate2_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.03.24'

stop_time = 1201
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b2D4',
	3:'b2D8',
	4:'b2E6',
	5:'b2E7',
	6:'b2G5',
	7:'b2G6',
	8:'b2G9',
	9:'b2H1',
	10:'b2H5',
	11:'b1C2',
	12:'b1C4',
	13:'b1C5',
	14:'b1D2',
	15:'b1D5'
	}
	
num_strains = max(cultures_dict.keys())
	
culture_locations_string = '''
1	0	5	0	9	0	13	0	0	0	0	0
1	3	5	7	9	11	13	15	0	0	0	0
1	3	5	7	9	11	13	15	0	0	0	0
0	3	0	7	0	11	0	15	0	0	0	0
2	0	6	0	10	0	14	0	0	0	0	0
2	4	6	8	10	12	14	0	0	0	0	0
2	4	6	8	10	12	14	0	0	0	0	0
0	4	0	8	0	12	0	0	0	0	0	0
'''

camp_concentrations_string = '''
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
'''

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

