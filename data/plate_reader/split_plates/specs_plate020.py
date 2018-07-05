growth_file = "plate020_600.txt"
miller_csv = "plate020_420.txt"
optics_csv = "plate020_550.txt"

raw_growth = "2017.4.14_lacstar_c61&71r18L10lib_uninduced_plate1_repeat2_OD600.txt"
raw_miller = "2017.4.14_lacstar_c61&71r18L10lib_uninduced_plate1_repeat2_miller_420.txt"
raw_optics = "2017.4.14_lacstar_c61&71r18L10lib_uninduced_plate1_repeat2_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.04.14'

stop_time = 1201
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b2C6',
	3:'b2D2',
	4:'b2D6',
	5:'b2D7',
	6:'b2D9',
	7:'b2E2',
	8:'b2E3',
	9:'b2E4',
	10:'b2F1',
	11:'b2H1',
	12:'b2H2'
	}
	
num_strains = max(cultures_dict.keys())
	
culture_locations_string = '''
1	0	5	0	9	0	1	0	5	0	9	0
1	3	5	7	9	11	1	3	5	7	9	11
1	3	5	7	9	11	1	3	5	7	9	11
0	3	0	7	0	11	0	3	0	7	0	11
2	0	6	0	10	0	2	0	6	0	10	0
2	4	6	8	10	12	2	4	6	8	10	12
2	4	6	8	10	12	2	4	6	8	10	12
0	4	0	8	0	12	0	4	0	8	0	12
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

