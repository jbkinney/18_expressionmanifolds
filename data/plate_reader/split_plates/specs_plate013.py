growth_file = "plate013_600.txt"
miller_csv = "plate013_420.txt"
optics_csv = "plate013_550.txt"

raw_growth = "2017.4.6_lacstar_c61&71r18L10lib_plate3_repeat2_OD600.txt"
raw_miller = "2017.4.6_lacstar_c61&71r18L10lib_plate3_repeat2_miller_420.txt"
raw_optics = "2017.4.6_lacstar_c61&71r18L10lib_plate3_repeat2_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.04.06'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	2:'b2F1',
	3:'b2F2',
	4:'b2F3',
	5:'b2F4',
	6:'b2F5',
	7:'b2F6',
	8:'b2F7',
	9:'b2F8',
	10:'b2F9',
	11:'b2G1',
	12:'b2G2',
	1:'b1A5'
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
0	0	0	0	0	0	250	250	250	250	250	250
0	0	0	0	0	0	250	250	250	250	250	250
0	0	0	0	0	0	250	250	250	250	250	250
0	0	0	0	0	0	250	250	250	250	250	250
0	0	0	0	0	0	250	250	250	250	250	250
0	0	0	0	0	0	250	250	250	250	250	250
0	0	0	0	0	0	250	250	250	250	250	250
0	0	0	0	0	0	250	250	250	250	250	250
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

