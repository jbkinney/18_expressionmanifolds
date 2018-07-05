growth_file = "plate011_600.txt"
miller_csv = "plate011_420.txt"
optics_csv = "plate011_550.txt"

raw_growth = "2017.4.4_lacstar_c61r18L10lib_plate1_repeat2_OD600.txt"
raw_miller = "2017.4.4_lacstar_c61r18L10lib_plate1_repeat2_miller_420.txt"
raw_optics = "2017.4.4_lacstar_c61r18L10lib_plate1_repeat2_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.04.04'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	2:'b2C6',
	3:'b2C7',
	4:'b2C8',
	5:'b2C9',
	6:'b2D1',
	7:'b2D2',
	8:'b2D3',
	9:'b2D4',
	10:'b2D5',
	11:'b2D6',
	12:'b2D7',
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

