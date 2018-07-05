growth_file = "plate009_600.txt"
miller_csv = "plate009_420.txt"
optics_csv = "plate009_550.txt"

raw_growth = "2017.3.20_lacstar_c71r18L35lib_plate3_OD600.txt"
raw_miller = "2017.3.21_lacstar_c71r18L35lib_plate3_miller_420.txt"
raw_optics = "2017.3.21_lacstar_c71r18L35lib_plate3_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.03.21'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	2:'b2B3',
	3:'b2B4',
	4:'b2B5',
	5:'b2B6',
	6:'b2B7',
	7:'b2B8',
	8:'b2B9',
	9:'b2C1',
	1:'b1A5'
	}
	
num_strains = max(cultures_dict.keys())
	
culture_locations_string = '''
1	0	4	0	7	0	1	0	4	0	7	0
1	0	4	0	7	0	1	0	4	0	7	0
1	3	4	6	7	9	1	3	4	6	7	9
0	3	0	6	0	9	0	3	0	6	0	9
2	3	5	6	8	9	2	3	5	6	8	9
2	0	5	0	8	0	2	0	5	0	8	0
2	0	5	0	8	0	2	0	5	0	8	0
0	0	0	0	0	0	0	0	0	0	0	0
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

