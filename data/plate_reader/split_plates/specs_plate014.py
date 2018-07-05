growth_file = "plate014_600.txt"
miller_csv = "plate014_420.txt"
optics_csv = "plate014_550.txt"

raw_growth = "2017.4.6_lacstar_c71r18L10lib_plate4_repeat2_OD600.txt"
raw_miller = "2017.4.6_lacstar_c71r18L10lib_plate4_repeat2_miller_420.txt"
raw_optics = "2017.4.6_lacstar_c71r18L10lib_plate4_repeat2_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.04.06'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	2:'b2G3',
	3:'b2G4',
	4:'b2G5',
	5:'b2G6',
	6:'b2G7',
	7:'b2G8',
	8:'b2G9',
	9:'b2H1',
	10:'b2H2',
	11:'b2H3',
	12:'b2H4',
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

