growth_file = "plate024_600.txt"
miller_csv = "plate024_420.txt"
optics_csv = "plate024_550.txt"

raw_growth = "2017.5.4_lacstar_minus10nulls_plate1_repeat1_OD600.txt"
raw_miller = "2017.5.4_lacstar_minus10nulls_plate1_repeat1_miller_420.txt"
raw_optics = "2017.5.4_lacstar_minus10nulls_plate1_repeat1_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.05.04'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b2I4',
	3:'b2I8',
	4:'b3A1',
	5:'b3A3',
	6:'b3A5',
	7:'b3A7',
	8:'b3A9',
	9:'b3B2',
	10:'b3B4',
	11:'b3B6',
	12:'b3B8'
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

