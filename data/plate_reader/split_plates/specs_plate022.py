growth_file = "plate022_600.txt"
miller_csv = "plate022_420.txt"
optics_csv = "plate022_550.txt"

raw_growth = "2017.4.20_lacstar_c-r18L35lib_plate2_repeat1_OD600.txt"
raw_miller = "2017.4.20_lacstar_c-r18L35lib_plate2_repeat1_miller_420.txt"
raw_optics = "2017.4.20_lacstar_c-r18L35lib_plate2_repeat1_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.04.20'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	2:'b1F5',
	3:'b1F6',
	4:'b1F7',
	5:'b1F8',
	6:'b1F9',
	7:'b1G1',
	8:'b1G2',
	9:'b1G3',
	10:'b1G4',
	11:'b1G5',
	12:'b1G6',
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

