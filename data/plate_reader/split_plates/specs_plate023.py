growth_file = "plate023_600.txt"
miller_csv = "plate023_420.txt"
optics_csv = "plate023_550.txt"

raw_growth = "2017.4.21_lacstar_c61&71r18L10&c61r18L35lib_uninduced_plate2_repeat2_OD600.txt"
raw_miller = "2017.4.21_lacstar_c61&71r18L10&c61r18L35lib_uninduced_plate2_repeat2_420.txt"
raw_optics = "2017.4.21_lacstar_c61&71r18L10&c61r18L35lib_uninduced_plate2_repeat2_550.txt"

num_rows = 8
num_cols = 12
date = '17.04.21'

stop_time = 60
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b2D4',
	3:'b2D8',
	4:'b2E6',
	5:'b2E7',
	6:'b2F6',
	7:'b2G2',
	8:'b2G5',
	9:'b2G6',
	10:'b2G9',
	11:'b2H1',
	12:'b2H5',
	13:'b1B6',
	14:'b1B8',
	15:'b1B9',
	16:'b1C1',
	17:'b1C2',
	18:'b1C4',
	19:'b1C6',
	20:'b1C7',
	21:'b1C8',
	22:'b1C9',
	23:'b1D1',
	24:'b1D4'
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

