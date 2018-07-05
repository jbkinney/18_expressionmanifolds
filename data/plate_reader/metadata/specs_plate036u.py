plate_name = 'plate036u'

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2017.5.24_lacstar_c61&71r18l10lib_plate1_uninduced_315_OD600.txt"
raw_miller = "2017.5.24_lacstar_c61&71r18l10lib_plate1_uninduced_miller_420.txt"
raw_optics = "2017.5.24_lacstar_c61&71r18l10lib_plate1_uninduced_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.05.24'

stop_time = 600
start_time = 5

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b2E6',
	3:'b2E7',
	4:'b2E8',
	5:'b2E9',
	6:'b2F1',
	7:'b2F2',
	8:'b2F3',
	9:'b2F4',
	10:'b2F5',
	11:'b2F6',
	12:'b2F7',
	13:'b2F8',
	14:'b2F9',
	15:'b2G1',
	16:'b2G2',
	17:'b2G3',
	18:'b2G4',
	19:'b2G5',
	20:'b2G6',
	21:'b2G7',
	22:'b2G8',
	23:'b2G9',
	24:'b1F1'
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

