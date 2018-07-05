plate_name = 'plate034u'

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2017.5.19_lacstar_c71r18l35lib_plate1_uninduced_305_OD600.txt"
raw_miller = "2017.5.19_lacstar_c71r18l35lib_plate1_uninduced_miller_420.txt"
raw_optics = "2017.5.19_lacstar_c71r18l35lib_plate1_uninduced_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.05.19'

stop_time = 600
start_time = 30

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b1I5',
	3:'b1I6',
	4:'b1I7',
	5:'b1I8',
	6:'b1I9',
	7:'b2A1',
	8:'b2A2',
	9:'b2A3',
	10:'b2A4',
	11:'b2A5',
	12:'b2A6',
	13:'b2A7',
	14:'b2A8',
	15:'b2A9',
	16:'b2B1',
	17:'b2B2',
	18:'b2B3',
	19:'b2B4',
	20:'b2B5',
	21:'b2B6',
	22:'b2B7',
	23:'b2B8',
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

