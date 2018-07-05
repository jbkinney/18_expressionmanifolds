plate_name = 'plate037u'

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2017.5.25_lacstar_c71r18l10lib_misc_plate1_uninduced_315_OD600.txt"
raw_miller = "2017.5.25_lacstar_c71r18l10lib_misc_plate1_uninduced_miller_420.txt"
raw_optics = "2017.5.25_lacstar_c71r18l10lib_misc_plate1_uninduced_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.05.25'

stop_time = 600
start_time = 15

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b2H1',
	3:'b2H2',
	4:'b2H3',
	5:'b2H4',
	6:'b2H5',
	7:'b2H6',
	8:'b2H7',
	9:'b2H8',
	10:'b1B2',
	11:'b1E8',
	12:'b1E6',
	13:'b1A3',
	14:'b1E3',
	15:'b1E5',
	16:'b1A6',
	17:'b1A9',
	18:'b3C7',
	19:'b3C9',
	20:'b1A7',
	21:'b1E7',
	22:'b1A1',
	23:'b1I1',
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

