plate_name = 'plate094_2'

induced = False

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2018_02_20_lacstar_c61c71selectlib_plate_2uM_induced_3h_OD600.txt"
raw_miller = "2018_02_20_lacstar_c61c71selectlib_plate_2uM_induced_miller_10h_420.txt"
raw_optics = "2018_02_20_lacstar_c61c71selectlib_plate_2uM_induced_miller_10h_550.txt"

num_rows = 8
num_cols = 12
date = '18.02.20'

stop_time = 600
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b1C3',
	3:'b1E2',
	4:'b2C7',
	5:'b2C8',
	6:'b2C9',
	7:'b2D1',
	8:'b2D3',
	9:'b2D5',
	10:'b2E1',
	11:'b2E8',
	12:'b2F2',
	13:'b1I8',
	14:'b2A5',
	15:'b2B1',
	16:'b2F3',
	17:'b2F5',
	18:'b2G3',
	19:'b2G4',
	20:'b2G7',
	21:'b2G8',
	22:'b2H3',
	23:'b2H8',
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

camp_concentrations_string_u = '''
2	2	2	2	2	2	2	2	2	2	2	2
2	2	2	2	2	2	2	2	2	2	2	2
2	2	2	2	2	2	2	2	2	2	2	2
2	2	2	2	2	2	2	2	2	2	2	2
2	2	2	2	2	2	2	2	2	2	2	2
2	2	2	2	2	2	2	2	2	2	2	2
2	2	2	2	2	2	2	2	2	2	2	2
2	2	2	2	2	2	2	2	2	2	2	2
'''

camp_concentrations_string_i = '''
109	109	109 109	109	109	109 109	109	109	109 109
109	109	109 109	109	109	109 109	109	109	109 109
109	109	109 109	109	109	109 109	109	109	109 109
109	109	109 109	109	109	109 109	109	109	109 109
109	109	109 109	109	109	109 109	109	109	109 109
109	109	109 109	109	109	109 109	109	109	109 109
109	109	109 109	109	109	109 109	109	109	109 109
109	109	109 109	109	109	109 109	109	109	109 109
'''

if induced: camp_concentrations_string = camp_concentrations_string_i
else: camp_concentrations_string = camp_concentrations_string_u

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

