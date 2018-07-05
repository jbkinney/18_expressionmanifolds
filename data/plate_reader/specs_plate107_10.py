plate_name = 'plate107_10'

induced = False

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2018_03_18_lacstar_c60c82_10uM_induced_plate_OD600_3h.txt"
raw_miller = "2018_03_18_lacstar_c60c82_10uM_induced_plate_miller_10h_420.txt"
raw_optics = "2018_03_18_lacstar_c60c82_10uM_induced_plate_miller_10h_550.txt"

num_rows = 8
num_cols = 12
date = '18.03.19'

stop_time = 600
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b7D8',
	3:'b7D9',
	4:'b7E1',
	5:'b7E2',
	6:'b7E3',
	7:'b7E4',
	8:'b7E5',
	9:'b7E6',
	10:'b7E7',
	11:'b7E8',
	12:'b7E9',
	13:'b13H4',
	14:'b13H5',
	15:'b13H6',
	16:'b13H7',
	17:'b13H8',
	18:'b13H9',
	19:'b13I1',
	20:'b13I2',
	21:'b13I3',
	22:'b13I4',
	23:'b13I5',
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

camp_concentrations_string_i = '''
25	25	25	25	25	25	25	25	25	25	25	25
25	25	25	25	25	25	25	25	25	25	25	25
25	25	25	25	25	25	25	25	25	25	25	25
25	25	25	25	25	25	25	25	25	25	25	25
25	25	25	25	25	25	25	25	25	25	25	25
25	25	25	25	25	25	25	25	25	25	25	25
25	25	25	25	25	25	25	25	25	25	25	25
25	25	25	25	25	25	25	25	25	25	25	25
'''

camp_concentrations_string_u = '''
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
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

