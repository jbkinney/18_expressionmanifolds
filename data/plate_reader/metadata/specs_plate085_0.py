plate_name = 'plate085_0'

induced = 0

growth_file = plate_name + "_600.txt"
miller_csv = plate_name + "_420.txt"
optics_csv = plate_name + "_550.txt"

raw_growth = "2018_01_25_lacstar_crp_titration_plate1_uninduced_3h_OD600.txt"
raw_miller = "2018_01_25_lacstar_crp_titration_plate1_uninduced_miller_230h_420.txt"
raw_optics = "2018_01_25_lacstar_crp_titration_plate1_uninduced_miller_230h_550.txt"

num_rows = 8
num_cols = 12
date = '18.01.25'

stop_time = 150
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b5H6',
	2:'b5H7',
	3:'b5G9',
	4:'b5H3',
	5:'b5H2',
	6:'b5H1',
	7:'b5G4',
	8:'b5F2',
	9:'b1D9',
	10:'b1C2',
	11:'b1C4',
	12:'b1D1',
	13:'b1D3',
	14:'b1B9',
	15:'b1D5',
	16:'b1D6',
	17:'b2A6',
	18:'b2B7',
	19:'b2B4',
	20:'b2B5',
	21:'b2B6',
	22:'b2B8',
	23:'b2B2',
	24:'b2A2'
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

camp_concentrations_string_0 = '''
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0
'''

camp_concentrations_string_10 = '''
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
10	10	10	10	10	10	10	10	10	10	10	10
'''

camp_concentrations_string_5 = '''
5	5	5	5	5	5	5	5	5	5	5	5
5	5	5	5	5	5	5	5	5	5	5	5
5	5	5	5	5	5	5	5	5	5	5	5
5	5	5	5	5	5	5	5	5	5	5	5
5	5	5	5	5	5	5	5	5	5	5	5
5	5	5	5	5	5	5	5	5	5	5	5
5	5	5	5	5	5	5	5	5	5	5	5
5	5	5	5	5	5	5	5	5	5	5	5
'''

if induced==0: camp_concentrations_string = camp_concentrations_string_0
elif induced==5: camp_concentrations_string = camp_concentrations_string_5
elif induced==10: camp_concentrations_string = camp_concentrations_string_10


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

