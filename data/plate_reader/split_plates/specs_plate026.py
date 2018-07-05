growth_file = "plate026_600.txt"
miller_csv = "plate026_420.txt"
optics_csv = "plate026_550.txt"

raw_growth = "2017.5.5_lacstar_uninducedsamples_plate1_repeat1_OD600.txt"
raw_miller = "2017.5.5_lacstar_uninduced_plate1_repeat1_miller_420.txt"
raw_optics = "2017.5.5_lacstar_uninduced_plate1_repeat1_miller_550.txt"

num_rows = 8
num_cols = 12
date = '17.05.05'

stop_time = 1201
start_time = 0

cultures_dict = {
	0:'RDM',
	1:'b1A5',
	2:'b2D6',
	3:'b1C6',
	4:'b2H2',
	5:'b2E2',
	6:'b2E4',
	7:'b2E3',
	8:'b2F9',
	9:'b2F1',
	10:'b2G2',
	11:'b1D4',
	12:'b2D7',
	13:'b2D2',
	14:'b1C8',
	15:'b2F6',
	16:'b2D9',
	17:'b1B6',
	18:'b1C1',
	19:'b1C7',
	20:'b1B8',
	21:'b1B9',
	22:'b1C9',
	23:'b2C6',
	24:'b1D1'
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

