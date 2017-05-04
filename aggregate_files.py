import os

model = 1
min = 1.8
x_0 = 40
n = 5000

flavor = '/Users/Trevor/Desktop/figure data/unconstrainedNicheRuns/extant_m{}_{}_{}_{}'.format(model, min, x_0, n)

# Move the data
batch = 0
while True:
	try:
		# 'C:/Users/Trevor/Desktop/Tutorial Data/{}_b{}.csv'
		curr = pd.read_csv('{}_b{}.csv'.format(flavor, batch),
			header=None,
			names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])
		t = curr[curr['niche'] == 0] # terrestrial
		a = curr[curr['niche'] == 1] # aquatic
		data_t.append(list(map(float, t['mass'].values.tolist())))
		data_a.append(list(map(float, a['mass'].values.tolist())))
	except :
		print("Unexpected error:", sys.exc_info()[0])
		break