import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from random import uniform
from datetime import datetime
from scipy import stats
from collections import defaultdict
import itertools
import sys

sns.set(context='poster', style='white', palette='deep', color_codes=True) #, font_scale=1.5)

# t_edges = np.linspace(0, t_max, num=1000)

m_edges = np.logspace(0, 10)

# data = defaultdict(list)
data = []
# for p_r in [0.00005, 0.0005, 0.005, 0.05]:
model = 1
min = 1.8
x_0 = 40
n = 5000
# n = '{:f}'.format(p_r).rstrip('0')

# flavor = 'extant_m{}_{}_{}_{}'.format(model, min, x_0, n)

# flavor = 'min_cope_ext'

# print('Grabbing the {} data.'.format(flavor))

# batch = 0
# # Read in the data
# while True:
# 	try:
# 		curr = pd.read_csv('/Users/Trevor/Desktop/figure data/Tutorial Data/{}_b{}.csv'.format(flavor, batch),
# 			header=None,
# 			names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])
# 		# data[p_r].append(list(map(float, curr['mass'].values.tolist())))
# 		data.append(list(map(float, curr['mass'].values.tolist())))
# 		batch += 1
# 	except FileNotFoundError:
# 		print('Loaded {} files.'.format(batch))
# 		break
# 	except:
# 		print("Unexpected error:", sys.exc_info()[0])
# 		break



# data_b = []
# flavor = 'extant_m{}_{}_{}_{}'.format(model, min, x_0, n)

# print('Grabbing the {} data.'.format(flavor))

# batch = 0
# # Read in the data
# while True:
# 	try:
# 		# curr = pd.read_csv('/Users/Trevor/Desktop/figure data/unconstrainedNicheRuns/{}_b{}.csv'.format(flavor, batch),
# 		curr = pd.read_csv('{}_b{}.csv'.format(flavor, batch),
# 			header=None,
# 			names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])
# 		# data[p_r].append(list(map(float, curr['mass'].values.tolist())))
# 		data_b.append(list(map(float, curr['mass'].values.tolist())))
# 		batch += 1
# 	except FileNotFoundError:
# 		print('Loaded {} files.'.format(batch))
# 		break
# 	except:
# 		print("Unexpected error:", sys.exc_info()[0])
# 		break

# biggest = 0
# for ed in data_b:
# 	for e in ed:
# 		if e > biggest:
# 			biggest = e

# print('Biggest m_min: {}'.format(biggest))

########## For getting aggregate results from MOMsim
data_t = []
data_a = []

model = 1
min = 1.8
x_0 = 40
n = 5000

# flavor = '/Users/Trevor/Desktop/figure data/MOMsim/run 3/extant_m{}_{}_{}_{}'.format(model, min, x_0, n)
flavor = 'extant_m{}_{}_{}_{}no'.format(model, min, x_0, n)

# Read in the data
batch = 0
while True:
	try:
		# 'C:/Users/Trevor/Desktop/Tutorial Data/{}_b{}.csv'
		curr = pd.read_csv('{}_b{}.csv'.format(flavor, batch),
			header=None,
			names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])
		t = curr[curr['niche'] == 0]['mass'] # terrestrial
		a = curr[curr['niche'] == 1]['mass'] # aquatic
		data_t.append(list(map(float, t.values.tolist())))
		data_a.append(list(map(float, a.values.tolist())))
		batch += 1
	except FileNotFoundError:
		print('Loaded {} files.'.format(batch))
		break
	except:
		print("Unexpected error:", sys.exc_info()[0])
		break



def plot_MOM( ):
	fig = sns.plt.figure(figsize=(9, 4.75), dpi=80)
	ax = fig.add_subplot(111)

	mom = pd.read_csv('MOM_data_full.txt', sep=', ', engine='python')

	x_t = mom[mom['land'] == 1]['mass']
	x_a = mom[mom['land'] == 0]['mass']
	x_mom = mom['mass']

	n_mom = len(x_t) + len(x_a)
	m_t_dist = np.histogram(x_t, bins=m_edges, density=True)[0]
	m_a_dist = np.histogram(x_a, bins=m_edges, density=True)[0]
	m_mom_dist = np.histogram(x_mom, bins=m_edges, density=True)[0]

	p_t = ax.scatter(m_edges[:-1], m_t_dist, marker='x', s=42, color='g')
	p_a = ax.scatter(m_edges[:-1], m_a_dist, marker='o', s=42, color='b')
	p_mom = ax.scatter(m_edges[:-1], m_mom_dist, marker='+', s=30, color='r')

	ax.legend([p_t, p_a, p_mom], ['terrestrial', 'aquatic', 'all mammals'])

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	sns.plt.ylim(0.000000000001, sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('frequency (log scale), arb.')
	sns.plt.yticks([])
	sns.plt.tight_layout()
	sns.despine()


def plot_clade_largest( data ):
	big_start = datetime.now()
	print('Starting to plot largest extant masses.')

	largest = []
	for extant_dist in data:
		# print(np.histogram(extant_dist, bins=m_edges)[0])
		largest.append(np.amax(extant_dist))

	fig = sns.plt.figure(figsize=(9, 4.75), dpi=80)
	ax = fig.add_subplot(111)

	ax.scatter(range(len(largest)), largest, s=42, marker='x')

	sns.plt.xlabel('model time')
	sns.plt.xticks([])
	sns.plt.yscale('log')
	sns.plt.ylabel('largest mass seen, g')
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_extant_dist( data, ax, leg='simulation' ):
	big_start = datetime.now()
	print('Starting to plot extant mass distribution.')

	dists = []
	for extant_dist in data:
		# print(np.histogram(extant_dist, bins=m_edges)[0])
		dists.append(np.histogram(extant_dist, bins=m_edges)[0] / len(extant_dist))

	dists = np.array(dists)

	# Take the mean (and quartiles) between the histogram counts
	m_dist = np.mean(dists, axis=0)

	# Grab 95% confidence intervals
	sigma = np.std(dists, axis=0)
	# ci_low = [stats.norm.interval(0.95, loc=m, scale=s)[0] for (m, s) in list(zip(m_dist, sigma))]
	# ci_high = [stats.norm.interval(0.95, loc=m, scale=s)[1] for (m, s) in list(zip(m_dist, sigma))]

	mom = pd.read_csv('MOM_data_full.txt', sep=', ', engine='python')
	x_t = mom[mom['land'] == 1]['mass']
	x_a = mom[mom['land'] == 0]['mass']
	# x_mom = mom['mass']
	n_mom = len(x_t) + len(x_a)
	m_t_dist = np.histogram(x_t, bins=m_edges)[0] / n_mom
	m_a_dist = np.histogram(x_a, bins=m_edges)[0] / n_mom
	# m_mom_dist = np.histogram(x_mom, bins=m_edges)[0]
	p_t = ax.scatter(m_edges[:-1], m_t_dist, marker='x', s=30, color='g')
	p_a = ax.scatter(m_edges[:-1], m_a_dist, marker='o', s=30, color='b')
	# p_mom = ax.scatter(m_edges[:-1], m_mom_dist, marker='x', s=30, color='r')

	# p_m = ax.scatter(m_edges[:-1], m_dist, marker='D', s=42, color='k')

	ax.legend([p_t, p_a], ['terrestrial', 'aquatic', leg])

	# ax.plot(m_edges[1:-1], ci_low[1:], linewidth=0.75, ls='--', color='k')
	# ax.plot(m_edges[1:-1], ci_high[1:], linewidth=0.75, ls='--', color='k')
	

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	# (1.0 / 5000.0)
	# 0.000000000001
	# 0.9
	sns.plt.ylim((1.0 / 5000.0), sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('frequency (log scale), arb.')
	sns.plt.yticks([])
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_extant_dist_density( data, ax, leg='simulation' ):
	big_start = datetime.now()
	print('Starting to plot extant mass distribution.')

	dists = []
	for extant_dist in data:
		# print(np.histogram(extant_dist, bins=m_edges)[0])
		dists.append(np.histogram(extant_dist, bins=m_edges, density=True)[0])

	dists = np.array(dists)

	# Take the mean (and quartiles) between the histogram counts
	m_dist = np.mean(dists, axis=0)

	# Grab 95% confidence intervals
	sigma = np.std(dists, axis=0)
	ci_low = [stats.norm.interval(0.95, loc=m, scale=s)[0] for (m, s) in list(zip(m_dist, sigma))]
	ci_high = [stats.norm.interval(0.95, loc=m, scale=s)[1] for (m, s) in list(zip(m_dist, sigma))]

	mom = pd.read_csv('MOM_data_full.txt', sep=', ', engine='python')
	x_t = mom[mom['land'] == 1]['mass']
	# x_a = mom[mom['land'] == 0]['mass']
	# x_mom = mom['mass']
	# n_mom = len(x_t) + len(x_a)
	m_t_dist = np.histogram(x_t, bins=m_edges, density=True)[0]
	# m_a_dist = np.histogram(x_a, bins=m_edges, density=True)[0]
	# m_mom_dist = np.histogram(x_mom, bins=m_edges)[0]
	p_t = ax.scatter(m_edges[:-1], m_t_dist, marker='x', s=30, color='g')
	# p_a = ax.scatter(m_edges[:-1], m_a_dist, marker='o', s=30, color='b')
	# p_mom = ax.scatter(m_edges[:-1], m_mom_dist, marker='x', s=30, color='r')

	p_m = ax.scatter(m_edges[:-1], m_dist, marker='D', s=42, color='k')

	ax.legend([p_t, p_m], ['MOM terrestrial', leg])

	ax.plot(m_edges[1:-1], ci_low[1:], linewidth=0.75, ls='--', color='k')
	ax.plot(m_edges[1:-1], ci_high[1:], linewidth=0.75, ls='--', color='k')
	

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	# (1.0 / 5000.0)
	# 0.000000000001
	# 0.9
	sns.plt.ylim(0.000000000001, sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('density (log scale), arb.')
	sns.plt.yticks([])
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_extant_dist_MOMsim( data_t, data_a, ax, leg='simulation' ):
	big_start = datetime.now()
	print('Starting to plot extant mass distribution.')

	dists_t = []
	dists_a = []
	for (ed_t, ed_a) in list(zip(data_t, data_a)):
		# print(np.histogram(ed_t, bins=m_edges)[0])
		# print(np.histogram(ed_a, bins=m_edges)[0])
		n_tot = len(ed_a) + len(ed_t)
		# print(n_tot)
		dists_t.append(np.histogram(ed_t, bins=m_edges)[0] / n_tot)
		dists_a.append(np.histogram(ed_a, bins=m_edges)[0] / n_tot)
	
	dists_t = np.array(dists_t)
	dists_a = np.array(dists_a)

	# Take the mean (and quartiles) between the histogram counts
	m_dist_t = np.mean(dists_t, axis=0)
	m_dist_a = np.mean(dists_a, axis=0)

	# Grab 95% confidence intervals
	sigma_t = np.std(dists_t, axis=0)
	ci_low_t = [stats.norm.interval(0.95, loc=m, scale=s)[0] for (m, s) in list(zip(m_dist_t, sigma_t))]
	ci_high_t = [stats.norm.interval(0.95, loc=m, scale=s)[1] for (m, s) in list(zip(m_dist_t, sigma_t))]

	sigma_a = np.std(dists_a, axis=0)
	ci_low_a = [stats.norm.interval(0.95, loc=m, scale=s)[0] for (m, s) in list(zip(m_dist_a, sigma_a))]
	ci_high_a = [stats.norm.interval(0.95, loc=m, scale=s)[1] for (m, s) in list(zip(m_dist_a, sigma_a))]

	mom = pd.read_csv('MOM_data_full.txt', sep=', ', engine='python')
	x_t = mom[mom['land'] == 1]['mass']
	x_a = mom[mom['land'] == 0]['mass']
	# x_mom = mom['mass']
	n_mom = len(x_t) + len(x_a)
	m_t_dist = np.histogram(x_t, bins=m_edges)[0] / n_mom
	m_a_dist = np.histogram(x_a, bins=m_edges)[0] / n_mom
	# m_mom_dist = np.histogram(x_mom, bins=m_edges)[0]
	p_t = ax.scatter(m_edges[:-1], m_t_dist, marker='x', s=30, color='g')
	p_a = ax.scatter(m_edges[:-1], m_a_dist, marker='o', s=30, color='b')
	# p_mom = ax.scatter(m_edges[:-1], m_mom_dist, marker='x', s=30, color='r')

	p_m_t = ax.scatter(m_edges[:-1], m_dist_t, marker='D', s=42, color='k')
	p_m_a = ax.scatter(m_edges[:-1], m_dist_a, marker='D', s=42, color='k')

	ax.legend([p_t, p_a, p_m_t, p_m_a], ['MOM terrestrial', 'MOM aquatic', 'sim terrestrial', 'sim aquatic'])

	ax.plot(m_edges[1:-1], ci_low_t[1:], linewidth=0.75, ls='--', color='0.75')
	ax.plot(m_edges[1:-1], ci_high_t[1:], linewidth=0.75, ls='--', color='0.75')

	ax.plot(m_edges[1:-1], ci_low_a[1:], linewidth=0.75, ls='--', color='0.75')
	ax.plot(m_edges[1:-1], ci_high_a[1:], linewidth=0.75, ls='--', color='0.75')
	

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	# (1.0 / 5000.0)
	# 0.000000000001
	# 0.9
	sns.plt.ylim((1.0 / 5000.0), sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('density (log scale), arb.')
	sns.plt.yticks([])
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_extant_niche_dists( data ):
	big_start = datetime.now()
	print('Starting to plot extant mass distribution.')

	dists = []
	for row in list(zip(data['mass'], data['niche'])):
		d = pd.DataFrame(data=np.asarray(row).T.tolist(), columns=['mass', 'niche'])
		niches = d['niche'].unique()
		dists.append([np.histogram(d[d['niche'] == niche], bins=m_edges, density=True)[0] for niche in niches])

	dists = np.array(dists)

	# Take the mean (and quartiles) between the histogram counts
	m_dist = np.mean(dists, axis=0)
	print(m_dist)

	# Grab 95% confidence intervals
	sigma = np.std(dists, axis=0)
	ci_low = [stats.norm.interval(0.95, loc=m, scale=s)[0] for (m, s) in list(zip(m_dist, sigma))]
	ci_high = [stats.norm.interval(0.95, loc=m, scale=s)[1] for (m, s) in list(zip(m_dist, sigma))]


	fig = sns.plt.figure(figsize=(9, 4.75), dpi=80)
	ax = fig.add_subplot(111)

	mom = pd.read_csv('MOM_data_full.txt', sep=', ', engine='python')
	# x_t = mom[mom['land'] == 1]['mass']
	x_a = mom[mom['land'] == 0]['mass']
	# x_mom = mom['mass']
	# n_mom = len(x_t) + len(x_a)
	# m_t_dist = np.histogram(x_t, bins=m_edges, density=True)[0]
	m_a_dist = np.histogram(x_a, bins=m_edges, density=True)[0]
	# m_mom_dist = np.histogram(x_mom, bins=m_edges)[0]
	# p_t = ax.scatter(m_edges[:-1], m_t_dist, marker='x', s=30, color='g')
	p_a = ax.scatter(m_edges[:-1], m_a_dist, marker='o', s=30, color='b')
	# p_mom = ax.scatter(m_edges[:-1], m_mom_dist, marker='x', s=30, color='r')

	p_m = ax.scatter(m_edges[:-1], m_dist, marker='D', s=42, color='k')

	# ax.legend([p_t, p_m], ['MOM terrestrial', flavor])

	ax.plot(m_edges[1:-1], ci_low[1:], linewidth=0.75, ls='--', color='k')
	ax.plot(m_edges[1:-1], ci_high[1:], linewidth=0.75, ls='--', color='k')
	

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	# (1.0 / 5000.0)
	# 0.000000000001
	# 0.9
	sns.plt.ylim(0.000000000001, sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('density (log scale), arb.')
	sns.plt.yticks([])
	sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_lambda_dist( ):
	x = pd.read_csv('lambda_dist.csv', header=None, names=['x'])

	(dist, edges) = np.histogram(x, bins=42)

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)

	ax.scatter(edges[:-1], dist)
	sns.plt.xlabel('lambda draw')
	sns.plt.ylabel('frequency')
	# sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.despine()


def plot_group_lifetimes( data ):
	big_start = datetime.now()
	print('Starting to plot lifetimes of {} clades.'.format(len(data)))

	mins = []
	lifetimes = []
	first = 0
	for m in data:
		x_mins = m['m_min'].unique()
		for x_min in x_mins[first:]:
			first = 1
			clade = m[m['m_min'] == x_min]

			# first birth
			fb = np.amin(clade['birth'].values.tolist())
			ld = np.amax(clade['death'].values.tolist())
			mins.append(x_min)
			lifetimes.append(ld - fb)

	g = sns.JointGrid(mins, lifetimes, size=8, dropna=False, xlim=(1, 10**7), ylim=(1000, 10**5))
	g = g.plot_joint(plt.scatter, color='m', edgecolor='white')
	g.ax_marg_x.hist(mins, bins=m_edges, color='b', alpha=0.6)
	g.ax_marg_y.hist(lifetimes, orientation='horizontal', bins=m_edges, color='r', alpha=0.6)
	sns.plt.xscale('log')
	sns.plt.yscale('log')
	# ax.scatter(range(len(lifetimes)), lifetimes, marker='x', s=42)
	sns.plt.xlabel('group m_min')
	sns.plt.ylabel('group lifetime, model steps')
	sns.plt.tight_layout()
	sns.despine()


	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_extant_scatter( data, ax, leg='simulation' ):
	big_start = datetime.now()
	print('Starting to plot extant mass distribution.')

	dists = []
	for extant_dist in data:
		# print(np.histogram(extant_dist, bins=m_edges)[0])
		dists.append(np.histogram(extant_dist, bins=m_edges)[0] / len(extant_dist))

	dists = np.array(dists)

	# Take the mean (and quartiles) between the histogram counts
	m_dist = np.mean(dists, axis=0)

	# Grab 95% confidence intervals
	sigma = np.std(dists, axis=0)
	ci_low = [stats.norm.interval(0.95, loc=m, scale=s)[0] for (m, s) in list(zip(m_dist, sigma))]
	ci_high = [stats.norm.interval(0.95, loc=m, scale=s)[1] for (m, s) in list(zip(m_dist, sigma))]

	mom = pd.read_csv('MOM_data_full.txt', sep=', ', engine='python')
	x_t = mom[mom['land'] == 1]['mass']
	x_a = mom[mom['land'] == 0]['mass']
	# x_mom = mom['mass']
	n_mom = len(x_t) + len(x_a)
	m_t_dist = np.histogram(x_t, bins=m_edges)[0] / n_mom
	# m_a_dist = np.histogram(x_a, bins=m_edges, density=True)[0]
	# m_mom_dist = np.histogram(x_mom, bins=m_edges)[0]
	p_t = ax.scatter(m_edges[:-1], m_t_dist, marker='x', s=30, color='g')
	# p_a = ax.scatter(m_edges[:-1], m_a_dist, marker='o', s=30, color='b')
	# p_mom = ax.scatter(m_edges[:-1], m_mom_dist, marker='x', s=30, color='r')

	p_m = ax.scatter(m_edges[:-1], m_dist, marker='D', s=42, color='k')

	ax.legend([p_t, p_m], ['MOM terrestrial', leg])

	ax.plot(m_edges[1:-1], ci_low[1:], linewidth=0.75, ls='--', color='k')
	ax.plot(m_edges[1:-1], ci_high[1:], linewidth=0.75, ls='--', color='k')
	

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	# (1.0 / 5000.0)
	# 0.000000000001
	# 0.9
	sns.plt.ylim((1.0 / 5000.0), sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('frequency (log scale), arb.')
	sns.plt.yticks([])
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


# plot_group_lifetimes(data[0.0005])

fig = sns.plt.figure(figsize=(9, 4.75), dpi=80)
ax = fig.add_subplot(111)
# plot_extant_dist_MOMsim(data_t, data_a, ax)
plot_extant_dist(data, ax, leg=flavor)

# ax = fig.add_subplot(212)
# plot_extant_dist_density(data_b, ax, leg='expanding nichespace')

# cnt = 1
# for p_r in [0.00005, 0.0005, 0.005, 0.05]:
# 	ax = fig.add_subplot(220+cnt)
# 	pdata = data[p_r]
# 	plot_extant_dist_density(pdata, ax, leg='p_r = {:f}'.format(p_r).rstrip('0'))
# 	cnt += 1




# plot_MOM()
# plot_schematic()
# plot_extant_dist(data)
# plot_extant_niche_dists(data)
# plot_clade_largest(data)
# plot_lambda_dist()
# plot_extant_dist_w_MOM()

sns.plt.show()
# sns.plt.savefig("../../Thesis/figs/tutorial{}.pdf".format(flavor))