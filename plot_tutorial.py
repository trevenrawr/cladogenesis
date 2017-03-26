import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from random import uniform
from datetime import datetime
from scipy import stats
import itertools
import sys

sns.set(context='poster', style='white', palette='deep', color_codes=True) #, font_scale=1.5)

# n_max = len(all_species.index)
# t_max = (n_max - 1) / 2

# t_edges = np.linspace(0, t_max, num=1000)

m_edges = np.logspace(0, 10)

flavor = 'basic'
data = [];

# Read in the data
batch = 0
while True:
	try:
		curr = pd.read_csv('/Users/Trevor/Desktop/Tutorial Data/{}_b{}.csv'.format(flavor, batch),
			header=None,
			names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])
		data.append(list(map(float, curr['mass'].values.tolist())))
		batch += 1
	except:
		# if batch == 0:
			# raise ValueError('Did not load any files.')
			# print('Did not load any files.')
		break

print('Finished loading data!  Loaded {} files.'.format(batch))


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


def plot_schematic( ):
	fig = sns.plt.figure(figsize=(9, 4.75), dpi=80)
	ax = fig.add_subplot(111)

	x = np.array(range(1000))

	ax.plot(x, 100 - np.exp(-x), lw=2.0, color='c')
	# ax.plot(x, y, lw=2.0, color='m')
	# ax.plot(x, y, lw=2.0, color='y')

	sns.plt.yscale('log')
	sns.plt.ylim((1.0 / 5000.0), sns.plt.ylim()[1])
	sns.plt.xlabel('time')
	sns.plt.ylabel('mass (log scale)')
	sns.plt.yticks([])
	sns.plt.xticks([])
	sns.plt.tight_layout()
	sns.despine()


def plot_extant_dist( data ):
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


	fig = sns.plt.figure(figsize=(9, 4.75), dpi=80)
	ax = fig.add_subplot(111)

	mom = pd.read_csv('MOM_data_full.txt', sep=', ', engine='python')
	x_t = mom[mom['land'] == 1]['mass']
	x_a = mom[mom['land'] == 0]['mass']
	# x_mom = mom['mass']
	n_mom = len(x_t) + len(x_a)
	m_t_dist = np.histogram(x_t, bins=m_edges, density=True)[0]
	# m_a_dist = np.histogram(x_a, bins=m_edges)[0] / n_mom
	# m_mom_dist = np.histogram(x_mom, bins=m_edges)[0] / n_mom
	p_t = ax.scatter(m_edges[:-1], m_t_dist, marker='x', s=30, color='g')
	# p_a = ax.scatter(m_edges[:-1], m_a_dist, marker='o', s=30, color='b')
	# p_mom = ax.scatter(m_edges[:-1], m_mom_dist, marker='x', s=30, color='r')

	p_m = ax.scatter(m_edges[:-1], m_dist, marker='D', s=42, color='k')

	ax.legend([p_t, p_m], ['MOM terrestrial', flavor])

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


def plot_extant_dist_w_MOM( ):
	big_start = datetime.now()
	print('Starting to plot terrestrial and aquatic MOM mass distribution.')


	mom = pd.read_csv('MOM_data_full.txt', sep=None)

	extant_species = pd.read_csv('extant_m{}_{}_{}_{}.csv'.format(model, min, x_0, n),
		header=None,
		names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])

	x = extant_species['mass']

	x_t = mom[mom['land'] == 1]['mass']
	x_a = mom[mom['land'] == 0]['mass']

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)
	pal = sns.color_palette()

	m_t_dist = np.histogram(x_t, bins=m_edges)[0]
	m_a_dist = np.histogram(x_a, bins=m_edges)[0]
	m_dist = np.histogram(x, bins=m_edges)[0]
	p_t = ax.scatter(m_edges[:-1], m_t_dist, marker='D', s=42, color=pal[1])
	p_a = ax.scatter(m_edges[:-1], m_a_dist, marker='o', s=42, color=pal[0])
	p_x = ax.scatter(m_edges[:-1], m_dist, marker='s', s=36, color=pal[2])

	ax.legend([p_t, p_a, p_x], ['terrestrial', 'aquatic', 'simulation'])

	x_mins = extant_species['m_min'].unique()
	sns.plt.text(0.1, 0.1, '{} groups represented in simulation results'.format(len(x_mins)))
	print('{} groups represented in simulation results'.format(len(x_mins)))

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	sns.plt.ylim(0.9, sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('frequency')
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


# plot_MOM()
# plot_schematic()
plot_extant_dist(data)
# plot_lambda_dist()
# plot_extant_dist_w_MOM()

sns.plt.show()
# sns.plt.savefig("../../Thesis/figs/tutorial{}.pdf".format(flavor))