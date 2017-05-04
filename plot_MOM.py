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




########## For getting aggregate results from MOMsim

model = 1
min = 1.8
x_0 = 40
n = 5000

all_species = pd.read_csv('simMOMFig_all.csv'.format(model, min, x_0, n),
	header=None,
	names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])

t_max = np.amax(all_species['birth'])
m_edges = np.logspace(0, 10)
t_edges = np.linspace(0, t_max, num=1000)

# flavor = '/Users/Trevor/Desktop/figure data/MOMsim/run 3/extant_m{}_{}_{}_{}'.format(model, min, x_0, n)
flavor = 'extant_m{}_{}_{}_{}'.format(model, min, x_0, n)

# Read in the batched MOM data
data_t = []
data_a = []

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

	p_m_t = ax.scatter(m_edges[:-1], m_dist_t, marker='D', s=42, color='#06470c')
	p_m_a = ax.scatter(m_edges[:-1], m_dist_a, marker='D', s=42, color='#0504aa')

	ax.legend([p_t, p_a, p_m_t, p_m_a], ['MOM terrestrial', 'MOM aquatic', 'sim terrestrial', 'sim aquatic'])

	ax.plot(m_edges[1:-1], ci_low_t[1:], linewidth=0.75, ls='--', color='0.5')
	ax.plot(m_edges[1:-1], ci_high_t[1:], linewidth=0.75, ls='--', color='0.5')

	ax.plot(m_edges[1:-1], ci_low_a[1:], linewidth=0.75, ls='--', color='0.5')
	ax.plot(m_edges[1:-1], ci_high_a[1:], linewidth=0.75, ls='--', color='0.5')
	

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


def plot_clade_largest_extant( m, ax ):
	big_start = datetime.now()
	print('Starting to plot largest extant mass in clade over time.')

	x_mins = m['m_min'].unique()

	c = 'g'
	counter = 0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]
		clade_indexes = list(clade.index.values)

		print('Starting clade {} of {}.  Contains {} species.'.format(counter, len(x_mins), len(clade)))

		file_found = False
		try:
			time_saver = pd.read_csv('largest_extant_m{}_{}_{}_{}_niche{}.csv'.format(model, min, x_0, n, counter),
				header=None,
				names=['events', 'largest'])
		except FileNotFoundError:
			print("Did not load any files.")
		except:
			print("Unexpected error:", sys.exc_info()[0])
		else:
			print("Found a time-saving file!  Loaded for niche {}.".format(counter))
			events = time_saver['events'].values.tolist()
			largest = time_saver['largest'].values.tolist()
			file_found = True

		counter += 1

		# Plot the largest seen to date version
		largest_masses = clade['mass'].tolist()
		for ii in range(len(largest_masses) - 1):
			if largest_masses[ii + 1] < largest_masses[ii]:
				largest_masses[ii + 1] = largest_masses[ii]
		ax.plot(clade['birth'], largest_masses, linewidth=2.0, c=c)

		# Plot all species lines
		# for _, species in clade.iterrows():
		# 	ax.plot([species['birth'], species['death']], [species['mass'], species['mass']], linewidth=1.0)
		
		if not file_found:
			masses = clade['mass'].tolist()
			births = clade['birth'].tolist()
			deaths = clade['death'].tolist()
			events = [births[0]]
			largest = [masses[0]]
			last_largest_death = deaths[0]
			last_largest_extant_fallback = 0

			# TODO: Add in reading from file to save calculation time.

			# Iterate chronologically through the clade
			for ii in range(len(masses) - 1):
				sys.stdout.write('\r')
				percent = ii / len(masses)
				sys.stdout.write('[{:50}] {:06.2f}%'.format('='*int(percent*50), 100*percent))
				sys.stdout.flush()

				largest_extant = 0
				
				# If the next species evolves after the currently largest one has died,
				#  then we need to find out what to do in the intervening time,
				#  that is to say, find the largest extant species at the time of death
				if births[ii + 1] >= last_largest_death:
					# print('Currently {}; Largest species died at {}'.format(births[ii + 1], last_largest_death))
					# Drop in a point to signify the end of the last largest's life
					events.append(last_largest_death - 1)
					largest.append(largest[-1])

					first_catch = True
					largest_extant_index = 0
					last_index = 0
					for jj in range(last_largest_extant_fallback, len(masses) - 1):
						if clade['birth'].iloc[jj] > last_largest_death:
							# We've passed all the species that could satisfy the requirements
							break

						# print('{} <? {} <? {}'.format(clade['birth'].iloc[jj], last_largest_death, clade['death'].iloc[jj]))
						if clade['birth'].iloc[jj] <= last_largest_death and clade['death'].iloc[jj] > last_largest_death:
							if clade['mass'].iloc[jj] > largest_extant:
								# print('YES!')
								# print('Species {} is larger than {}g, at a mass of {}.'.format(index, largest_extant, species['mass']))
								largest_extant = clade['mass'].iloc[jj]
								largest_extant_index = jj
								if first_catch:
									# This was the first index that was alive the last time we needed to search for the
									#  largest extant species after one dies off.  We can start our search here instead of
									#  at the beginning next time and save LOTS of time.  (Hopefully!)
									last_largest_extant_fallback = jj
									# f = open('fallbacks.txt', 'a')
									# f.write(str(clade['id'].iloc[jj]) + '\n')
									# f.close()
									# print('First hit: {}'.format(clade['id'].iloc[jj]))
									first_catch = False

					if largest_extant_index > 0:
						events.append(last_largest_death)
						largest.append(largest_extant)
						last_largest_death = clade['death'].iloc[largest_extant_index]
						# print('\nLast largest died at {}.  Next largest extant index: {}, dies at {}'.format(events[-1], largest_extant_index, last_largest_death))
					else:
						# events.append(last_largest_death)
						# largest.append(0)
						print('\nClade has died off???  Currently {}; last largest died: {}'.format(births[ii + 1], last_largest_death))


				# If the next species is larger than the last lasrgest, insert a point
				#  at its birthdate and its mass
				if masses[ii + 1] >= largest[-1]:
					last_largest_death = deaths[ii + 1]
					events.append(births[ii + 1] - 1)
					largest.append(largest[-1])
					events.append(births[ii + 1])
					largest.append(masses[ii + 1])
				elif largest_extant == 0:
					events.append(births[ii + 1])
					largest.append(largest[len(largest) - 1])

			with open('largest_extant_m{}_{}_{}_{}_niche{}.csv'.format(model, min, x_0, n, ii), 'w') as f:
				writer = csv.writer(f, delimiter=',', lineterminator='\n')
				for row in list(zip(events, largest)):
					writer.writerow(row)
			print("Saved calculated largest_extant to csv.")


		ax.plot(events, largest, linewidth=1.0, c=c)
		c = 'b'

	sns.plt.xlabel('model time')
	sns.plt.xticks([])
	sns.plt.yscale('log')
	sns.plt.ylabel('largest mass seen, g')
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_n_from_rate( m, ax ):
	big_start = datetime.now()
	print('Starting to plot clade extant counts from spawn and death rates.')

	x_mins = m['m_min'].unique()

	lw = 2.0
	c = 'g'
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]

		deaths = list(map(float, clade['death'].values.tolist()))
		d_dist  = np.histogram(deaths, bins=t_edges)[0]

		births = list(map(float, clade['birth'].values.tolist()))
		b_dist  = np.histogram(births, bins=t_edges)[0]

		bd_dist = b_dist - d_dist

		last_n = 0
		n_extant = []
		for bd in bd_dist:
			n_extant.append(last_n + bd)
			last_n = last_n + bd

		ax.plot(t_edges[:-1], n_extant, linewidth=lw, c=c)
		lw = 1.0
		c = 'b'

		# print('Clade with m_min = {} had at most {} species alive at once.'.format(x_min, max(n_extant)))

	sns.plt.xlabel('model time')
	sns.plt.xticks([])
	sns.plt.ylabel('number of extant species'.format(t_edges[1]))
	# sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def get_axis_limits(ax, xlog=False, ylog=False):
	scale = 0.9

	if xlog:
		# x = ax.get_xlim()[0] + ((1.0 - scale) ** ax.get_xlim()[1] * 0.25)
		x = 1
	else:
		x = ax.get_xlim()[0] + ((1.0 - scale) * ax.get_xlim()[1] * 0.25)

	if ylog:
		y = ax.get_ylim()[1] ** scale
	else:
		y = ax.get_ylim()[1] * scale

	return x, y


fig = sns.plt.figure(figsize=(9, 14.25), dpi=80)

ax = fig.add_subplot(311)
plot_extant_dist_MOMsim(data_t, data_a, ax)
ax.annotate('a)', xy=(0.5, 0.1))

ax = fig.add_subplot(312)
plot_clade_largest_extant(all_species, ax)
ax.annotate('b)', xy=get_axis_limits(ax, ylog=True))

ax = fig.add_subplot(313)
plot_clade_n_from_rate(all_species, ax)
ax.annotate('c)', xy=get_axis_limits(ax))


sns.plt.show()