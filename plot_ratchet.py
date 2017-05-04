import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from random import uniform
from datetime import datetime
import itertools
import sys
import csv

sns.set(context='poster', style='white', palette='muted', color_codes=True) #, font_scale=1.5)

model = 1
min = 1.8
x_0 = 40
n = 5000

all_species = pd.read_csv('simMOMFig_all.csv'.format(model, min, x_0, n),
# all_species = pd.read_csv('simMeteorBiasedFig_all.csv',
	header=None,
	names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])

# all_species_b = pd.read_csv('all_m{}_{}_{}_{}.csv'.format(model, min, x_0, n),
# # all_species_b = pd.read_csv('simRadiationFig_p0.5_all.csv',
# 	header=None,
# 	names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])

n_max = len(all_species.index)
t_max = np.amax(all_species['birth'])

m_edges = np.logspace(0, 10)
t_edges = np.linspace(0, t_max, num=1000)

def plot_random_walks( walks, steps, ax ):
	print('Plotting {} random walks of {} steps each'.format(walks, steps))

	def step_gen( tot ):
		x = 0
		for _ in range(tot):
			yield x
			x += uniform(-1, 1)

	ax.plot([0, steps], [0, 0], linewidth=0.5, color="k", linestyle="--")

	for _ in range(walks):
		ax.plot(range(steps), list(step_gen(steps)), linewidth=1.0)

	sns.plt.xlabel('iteration')
	sns.plt.xticks([])
	sns.plt.ylabel('position, arb.')
	sns.plt.yticks([])
	sns.plt.tight_layout()
	sns.despine()


def plot_random_walks_from_data( m, walks ):
	print('Plotting {} random walks from data'.format(walks))

	# Grab extant species for choosing which one to trace lineage of
	extant_species = pd.read_csv('simMeteorBiasedFig_extant.csv'.format(model, min, x_0, n),
		header=None,
		names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])

	fig = sns.plt.figure(figsize=(9, (4.75 * 0.4 * walks)), dpi=80)
	pal = itertools.cycle(sns.color_palette())

	selection = extant_species.sample(n=walks)
	print('{}'.format(selection))

	for ii in range(walks):
		ax = fig.add_subplot(walks, 1, (ii + 1))
		species = selection.iloc[ii]
		x = []

		while species['ancestor'] > 0:
			# Build the list backwards, since we're starting at the end
			x.insert(0, species['mass'])
			species = m[m['id'] == species['ancestor']].iloc[0]

		x.insert(0, x_0)
		ax.plot(range(len(x)), x, linewidth=1.0, color=next(pal))
		sns.plt.ylabel('mass, g')
		sns.plt.yscale('log')
		sns.plt.xticks([])


	sns.plt.xlabel('model time')
	sns.plt.tight_layout()
	sns.despine()


def plot_extant_dist( ax ):
	big_start = datetime.now()
	print('Starting to plot extant mass distribution.')

	extant_species = pd.read_csv('extant_m{}_{}_{}_{}.csv'.format(model, min, x_0, n),
		header=None,
		names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])

	x = list(map(float, extant_species['mass'].values.tolist()))

	m_dist = np.histogram(x, bins=m_edges, density=True)[0]
	ax.scatter(m_edges[:-1], m_dist, marker='D', s=42)
	
	sns.plt.xscale('log')
	sns.plt.yscale('log')
	sns.plt.ylim(0.000000000001, sns.plt.ylim()[1])
	sns.plt.yticks([])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('density (log scale)')
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_niche_extant_dists( ax ):
	big_start = datetime.now()
	print('Starting to plot niche-wise extant mass distributions.')

	extant_species = pd.read_csv('extant_m{}_{}_{}_{}.csv'.format(model, min, x_0, n),
		header=None,
		names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])

	niches = extant_species['niche'].unique()

	pal = itertools.cycle(sns.color_palette())

	lw = 2.0
	for niche in niches:
		clade = extant_species[extant_species['niche'] == niche]
		
		m_dist  = np.histogram(clade['mass'], bins=m_edges)[0] / len(extant_species.index)
		ax.plot(m_edges[:-1], m_dist, marker='D', ms=4, lw=lw, color=next(pal))
		lw = 1.0

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	sns.plt.yticks([])
	sns.plt.ylim((1 / len(extant_species.index)), sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('frequency (log scale), arb.')
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_extant_dist_w_MOM( ax ):
	big_start = datetime.now()
	print('Starting to plot terrestrial and aquatic MOM mass distribution.')

	extant_species = pd.read_csv('simMOMFig_extant.csv'.format(model, min, x_0, n),
		header=None,
		names=['id', 'birth', 'mass', 'm_min', 'death', 'ancestor', 'niche'])

	niches = extant_species['niche'].unique()
	pal = itertools.cycle(sns.color_palette())

	mom = pd.read_csv('MOM_data_full.txt', sep=', ', engine='python')
	x_t = mom[mom['land'] == 1]['mass']
	x_a = mom[mom['land'] == 0]['mass']
	n_mom = len(mom)

	m_t_dist = np.histogram(x_t, bins=m_edges)[0] / n_mom
	m_a_dist = np.histogram(x_a, bins=m_edges)[0] / n_mom
	p_t = ax.scatter(m_edges[:-1], m_t_dist, marker='x', s=42, color='g')
	p_a = ax.scatter(m_edges[:-1], m_a_dist, marker='o', s=42, color='b')

	for niche in niches:
		clade = extant_species[extant_species['niche'] == niche]
		m_dist = np.histogram(clade['mass'], bins=m_edges)[0] / len(extant_species)
		p_x = ax.scatter(m_edges[:-1], m_dist, marker='D', s=36, color='k')

	ax.legend([p_t, p_a, p_x], ['terrestrial', 'aquatic', 'simulation'])

	x_mins = extant_species['m_min'].unique()
	# sns.plt.text(0.1, 0.1, '{} groups represented in simulation results'.format(len(x_mins)))
	print('{} groups represented in simulation results'.format(len(x_mins)))

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	sns.plt.ylim((1.0 / 5000.0), sns.plt.ylim()[1])
	sns.plt.yticks([])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('frequency (log scale), arb.')
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_dist( x ):
	big_start = datetime.now()
	print('Starting to plot all_species mass distribution.')

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)

	m_dist = np.histogram(x, bins=m_edges)[0]
	ax.scatter(m_edges[:-1], m_dist, marker='D', s=42)

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	sns.plt.ylim(0.9, sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('frequency')
	sns.plt.title('all species that ever lived')
	sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_dists( m ):
	big_start = datetime.now()
	print('Starting to plot clade-wise mass distribution.')

	x_mins = m['m_min'].unique()
	cutoff = 500

	print('There are {} clades total.'.format(len(x_mins)))

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)
	pal = itertools.cycle(sns.color_palette())

	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]
		if len(clade.index) > cutoff:
			m_dist  = np.histogram(clade['mass'], bins=m_edges)[0]
			ax.scatter(m_edges[:-1], m_dist, marker='D', color=next(pal))

	sns.plt.xscale('log')
	sns.plt.yscale('log')
	sns.plt.ylim(0.9, sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('proportion')
	sns.plt.title('clades with over {} species'.format(cutoff))
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_largest( m, ax ):
	big_start = datetime.now()
	print('Starting to plot largest mass ever seen in clade over time.')

	x_mins = m['m_min'].unique()

	lw = 2.0
	c = 'g'
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]
		masses = clade['mass'].tolist()
		for ii in range(len(masses) - 1):
			if masses[ii + 1] < masses[ii]:
				masses[ii + 1] = masses[ii]

		ax.plot(clade['birth'], masses, linewidth=lw)
		lw = 1.0
		c = 'b'

	sns.plt.xlabel('model time')
	sns.plt.xticks([])
	sns.plt.yscale('log')
	sns.plt.ylabel('largest mass seen, g')
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


def plot_clade_largest_extant_downsample( m ):
	big_start = datetime.now()
	print('Starting to plot largest extant over time, through downsampling.')

	x_mins = m['m_min'].unique()

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)
	
	lw = 2.0
	for x_min in x_mins[0:1]:
		start = datetime.now()

		clade = m[m['m_min'] == x_min]
		n_clade = len(clade.index)
		print('Working on x_min = {}; Clade has {} species.'.format(x_min, n_clade))

		largest = np.zeros((len(t_edges), 1))
		last_largest_extant_fallback = 0
		seen = False
		ii = 0
		for t in t_edges:
			sys.stdout.write('\r')
			percent = ii / len(t_edges)
			sys.stdout.write('[{:50}] {:05.2f}%'.format('='*int(percent*50), 100*percent))
			sys.stdout.flush()

			first_catch = True
			largest_extant = 0
			for jj in range(last_largest_extant_fallback, n_clade - 1):
				if clade['birth'].iloc[jj] <= t and clade['death'].iloc[jj] >= t:
					if first_catch:
						last_largest_extant_fallback = jj
						first_catch = False

					if clade['mass'].iloc[jj] > largest_extant:
						largest_extant = clade['mass'].iloc[jj]

				# If this species' birthday is past time t, then all following b-days will be too
				elif clade['birth'].iloc[jj] > t:
					break

			largest[ii] = largest_extant

			# If we've seen a member of the clade, its reign has begun
			if largest_extant > 0:
				seen = True

			# If we have seen the clade and now it had died off, move on to next clade
			if seen and largest_extant == 0:
				break

			ii += 1

		print('Plotting x_min = {}.  {} seconds elapsed for this clade'.format(x_min, (datetime.now() - start).total_seconds()))
		ax.plot(t_edges, largest, linewidth=lw)
		lw = 1.0

	sns.plt.xlabel('model time')
	sns.plt.ylabel('largest of extant species')
	sns.plt.title('largest of extant species over time, by clade')
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_largest_wdist( m ):
	big_start = datetime.now()
	print('Starting to plot largest mass in clade over time, with x_min distribution.')

	x_mins = m['m_min'].unique()

	g = sns.JointGrid(x='id', y='m_min', data=x_mins)

	lw = 2.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]
		masses = clade['mass'].tolist()
		for ii in range(len(masses) - 1):
			if masses[ii + 1] < masses[ii]:
				masses[ii + 1] = masses[ii]

		g.ax_joint.plot(clade['birth'], masses, linewidth=lw)
		lw = 1.0

	g.ax_marg_x.set_axis_off()
	g.ax_marg_y.hist(x_mins, orientation='horizontal', bins=m_edges, color='b', alpha=.6)
	sns.plt.xlabel('x_mins')
	sns.plt.yscale('log')
	g.set_axis_labels('model time', 'largest mass in clade, g')
	sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_n_extant( m, ax ):
	big_start = datetime.now()
	print('Starting to plot n extant over time.')

	x_mins = m['m_min'].unique()
	
	lw = 2.0
	for x_min in x_mins:
		start = datetime.now()

		n_extant = np.zeros((len(t_edges), 1))
		clade = m[m['m_min'] == x_min]
		n_clade = len(clade.index)
		print('Working on x_min = {}; Clade has {} species.'.format(x_min, n_clade))

		bds = list(zip(clade['birth'].tolist(), clade['death'].tolist()))

		seen = False
		for (count, t) in zip(n_extant, t_edges):
			if n_clade > 5000:
				if t == t_edges[250]:
					print('Running time {}, 25% done.  {} seconds elapsed for this clade'.format(t, (datetime.now() - start).total_seconds()))
				if t == t_edges[500]:
					print('Running time {}, 50% done.  {} seconds elapsed for this clade'.format(t, (datetime.now() - start).total_seconds()))
				if t == t_edges[750]:
					print('Running time {}, 75% done.  {} seconds elapsed for this clade'.format(t, (datetime.now() - start).total_seconds()))
			
			for (birth, death) in bds:
				if birth <= t and death >= t:
					count += 1
				# If this species' birthday is past time t, then all following b-days will be too
				elif birth > t:
					break

			# If we've seen a member of the clade, its reign has begun
			if count > 0:
				seen = True
			# If we have seen the clade and now it had died off, move on to next clade
			if seen and count == 0:
				break

			# With heuristics: 3981.3 sec
			# Without heuristics, with detailed status: 4108.3
			# Without heuristics: 4102.9

		print('Plotting x_min = {}.  {} seconds elapsed for this clade'.format(x_min, (datetime.now() - start).total_seconds()))
		ax.plot(t_edges, n_extant, linewidth=lw)
		lw = 1.0

	sns.plt.xlabel('model time')
	sns.plt.xticks([])
	sns.plt.yscale('log')
	sns.plt.ylabel('largest mass seen, g')
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_sizes( m ):
	big_start = datetime.now()
	print('Starting to plot sizes of clades.')

	x_mins = m['m_min'].unique()

	clade_sizes = [len(m[m['m_min'] == x_min].index) for x_min in x_mins]

	g = sns.JointGrid(x_mins, clade_sizes, dropna=False, xlim=(1, 10**8), ylim=(0.8, 10**7))
	g = g.plot_joint(plt.scatter, color='m', edgecolor='white')
	g.ax_marg_x.hist(x_mins, bins=m_edges, color='b', alpha=0.6)
	g.ax_marg_y.hist(clade_sizes, orientation='horizontal', bins=m_edges, color='r', alpha=0.6)
	sns.plt.xscale('log')
	sns.plt.yscale('log')
	#sns.plt.subplots_adjust(top=0.9)
	g.fig.suptitle('sizes of {} clades'.format(len(x_mins)))
	g.set_axis_labels('minimum size for clade', 'total species in clade, over all time')

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_spawnrate( m, ax ):
	big_start = datetime.now()
	print('Starting to plot clade spawn rates.')

	x_mins = m['m_min'].unique()

	lw = 2.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]

		births = list(map(float, clade['birth'].values.tolist()))
		bd_dist  = np.histogram(births, bins=t_edges)[0]

		ax.plot(t_edges[:-1], bd_dist, linewidth=lw)
		lw = 1.0

	sns.plt.xlabel('model time')
	sns.plt.xticks([])
	sns.plt.ylabel('speciation rate (new spec per {:.0f} steps)'.format(t_edges[1]))
	# sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.plt.tight_layout()
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_extrate( m, ax ):
	big_start = datetime.now()
	print('Starting to plot clade extinction rates.')

	x_mins = m['m_min'].unique()

	lw = 2.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]

		deaths = list(map(float, clade['death'].values.tolist()))
		bd_dist  = np.histogram(deaths, bins=t_edges)[0]

		ax.plot(t_edges[:-1], bd_dist, linewidth=lw)
		lw = 1.0

	sns.plt.xlabel('model time')
	sns.plt.xticks([])
	sns.plt.ylabel('extinction rate (new spec per {:.0f} steps)'.format(t_edges[1]))
	# sns.plt.gcf().subplots_adjust(bottom=0.15)
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


def plot_group_lifetimes( m, ax ):
	big_start = datetime.now()
	print('Starting to plot sizes of clades.')

	x_mins = m['m_min'].unique()

	lifetimes = []
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]

		# first birth
		fb = np.amin(clade['birth'].values.tolist())
		ld = mp.amax(clade['death'].values.tolist())
		lifetimes.append(ld - fb)

	g = sns.JointGrid(x_mins, lifetimes, dropna=False, xlim=(1, 10**10), ylim=(0.8, 10**7))
	g = g.plot_joint(plt.scatter, color='m', edgecolor='white')
	g.ax_marg_x.hist(x_mins, bins=m_edges, color='b', alpha=0.6)
	g.ax_marg_y.hist(clade_sizes, orientation='horizontal', bins=m_edges, color='r', alpha=0.6)
	sns.plt.xscale('log')
	sns.plt.yscale('log')
	# ax.scatter(range(len(lifetimes)), lifetimes, marker='x', s=42)
	sns.plt.xlabel('group m_min')
	sns.plt.ylabel('group lifetime, model steps')
	sns.plt.tight_layout()
	sns.despine()


	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))



def get_axis_limits(ax, ylog=False):
	scale = 0.9

	x = ax.get_xlim()[0] + ((1.0 - scale) * ax.get_xlim()[1] * 0.25)

	if ylog:
		y = ax.get_ylim()[1] ** scale
	else:
		y = ax.get_ylim()[1] * scale

	return x, y


fig = sns.plt.figure(figsize=(9, 14.25), dpi=80)
# ax = fig.add_subplot(111)
# plot_niche_extant_dists(ax)


# plot_random_walks_from_data(all_species, 3)

ax = fig.add_subplot(311)
plot_extant_dist_w_MOM(ax)
ax.annotate('a)', xy=get_axis_limits(ax, ylog=True))

ax = fig.add_subplot(312)
plot_clade_largest_extant(all_species, ax)
ax.annotate('b)', xy=get_axis_limits(ax, ylog=True))

ax = fig.add_subplot(313)
plot_clade_n_from_rate(all_species, ax)
ax.annotate('c)', xy=get_axis_limits(ax))

# ax = fig.add_subplot(211)
# plot_clade_largest(all_species, ax)
# ax.annotate('b)', xy=get_axis_limits(ax, ylog=True))

# ax = fig.add_subplot(212)
# plot_clade_n_from_rate(all_species, ax)
# ax.annotate('c)', xy=get_axis_limits(ax))

### Figure choices ###
# plot_random_walks(3, 3000)
# plot_random_walks_from_data(all_species, 3)
# plot_dist(all_species['mass'])
# plot_extant_dist()
# plot_niche_extant_dists()
# plot_extant_dist_w_MOM()
# plot_clade_dists(all_species)
# plot_clade_largest(all_species)
# plot_clade_largest_extant(all_species)
# plot_clade_largest_extant_downsample(all_species)
# plot_clade_largest_wdist(all_species)
# plot_clade_n_extant(all_species)
# plot_clade_sizes(all_species)
# plot_clade_spawnrate(all_species)
# plot_clade_extrate(all_species)
# plot_clade_largest(all_species, ax)
# plot_clade_largest_extant(all_species)
# plot_clade_largest_extant_downsample(all_species)
# plot_clade_largest_wdist(all_species)
# plot_clade_n_extant(all_species)
# plot_clade_sizes(all_species)
# plot_clade_spawnrate(all_species)
# plot_clade_extrate(all_species)
# plot_clade_n_from_rate(all_species, ax)

sns.plt.show()
# sns.plt.savefig("something.pdf")