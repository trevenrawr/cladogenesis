import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import itertools
sns.set(context='poster', style='white', palette='muted', color_codes=True) #, font_scale=1.5)

model = 3
min = 1.8
x_0 = 40
n = 5000

all_species = pd.read_csv('all_ecm{}_{}_{}_{}.csv'.format(model, min, x_0, n),
	header=None,
	names=['birth', 'mass', 'm_min', 'death', 'ancestor'])

n_max = len(all_species.index)
t_max = (n_max - 1) / 2

m_edges = np.logspace(0, 10)
t_edges = np.linspace(0, t_max, num=1000)

to_time = lambda x: (x - 1) / 2


def plot_extant_dist( ):
	big_start = datetime.now()
	print('Starting to plot extant mass distribution.')

	extant_species = pd.read_csv('extant_ecm{}_{}_{}_{}.csv'.format(model, min, x_0, n),
		header=None,
		names=['birth', 'mass', 'm_min', 'death', 'ancestor'])

	x = extant_species['mass']

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)

	m_dist = np.histogram(x, bins=m_edges)[0]
	ax.scatter(m_edges[:-1], m_dist, marker='D', s=42)
	
	sns.plt.xscale('log')
	sns.plt.yscale('log')
	sns.plt.ylim(0.9, sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('frequency')
	sns.plt.title('extant species at simulation termination')
	sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))



def plot_extant_dist_w_MOM( ):
	big_start = datetime.now()
	print('Starting to plot terrestrial and aquatic MOM mass distribution.')


	mom = pd.read_csv('MOM_data_full.txt', sep=None)

	extant_species = pd.read_csv('extant_ecm{}_{}_{}_{}.csv'.format(model, min, x_0, n),
		header=None,
		names=['birth', 'mass', 'm_min', 'death', 'ancestor'])

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


def plot_clade_largest( m ):
	big_start = datetime.now()
	print('Starting to plot largest mass in clade over time.')

	x_mins = m['m_min'].unique()

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)

	lw = 2.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]
		masses = clade['mass'].tolist()
		for ii in range(len(masses) - 1):
			if masses[ii + 1] < masses[ii]:
				masses[ii + 1] = masses[ii]

		ax.plot(to_time(clade['birth']), masses, linewidth=lw)
		lw = 1.0

	sns.plt.yscale('log')
	sns.plt.xlabel('model time')
	sns.plt.ylabel('largest mass in clade')
	sns.plt.title('largest species seen to date, by clade')
	sns.plt.gcf().subplots_adjust(bottom=0.15)
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

		g.ax_joint.plot(to_time(clade['birth']), masses, linewidth=lw)
		lw = 1.0

	g.ax_marg_x.set_axis_off()
	g.ax_marg_y.hist(x_mins, orientation='horizontal', bins=m_edges, color='b', alpha=.6)
	sns.plt.xlabel('x_mins')
	sns.plt.yscale('log')
	g.set_axis_labels('model time', 'largest mass in clade, g')
	sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_n_extant( m ):
	big_start = datetime.now()
	print('Starting to plot n extant over time.')

	x_mins = m['m_min'].unique()

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)
	
	lw = 2.0
	for x_min in x_mins:
		start = datetime.now()

		n_extant = np.zeros((len(t_edges), 1))
		clade = m[m['m_min'] == x_min]
		n_clade = len(clade.index)
		print('Working on x_min = {}; Clade has {} species.'.format(x_min, n_clade))

		bds = list(zip(list(map(to_time, clade['birth'].tolist())), clade['death'].tolist()))

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
	sns.plt.ylabel('number of extant species')
	sns.plt.title('extant species over time, by clade')
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


def plot_clade_spawnrate( m ):
	big_start = datetime.now()
	print('Starting to plot clade spawn rates.')

	x_mins = m['m_min'].unique()

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)

	lw = 2.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]

		births = list(map(to_time, clade['birth'].tolist()))
		bd_dist  = np.histogram(births, bins=t_edges)[0]

		ax.plot(t_edges[:-1], bd_dist, linewidth=lw)
		lw = 1.0

	sns.plt.xlabel('model time')
	sns.plt.ylabel('speciation rate (new spec per {:.0f} steps)'.format(t_edges[1]))
	sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_extrate( m ):
	big_start = datetime.now()

	print('Starting to plot clade extinction rates.')

	x_mins = m['m_min'].unique()

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)

	lw = 2.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]

		deaths = list(clade['death'].tolist())
		bd_dist  = np.histogram(deaths, bins=t_edges)[0]


		ax.plot(t_edges[:-1], bd_dist, linewidth=lw)
		lw = 1.0

	sns.plt.xlabel('model time')
	sns.plt.ylabel('extinction rate (new spec per {:.0f} steps)'.format(t_edges[1]))
	sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))


def plot_clade_n_from_rate( m ):
	big_start = datetime.now()
	print('Starting to plot clade extant counts from spawn and death rates.')

	x_mins = m['m_min'].unique()

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)

	lw = 2.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]

		deaths = list(clade['death'].tolist())
		d_dist  = np.histogram(deaths, bins=t_edges)[0]

		births = list(map(to_time, clade['birth'].tolist()))
		b_dist  = np.histogram(births, bins=t_edges)[0]

		bd_dist = b_dist - d_dist

		last_n = 0
		n_extant = []
		for bd in bd_dist:
			n_extant.append(last_n + bd / 2)
			last_n = last_n + bd / 2

		ax.plot(t_edges[:-1], n_extant, linewidth=lw)
		lw = 1.0

	sns.plt.xlabel('model time')
	sns.plt.ylabel('number of extant species, by group'.format(t_edges[1]))
	sns.plt.gcf().subplots_adjust(bottom=0.15)
	sns.despine()

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))



# plot_dist(all_species['mass'])
# plot_extant_dist()
# plot_extant_dist_w_MOM()
# plot_clade_dists(all_species)
# plot_clade_largest(all_species)
# plot_clade_largest_wdist(all_species)
# plot_clade_n_extant(all_species)
# plot_clade_sizes(all_species)
plot_clade_spawnrate(all_species)
plot_clade_extrate(all_species)
plot_clade_n_from_rate(all_species)
sns.plt.show()