import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import itertools
sns.set(style='white', palette='muted', color_codes=True)

model = 3
min = 1.8
x_0 = 40
n = 5000

all_species = pd.read_csv('all_m{}_{}_{}_{}.csv'.format(model, min, x_0, n),
	header=None,
	names=['birth', 'mass', 'm_min', 'death', 'ancestor'])

n_max = len(all_species.index)
t_max = (n_max - 1) / 2

m_edges = np.logspace(0, 10)
t_edges = np.linspace(0, t_max, num=1000)

to_time = lambda x: (x - 1) / 2


def plot_dist( x ):
	sns.plt.figure()
	m_dist = np.histogram(x, bins=m_edges, density=True)[0]
	ax = sns.distplot(x, norm_hist=1, bins=m_edges, kde=False)
	ax.scatter(m_edges[:-1], m_dist, marker='D')
	sns.plt.xscale('log')
	sns.plt.yscale('log')
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('proportion')
	sns.plt.title('all species that ever lived')


def plot_clade_dists( m ):
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
	sns.plt.ylim(1, sns.plt.ylim()[1])
	sns.plt.xlabel('species mass, g')
	sns.plt.ylabel('proportion')
	sns.plt.title('clades with over {} species'.format(cutoff))


def plot_clade_largest( m ):
	x_mins = m['m_min'].unique()

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)

	lw = 1.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]
		masses = clade['mass'].tolist()
		for ii in range(len(masses) - 1):
			if masses[ii + 1] < masses[ii]:
				masses[ii + 1] = masses[ii]

		ax.plot(clade['birth'], masses, linewidth=lw)
		lw = 0.5

	sns.plt.yscale('log')
	sns.plt.xlabel('model time')
	sns.plt.ylabel('largest mass in clade')
	sns.plt.title('largest species seen to date, by clade')


def plot_clade_largest_wdist( m ):
	x_mins = m['m_min'].unique()

	g = sns.JointGrid(x='id', y='m_min', data=x_mins)

	lw = 1.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]
		masses = clade['mass'].tolist()
		for ii in range(len(masses) - 1):
			if masses[ii + 1] < masses[ii]:
				masses[ii + 1] = masses[ii]

		g.ax_joint.plot(clade['birth'], masses, linewidth=lw)
		lw = 0.5

	g.ax_marg_x.set_axis_off()
	g.ax_marg_y.hist(x_mins, orientation='horizontal', bins=m_edges, color='b', alpha=.6)
	sns.plt.xlabel('x_mins')
	sns.plt.yscale('log')
	g.set_axis_labels('model time', 'largest mass in clade, g')


def plot_clade_n_extant( m ):
	big_start = datetime.now()
	print('Starting to plot n extant.')

	x_mins = m['m_min'].unique()

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)
	
	lw = 1.0
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
		lw = 0.5

	print('All species took {} seconds.'.format((datetime.now() - big_start).total_seconds()))

	sns.plt.xlabel('model time')
	sns.plt.ylabel('number of extant species')
	sns.plt.title('extant species over time, by clade')


def plot_clade_sizes( m ):
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


def plot_clade_spawnrate( m ):
	x_mins = m['m_min'].unique()

	fig = sns.plt.figure()
	ax = fig.add_subplot(111)

	lw = 1.0
	for x_min in x_mins:
		clade = m[m['m_min'] == x_min]

		births = list(map(to_time, clade['birth'].tolist()))
		bd_dist  = np.histogram(births, bins=t_edges)[0]

		ax.plot(t_edges[:-1], bd_dist, linewidth=lw)
		lw = 0.5

	sns.plt.xlabel('model time')
	sns.plt.ylabel('speciation rate (spec appearing per {:.0f} model steps)'.format(t_edges[1]))


# def plot_n_extant( m ):
# 	extants = np.zeros((to_time(len(m.index)), 1))



#plot_dist(all_species['mass'])
#plot_clade_dists(all_species)
plot_clade_largest(all_species)
#plot_clade_largest_wdist(all_species)
#plot_clade_n_extant(all_species)
#plot_clade_sizes(all_species)
plot_clade_spawnrate(all_species)
sns.plt.show()