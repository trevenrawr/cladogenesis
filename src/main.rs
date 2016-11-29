extern crate rand;
extern crate time;
extern crate csv;
extern crate rustc_serialize;
#[macro_use]
extern crate clap;

use rand::random;
use std::collections::{HashMap, BinaryHeap};
use std::collections::hash_map::Entry;
use rand::distributions::normal::StandardNormal;
use csv::Writer;
use std::cmp::Ordering;

#[derive(Clone, Debug, RustcEncodable)]
struct Species {
	id: i64,          // ID, evolution order
	birth: f64,			// Birthdate
	mass: f64,          // mass (g)
	min_mass: f64,      // mass floor
	death: f64,       // iteration went (or will go) extinct
	parent: i64,      // parent's ID
}

impl PartialEq for Species {
	fn eq(&self, other: &Species) -> bool {
		self.id == other.id
	}
}

impl Eq for Species {}

impl PartialOrd for Species {
	fn partial_cmp(&self, other: &Species) -> Option<Ordering> {
		Some(self.cmp(other))
	}
}

impl Ord for Species {
	fn cmp(&self, other: &Species) -> Ordering {
		self.id.cmp(&other.id)
	}
}

#[derive(Clone, Debug, RustcEncodable, PartialEq, Eq, PartialOrd, Ord)]
enum Flavor {
	Speciation,
	Extinction,
}

#[derive(Clone, Debug, RustcEncodable)]
struct Event {
	time: f64,
	flavor: Flavor,
	species: Species,
}

impl PartialEq for Event {
	fn eq(&self, other: &Event) -> bool {
		self.time == other.time && 
			self.flavor == other.flavor &&
			self.species == other.species
	}
}

impl Eq for Event {}

impl PartialOrd for Event {
	fn partial_cmp(&self, other: &Event) -> Option<Ordering> {
		Some(self.cmp(other))
	}
}

impl Ord for Event {
	fn cmp(&self, other: &Event) -> Ordering {
		if self.time == other.time {
			// Check event flavors to break ties
			if self.flavor == other.flavor {
				// And then order by species (by id, arbitrary)
				return self.species.cmp(&other.species)
			}

			return self.flavor.cmp(&other.flavor)
		}

		// Note: We return opposite time ordering because we want smallest 
		//       times first in the BinaryHeap (priority queue).
		if self.time < other.time {
			Ordering::Greater
		// We already covered for equality.
		} else {
			Ordering::Less
		}
	}
}

#[derive(Clone)]
struct Params {
	x_min: f64,
	x_0: f64,
	n: usize,
	ratchet: bool,
	r_prob: f64,
	max_descendants: i64,

	// Timing parameters
	mean_lifespan: f64,
	t_max: f64,

	// Common rates
	spec_lambda: f64,
	ext_lambda: f64,
	rho: f64,
}


// Determines how many more rounds until a species goes extinct
// by drawing from a geometric distribution
fn doom_timer(ext_lambda:f64, rho: f64, mass: f64) -> f64 {
	let lambda = (10.0f64).powf(ext_lambda.log10()  - rho * mass.log10());
	let x: f64 = random();

	(1.0f64 - x).ln() / -lambda
}

// Determines how many more rounds until a species cladogenesizes
fn spec_timer(spec_lambda: f64) -> f64 {
	let x: f64 = random();

	(1.0f64 - x).ln() / -spec_lambda
}

// Draw new spawn times until one of them occurs after death.
fn to_spawn(t: f64, n_ext: f64, species: Species, spec_timer: fn(f64) -> f64, params: Params) -> Vec<Event> {
	let mut spec_time = spec_timer(params.spec_lambda) + t;
	let mut spec_events: Vec<Event> = Vec::new();
	let mut spawned: i64 = 0;

	while spec_time < species.death && spawned < params.max_descendants {
		spawned += 1;

		spec_events.push(Event{
			time: spec_time,
			flavor: Flavor::Speciation,
			species: species.clone(),
		});

		spec_time = spec_timer(params.spec_lambda * (params.n as f64 / (params.n as f64 + n_ext)));
	}

	spec_events
}

// Takes ancestor mass and returns a new descendant species mass
fn new_mass(mass_a: f64, x_min: f64) -> f64 {
	// Cope's Rule parameters
	let c1 = 0.33;			// log-lambda intercept
	let c2 = 1.30;			// log-size intercept
	let delta = 0.04;		// systematic bias

	// Monte Carlo distribution parameters
	let sigma = 0.63;		// variance
	let alpha = 0.30;		// power-law tail

	// Determines the groth factor, mu, determined by Cope's rule
	// and strengthened for smaller species
	let mu: f64;
	if mass_a.log10() < c2 {
		mu = (-c1 / c2) * mass_a.log10() + c1 + delta;
	} else {
		mu = delta;
	}

	// Draw Monte Carlo factor and enforce min size
	let l1 = mass_a / x_min;
	let mut tt: f64 = 0.0;
	while tt < 1.0 / l1 {
		let StandardNormal(r) = random();
		tt = (r * sigma + mu).exp(); //* 
		//((random::<f64>() * (1.0 - 1.0 / l1) + 1.0 / l1).powf(alpha)) /
		//(random::<f64>().powf(alpha));
	}

	mass_a * tt
}

// Takes a bunch of information and updates extant and event lists
fn speciate(events: &mut BinaryHeap<Event>,
				extant: &mut HashMap<i64, Species>,
				all_species: &mut Vec<Species>, 
				event: &Event,
				n_s: &mut i64,
				params: Params
			) {

	*n_s = *n_s + 1;

	let t = event.time;

	let mass_d = new_mass(event.species.mass, event.species.min_mass);

	// See if we've evolved a floor-raising characteristic
	let min_d: f64;
	if params.ratchet {
		if random::<f64>() < params.r_prob {
			// If we get a ratchet, new mass floor is ancestor mass
			min_d = mass_d;
		} else {
			// Else, the min remains the min of the ancestor
			min_d = event.species.min_mass;
		}
	} else {
		min_d = event.species.min_mass;
	}

	let doom = doom_timer(params.ext_lambda, params.rho, mass_d) + t;
	//println!("Species {} will die at {}", n_s, doom);

	// Take note of our new species
	let descendant = Species{
			id: *n_s,
			birth: t,
			mass: mass_d,
			min_mass: min_d,
			death: doom.clone(),
			parent: event.species.id
		};

	extant.insert(descendant.id, descendant.clone());
	all_species.push(descendant.clone());

	// Put their extinction event in the list
	events.push(Event{
		time: doom,
		flavor: Flavor::Extinction,
		species: descendant.clone(),
	});

	let n_ext = extant.len() as f64;

	let spec_events = to_spawn(t, n_ext, descendant, spec_timer, params);
	for spec_event in spec_events {
		events.push(spec_event.clone());
	}
}

fn extinate(extant: &mut HashMap<i64, Species>, 
				event: &Event) {

	extant.remove(&event.species.id);
}


fn main() {
	let args = clap::App::new("Cladogenesis")
	.version("0.9.27")
	.author("Trevor DiMartino")
	.about("A model of random walk evolution, based on Clauset and Erwin 2008 and modified to include some new capabilities.")
	.args_from_usage("\
		-m --min=[x_min]            'Minimum species mass, g (default 1.8)'
		-i --initial=[x_0]          'Initial species mass, g (default 40.0)'
		-n --nspecies=[n]           'Estimate of number of species at equilibrium (default 5000)'
		-r --ratchet=[ratchet]      'Turn on ratcheting capability (default true)'
		-p --r_prob=[r_prob]        'Probability of a ratchet trait evolving (default 0.0001)'
		-d --max_d=[max_d]          'Maximum number of descendants a species can have (default 2)' ")
	.get_matches();

	// Timing parameters
	let mean_lifespan = 2.012f64;				// Average species lifetime (My)

	// Common rates
	let spec_lambda = 2.0/mean_lifespan;		// Baseline speciation rate
	let ext_lambda = 1.0/mean_lifespan;		// Baseline extinction rate

	let params = Params{
		x_min:				value_t!(args.value_of("min"), f64).unwrap_or(1.8),
		x_0:					value_t!(args.value_of("initial"), f64).unwrap_or(40.0),
		n:						value_t!(args.value_of("nspecies"), usize).unwrap_or(5000),
		ratchet:				value_t!(args.value_of("ratchet"), bool).unwrap_or(true),
		r_prob:				value_t!(args.value_of("r_prob"), f64).unwrap_or(0.0001),
		max_descendants:	value_t!(args.value_of("max_d"), i64).unwrap_or(2),
		mean_lifespan: mean_lifespan,
		t_max: 10.0f64,
		spec_lambda: spec_lambda,
		ext_lambda: ext_lambda,
		rho: 0.025f64,
	};

	//write_distributions(spec_lambda, ext_lambda, params.rho);

	// Set up our data structures
	// extant a running set of who's alive
	// timer = BinaryHeap<usize> keeps the priority queue of death dates to come next
	let mut extant: HashMap<i64, Species> = HashMap::with_capacity((params.n as f64 * 1.5).ceil() as usize);
	let mut events = BinaryHeap::with_capacity((params.n * 2) as usize);

	// Set up initial extinction event
	let doom = doom_timer(params.ext_lambda, params.rho, params.x_0);

	let og = Species{
			id: 1, 
			birth: 0.0,
			mass: params.x_0,
			min_mass: params.x_min,
			death: doom,
			parent: 0,
		};

	extant.insert(og.id, og.clone());

	events.push(Event{
		time: doom.clone(),
		flavor: Flavor::Extinction,
		species: og.clone(),
	});

	let spec_events = to_spawn(0.0, 1.0, og.clone(), spec_timer, params.clone());
	println!("Upcoming speciations: {:?}", spec_events);
	for spec_event in spec_events {
		events.push(spec_event.clone());
	}


	let mut all_species = Vec::with_capacity(2 * params.t_max as usize + 1);
	all_species.push(og);


	// Start timing (for reflection)
	let start = time::precise_time_ns();

	// We start with our total number of species to be 1; the supreme event.species
	let mut n_s = 1;
	let mut t = 0.0;
	while t < params.t_max && !events.is_empty() {
		// Get the time from our timer
		let event = events.pop().unwrap();

		t = event.time;

		match event.flavor {
			Flavor::Speciation => speciate(
					&mut events,
					&mut extant,
					&mut all_species,
					&event,
					&mut n_s, 
					params.clone()
				),
			Flavor::Extinction => extinate(&mut extant, &event),
		}
	}

	// End timing
	let end = time::precise_time_ns();
	println!("Ran model for {} species in {} seconds.", n_s, (end - start) as f64 / 1000000000.0);

	// Print out our final set of extant species
	let path = format!("extant_ecm1_{}_{}_{}.csv", params.x_min, params.x_0, params.n);
	let mut writer = Writer::from_file(path).unwrap();
	for s in extant.into_iter() {
		writer.encode(s).ok().expect("CSV writer error");
	}


	let path = format!("all_ecm1_{}_{}_{}.csv", params.x_min, params.x_0, params.n);
	let mut writer = Writer::from_file(path).unwrap();
	for s in all_species.into_iter() {
		writer.encode(s).ok().expect("CSV writer error");
	}
}

fn write_distributions(spec_lambda: f64, ext_lambda: f64, rho: f64) {
	let mut deaths = Vec::with_capacity(10000);
	for _ in 0..10000 {
		deaths.push(doom_timer(ext_lambda, rho, 40.0));
	}

	let mut births = Vec::with_capacity(10000);
	for _ in 0..10000 {
		births.push(spec_timer(spec_lambda));
	}

	let mut writer = Writer::from_file("ec_death_dist.csv").unwrap();
	for s in deaths.into_iter() {
		writer.encode(s).ok().expect("CSV writer error");
	}

	let mut writer = Writer::from_file("ec_birth_dist.csv").unwrap();
	for s in births.into_iter() {
		writer.encode(s).ok().expect("CSV writer error");
	}

	println!("Wrote speciation and extinction timer distributions.");
}
