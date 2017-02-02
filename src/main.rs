extern crate rand;
extern crate time;
extern crate csv;
extern crate rustc_serialize;
extern crate pbr;
#[macro_use]
extern crate clap;

use rand::random;
use rand::distributions::normal::StandardNormal;
use csv::Writer;
//use pbr::ProgressBar;

#[derive(Clone, Debug, RustcEncodable)]
pub struct Species {
	id: usize,          // ID, evolution order
	birth: usize,		  // birthdate
	mass: f64,          // mass (g)
	min_mass: f64,      // mass floor
	death: usize,       // iteration went (or will go) extinct
	parent: usize,      // parent's ID
}


fn main() {
	let args = clap::App::new("Cladogenesis")
	.version("0.9.27")
	.author("Trevor DiMartino")
	.about("A model of random walk evolution, based on Clauset and Erwin 2008 and modified to include some new capabilities.")
	.args_from_usage("\
		-m --min=[x_min]						'Minimum species mass, g (default 1.8)'
		-i --initial=[x_0]					'Initial species mass, g (default 40.0)'
		-n --nspecies=[n]						'Estimate of number of species at equilibrium (default 5000)'
		-a --writeall=[write_all]			'Write out birth order, mass, and death time of all species (default true)'
		-r --ratchet=[ratchet]				'Turn on ratcheting capability (default true)'
		-p --r_prob=[r_prob]					'Probability of a ratchet trait evolving (default 0.0001)'
		-s --model_style=[model_style]	'Choose model style (default 1)' ")
	.get_matches();

	let x_min = value_t!(args.value_of("min"), f64).unwrap_or(1.8);
	let x_0 = value_t!(args.value_of("initial"), f64).unwrap_or(40.0);
	let n = value_t!(args.value_of("nspecies"), usize).unwrap_or(5000);
	let write_all = value_t!(args.value_of("writeall"), bool).unwrap_or(true);
	let ratchet = value_t!(args.value_of("ratchet"), bool).unwrap_or(true);
	let r_prob = value_t!(args.value_of("r_prob"), f64).unwrap_or(0.0001);
	let model_style = value_t!(args.value_of("model_style"), usize).unwrap_or(1);

	//// model_styles:
	// 1. Use "retain" to clean up extant
	// 2. Discard extant species upon random selection if they're already dead
	// 3. Drop a meteor on the extant set halfway through simulation
		// Multiplier for "meteor" (climate change?) bias against seed clade
		let m_bias = 0.9;
	// 4. Increase ratchet probability during the initial radiation phase of the model
		// Probability of ratchet during seed radiation phase
		let r_prob_radiation = 0.25;
	// 5. Promote radiation phase for recently ratcheted species
		// Probability of not accepting a non-recently ratcheted species during promotion phase
		let radiation_preference = 1.00;		// Default: 1.00
		// How long a recently ratcheted species gets preference, in model steps
		let radiation_duration = 100;			// Default: 100
	// 6. Create new niche-spaces for some ratchets
		// How likely it is that a ratchet lets the species (and therefore descendants) into a new (latent) nichespace
		let r_niche_prob = 0.1;					// Default: 0.1?


	//let r_magnitude: f64 = 0.1;   // x_min increase as result of ratchet (placeholder)

	// Determine how many n_ss to run
	let nu = 1.6;       // mean species lifetime (My)
	let tau = 500.0;     // total simulation time (My)
	// let t_max: usize = 25000;
	let t_max = ((tau / nu) * (n as f64)).ceil() as usize;

	// Cope's Rule parameters
	let c1 = 0.33;      // log-lambda intercept
	let c2 = 1.30;      // log-size intercept
	let delta = 0.04;   // systematic bias

	// Monte Carlo distribution parameters
	let sigma = 0.63;   // variance
	let alpha = 0.30;   // power-law tail

	// Extinction parameters
	let beta = 1.0/(n as f64);    // baseline extinction rate
	let rho = 0.025;              // rate of extinction increase

	// Determines how many more rounds until a species goes extinct
	// by drawing from a geometric distribution
	let doom_timer = |mass: f64| {
		let p = (10.0 as f64).powf(rho * mass.log10() + beta.log10());

		let x: f64 = random();
		(x.log10() / (1.0 - p).log10()).ceil() as usize
	};

	// Set up our data structures
	// extant: Vec<(id, mass, mass floor, extinction "date")>
	let mut extant = Vec::with_capacity((n as f64 * 1.5).ceil() as usize);
	let doom = doom_timer(x_0);
	extant.push(Species{
		id: 1,
		birth: 0,
		mass: x_0,
		min_mass: x_min,
		death: doom,
		parent: 0,
	});

	let mut all_species: Vec<Species> = Vec::with_capacity(0);
	if write_all {
		all_species = Vec::with_capacity(2 * t_max + 1);
		all_species.push(Species{
			id: 1, 
			birth: 0,
			mass: x_0, 
			min_mass: x_min, 
			death: doom, 
			parent: 0
		});
	}

	////////////////////////////////////////////////////////////////////////////
	/////   important functions   //////////////////////////////////////////////

	// Takes ancestor mass and returns a new descendant species mass
	let new_mass = |mass_a: f64, x_min: f64| -> f64 {
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
			// ((random::<f64>() * (1.0 - 1.0 / l1) + 1.0 / l1).powf(alpha)) /
			// (random::<f64>().powf(alpha));
		}

		mass_a * tt
	};

	let choose_anc = |step: usize, extant: &mut Vec<Species>, rec_rat: (f64, usize)| -> (usize, Species) {
		let mut ca1 = |extant: &mut Vec<Species>| {
			let ancestor_loc: usize = (random::<f64>() * extant.len() as f64).floor() as usize;
			(ancestor_loc, extant[ancestor_loc].clone())
		};

		let mut ca2 = |extant: &mut Vec<Species>| {
			loop {
				let ancestor_loc: usize = (random::<f64>() * extant.len() as f64).floor() as usize;
				let ancestor = extant[ancestor_loc].clone();
				if ancestor.death <= step {
					extant.remove(ancestor_loc);
				} else {
					return (ancestor_loc, ancestor);
				}
			}
		};

		// Use for model_style 5
		let mut ca5 = |extant: &mut Vec<Species>| {
			let mut choices = 0;
			loop {
				let ancestor_loc: usize = (random::<f64>() * extant.len() as f64).floor() as usize;
				let ancestor = extant[ancestor_loc].clone();

				// If the last ratchet wasn't recent, take whoever we chose
				if step > rec_rat.1 {
					return (ancestor_loc, ancestor);

				// If the last ratchet was recent and we chose a species of the promoted clade
				} else if ancestor.min_mass == rec_rat.0 {
					return (ancestor_loc, ancestor);

				// If the last ratchet was recent and we chose a species not of that clade
				} else {
					// Check to see if it squeaks by, probabalistically
					if random::<f64>() < radiation_preference {
						continue;
					} else {
						return (ancestor_loc, ancestor);
					}
				}
				choices += 1;

				if choices > extant.len() * 5 {
					println!("Can't find a suitable ancestor.");
					println!("rr_expire = {}", rec_rat.1);
					println!("rr_mass = {}", rec_rat.0);
					panic!("Throwing in the towel; ancestor must not exist.");
				}
			}
		};

		match model_style {
			2 => ca2(extant),
			5 => ca5(extant),
			_ => ca1(extant),
		}
	};

	let cleanup = |ancestor_loc: usize, extant: &mut Vec<Species>, step: usize| {
		let mut cu1 = |extant: &mut Vec<Species>, step: usize| {
			// Retain species whose doom timer is later than when we are now
			// ... that is to say, kill off the now-dead ones
			extant.retain(|species| species.death >= step);
		};

		let mut cu2 = |ancestor_loc: usize, extant: &mut Vec<Species>| {
			extant.remove(ancestor_loc);
		};

		match model_style {
			2 => cu2(ancestor_loc, extant),
			_ => cu1(extant, step),
		}
	};
	////////////////////////////////////////////////////////////////////////////


	let start = time::precise_time_ns();

	let mut n_s = 1;
	let mut step = 1;
	let mut recent_ratchet = (x_min, step);  // Used for model_style 5

	//let mut pb = ProgressBar::new(t_max as u64);

	while step <= t_max {
		//pb.inc();
		let (ancestor_loc, ancestor) = choose_anc(step, &mut extant, recent_ratchet);

		// Spawn two descendant species
		for _ in 0..2 {
			n_s += 1;

			let mass_d = new_mass(ancestor.mass, ancestor.min_mass);

			// See if we've evolved a floor-raising characteristic
			let min_d: f64;
			if ratchet {
				// Update ratchet probabilty during initial radiation, if proper model is selected
				let r_eff = if model_style == 4 && n_s < n {
					r_prob_radiation
				} else {
					r_prob
				};

				if random::<f64>() < r_eff {
					// If we get a ratchet, new mass floor is ancestor mass
					min_d = mass_d;
					if model_style == 5 {
						recent_ratchet = (min_d, step + radiation_duration)
					}

					// If we are allowing new nichespaces...
					if model_style == 6 {
						if random::<f64>() < r_niche_prob {
							
						}
					}
				} else {
					// Else, the min remains the min of the ancestor
					min_d = ancestor.min_mass;
				}
			} else {
				min_d = ancestor.min_mass;
			}

			let doom = doom_timer(mass_d) + step;

			extant.push(Species{
				id: n_s, 
				birth: step,
				mass: mass_d, 
				min_mass: min_d, 
				death: doom,
				parent: ancestor.id,
			});

			if write_all {
				all_species.push(Species{
					id: n_s,
					birth: step,
					mass: mass_d,
					min_mass: min_d,
					death: doom,
					parent: ancestor.id,
				});
			}
		}
		
		// Set the ancestor species to die in cleanup.
		extant[ancestor_loc] = Species{
			id: ancestor.id,
			birth: ancestor.birth,
			mass: ancestor.mass,
			min_mass: ancestor.min_mass,
			death: step,
			parent: ancestor.parent,
		};

		// If we are told to, check if we're halfway through the sim, and drop a meteor
		if model_style == 3 {
			if step == t_max / 2 {
				// Drop a meteor!
				let mut killed = 0;
				for ii in 0..extant.len() {
					let species = extant[ii].clone();

					let p_dead = if species.min_mass == x_min {
						// It's part of the seed clade
						0.5 * (1.0 + m_bias)
					} else {
						0.5 * (1.0 - m_bias)
					};

					if random::<f64>() < p_dead {
						extant[ii] = Species{
							id: species.id,
							birth: species.birth,
							mass: species.mass,
							min_mass: species.min_mass,
							death: step,
							parent: species.parent,
						};

						if write_all {
							// Log, in all_species, that they died now due to meteor
							let spec = all_species[species.id - 1].clone();

							// Sanity check!
							assert_eq!(species.mass, spec.mass);

							all_species[species.id - 1] = Species{
									id: spec.id,
									birth: spec.birth,
									mass: spec.mass,
									min_mass: spec.min_mass,
									death: step,
									parent: spec.parent
								};
						}

						killed += 1;
					}
				}

				println!("A meteor fell and killed off {} out of {} species.", killed, extant.len());
			}
		}

		// Clean up extant, depending on the model 
		cleanup(ancestor_loc, &mut extant, step);

		step += 1;
	}

	if model_style == 2 {
		// Go through the list once and clean it up by removing everything that's dead
		extant.retain(|species| species.death >= step);
	}

	// End timing
	let end = time::precise_time_ns();
	//pb.finish_print("Complete");
	println!("\nRan Model {} for {} species in {} seconds.", model_style, n_s, (end - start) as f64 / 1000000000.0);

	// Print out our final set of extant species
	let path = format!("extant_m{}_{}_{}_{}.csv", model_style, x_min, x_0, n);
	let mut writer = Writer::from_file(path).unwrap();
	for s in extant.into_iter() {
		writer.encode(s).ok().expect("CSV writer error");
	}

	if write_all {
		let start = time::precise_time_ns();

		// All species from the HashMap
		let path = format!("all_m{}_{}_{}_{}.csv", model_style, x_min, x_0, n);
		let mut writer = Writer::from_file(path).unwrap();
		for s in all_species.into_iter() {
			writer.encode(s).ok().expect("CSV writer error");
		}

		let end = time::precise_time_ns();
		println!("Saved all in {} seconds.", (end - start) as f64 / 1000000000.0);
	}

}