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
use pbr::ProgressBar;

#[derive(Clone, Debug, RustcEncodable)]
pub struct Species {
	id: usize,          // ID, evolution order
	birth: usize,		  // birthdate
	mass: f64,          // mass (g)
	min_mass: f64,      // mass floor
	death: usize,       // iteration went (or will go) extinct
	parent: usize,      // parent's ID
	niche: usize,		  // niche ID
}


// Compare a generated distribution to the one from the MOM data using a stabilized KS test
fn ks() {
	let mut rdr = csv::Reader::from_file("MOM_data_full.txt").unwrap();

	for row in rdr.decode() {
		let (land, spec_id, gen_id, fam_id, ord_id, mass, code):
			(usize, usize, usize, usize, usize, f64, String) = row.unwrap();
		println!("Species {} has mass {}.", spec_id, mass);
		break;
	}
}

// Generate a simple (multiplicative) random walk.
fn random_walk(x: Vec<f64>, to_step: usize) -> Vec<f64> {
	
	let ret = if to_step > 0 {
		let StandardNormal(r) = random();
		let n = x.last().unwrap() + (r as f64 - 0.5) * 2.0;

		random_walk([&x[..], &[n]].concat(), to_step - 1)
	} else {
		x
	};

	ret
}

#[test]
fn replace_test() {
	let mut extant: Vec<Species> = Vec::new();

	let x_0 = 40.0f64;
	let x_min = 1.8f64;

	let doom_timer = |mass: f64, n: usize| {
		// Extinction parameters
		let beta = 1.0/(n as f64);    // baseline extinction rate
		let rho = 0.025;              // rate of extinction increase

		let p = (10.0 as f64).powf(rho * mass.log10() + beta.log10());

		let x: f64 = random();
		(x.log10() / (1.0 - p).log10()).ceil() as usize
	};

	let ii = 0;

	extant.push(Species{
		id: 1,
		birth: 0,
		mass: x_0,
		min_mass: x_min,
		death: 100,
		parent: 0,
		niche: 0,
	});


	assert_eq!(extant[ii].death, 100);

	let species = extant[ii].clone();

	extant[ii] = Species{
		id: species.id,
		birth: species.birth,
		mass: species.mass,
		min_mass: species.min_mass,
		death: 0,
		parent: species.parent,
		niche: species.niche,
	};

	assert_eq!(extant[ii].death, 0);
}


fn main() {
	let args = clap::App::new("Cladogenesis")
	.version("0.9.27")
	.author("Trevor DiMartino")
	.about("A model of random walk evolution, based on Clauset and Erwin 2008 and modified to include some new capabilities.")
	.args_from_usage("\
		-m --min=[x_min]						'Minimum species mass, g (default 1.8)'
		-i --initial=[x_0]					'Initial species mass, g (default 40.0)'
		-t --tau=[tau]							'Maximum amount of sim time to run (default 250 My)'
		-n --nspecies=[n]						'Estimate of number of species at equilibrium (default 5000)'
		-a --writeall=[write_all]			'Write out birth order, mass, and death time of all species (default true)'
		-r --ratchet=[ratchet]				'Turn on ratcheting capability (default true)'
		-p --r_prob=[r_prob]					'Probability of a ratchet trait evolving (default 0.000005)'
		-s --model_style=[model_style]	'Choose model style (default 1)'
		-b --batches=[batches]				'Set number of runs to make at these settings (default 1)")
	.get_matches();

	let x_min = value_t!(args.value_of("min"), f64).unwrap_or(1.8);
	let x_0 = value_t!(args.value_of("initial"), f64).unwrap_or(40.0);
	let tau = value_t!(args.value_of("tau"), f64).unwrap_or(250.0);     // total simulation time (My)
	let n_0 = value_t!(args.value_of("nspecies"), usize).unwrap_or(5000);
	let write_all = value_t!(args.value_of("writeall"), bool).unwrap_or(true);
	let ratchet = value_t!(args.value_of("ratchet"), bool).unwrap_or(true);
	// Ratchets could be more likely than that.  Maybe 1/260000?
	let r_prob = value_t!(args.value_of("r_prob"), f64).unwrap_or(0.000005);
	let model_style = value_t!(args.value_of("model_style"), usize).unwrap_or(1);
	let batches = value_t!(args.value_of("batches"), usize).unwrap_or(1);

	let write_all = if batches > 1 && write_all {
		println!("WARNING: Writing all in a batch setting will use a lot of disk space.");
		println!("So I'll just turn off writing all species out for you.");
		false
	} else {
		write_all
	};

	//// model_styles:
	// 1. Use "retain" to clean up extant (geodist)
	// 2. Discard extant species upon random selection if they're already dead (slower than #1)
	// 3. Drop a meteor on the extant set halfway through simulation
		// Multiplier for "meteor" (climate change?) bias against seed clade
		// p_0 = 0.5 * (1.0 + m_bias)
		// p_else = 0.5 * (1.0 - m_bias)
		let m_bias = 0.0;  // Equal probabilities for all extant species to die
		// let m_bias = 0.999; // 99.95% extinction for seed clade, 0.05% extinction else
	// 4. Radiation--Increase ratchet probability during the initial radiation phase of the model
			// Until a number of species equal to n has been spawned
		// Probability of ratchet during seed radiation phase
		let r_prob_radiation = 0.5;
	// 5. Promotion--Promote radiation phase for recently ratcheted species
		// Probability of choosing only a recently ratcheted species during promotion phase
		let radiation_preference = 1.00;		// Default: 1.00
		// How long a recently ratcheted species gets preference, in model steps
		let radiation_duration = 25;			// Default: 100


	let mut n = vec![n_0];

	// How likely it is that a ratchet lets the species (and therefore descendants) into a new (latent) nichespace
	let r_niche_prob = 1.0;			// prob that a ratchet gets a niche.  Default: 1.0
	let min_niche_cap = 9;			// min_niche_cap < Min nichespace dimension capacity

	//let r_magnitude: f64 = 0.1;   // x_min increase as result of ratchet (placeholder)

	// Determine how many n_ss to run
	let nu = 1.6;       // mean species lifetime (My)
	// let t_max: usize = 100000;
	let t_max = ((tau / nu) * (n[0] as f64)).ceil() as usize;

	// Cope's Rule parameters
	let c1 = 0.33;      // log-lambda intercept
	let c2 = 1.30;      // log-size intercept
	let delta = 0.04;   // systematic bias

	// Monte Carlo distribution parameters
	let sigma = 0.63;   // variance
	let alpha = 0.30;   // power-law tail


	////////////////////////////////////////////////////////////////////////////
	/////   important functions   //////////////////////////////////////////////

	// Nichespace size determination
	let new_niche_n = |min_d: f64| {
		// x: log, y: log
		// For 69 equines; m_min when niche size = 1: 10^8.205 = 1.603e08
		// let nc0 = 3.7f64;			// Intercept (log)
		// let nc1 = -0.45f64;		// Slope for log/log fit (69 equines)

		// For 7 equines; m_min when niche size = 1: 10^6.3 \approx 2e06
		let nc0 = 3.8f64;
		let nc1 = -0.6f64;

		// To translate from actual niche size to model's n
		let scale_factor: f64 = n_0 as f64 / (10f64.powf(nc0 + x_min.log10() * nc1));

		((10f64.powf(nc0 + min_d.log10() * nc1)) * scale_factor).round() as usize

		// #MOMsim
		77
	};

	// Determines how many more rounds until a species goes extinct
	// by drawing from a geometric distribution
	let doom_timer = |mass: f64, n: usize| {
		// Extinction parameters
		let beta = 1.0/(n_0 as f64);    // baseline extinction rate
		let rho = 0.025;              // rate of extinction increase

		let mut p = 10f64.powf(rho * mass.log10() + beta.log10());
		if p.is_nan() {
			panic!("p is {}.  Inputs: mass = {}, n = {}", p, mass, n);
			// p = 1.0f64;
		} else if p > 1.0f64 {
			p = 1.0f64;
		}

		// (n_0 / n) scales the death length to the timeline of the seed group
		let d = (random::<f64>().log10() / (1.0 - p).log10()).floor() as usize;
		if d > 1000000 {
			panic!("HUGE lifespan: {}.  p_doom = {}", d, p)
		} else if d > 0 {
			d
		} else {
		   1usize
		}
	};

	// For testing, write out a bunch of death times.
	// let doom_mass = 21628432820627857924041177728811008f64;
	// let doom_mass = 26001975398382023000f64;
	// let doom_n = 3usize;
	// let mut writer = Writer::from_file(format!("doom_times_{}_{}.csv", doom_mass, doom_n)).unwrap();
	// for ii in 0..100000 {
	// 	writer.encode(doom_timer(doom_mass, doom_n)).ok().expect("CSV writer error");
	// }

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
		let x_max = 10f64.powi(21);	// Prevent overflow?

		let mut tt: f64 = 0.0;
		while tt < x_min / mass_a  || tt > x_max / mass_a {
			let StandardNormal(r) = random();
			tt = (r * sigma + mu).exp(); //* 
			// ((random::<f64>() * (1.0 - 1.0 / l1) + 1.0 / l1).powf(alpha)) /
			// (random::<f64>().powf(alpha));
		}

		if (mass_a * tt) > x_max {
			println!("New species is about to be larger than {} g: mass_a = {}, tt = {}", x_max, mass_a, tt);
		}

		if (mass_a * tt).is_nan() {
			panic!("New mass is about to be NaN.\nInputs: mass_a = {}, x_min = {}, tt = {}.", mass_a, x_min, tt);
		}

		mass_a * tt
	};

	let choose_anc = |step: usize, extant: &mut Vec<Species>, rec_rat: (f64, usize)| -> Option<(usize, Species)> {
		let mut ca1 = |extant: &mut Vec<Species>| {
			if extant.len() > 0 {
				let ancestor_loc: usize = (random::<f64>() * extant.len() as f64).floor() as usize;
				Some((ancestor_loc, extant[ancestor_loc].clone()))
			} else {
				None
			}
		};

		let mut ca2 = |extant: &mut Vec<Species>| {
			loop {
				let ancestor_loc: usize = (random::<f64>() * extant.len() as f64).floor() as usize;
				let ancestor = extant[ancestor_loc].clone();
				if ancestor.death <= step {
					extant.remove(ancestor_loc);
					return None;
				} else {
					return Some((ancestor_loc, ancestor));
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
					return Some((ancestor_loc, ancestor));

				// If the last ratchet was recent and we chose a species of the promoted clade
				} else if ancestor.min_mass == rec_rat.0 {
					return Some((ancestor_loc, ancestor));

				// If the last ratchet was recent and we chose a species not of that clade
				} else {
					// Check to see if it squeaks by, probabalistically
					if random::<f64>() < radiation_preference {
						continue;
					} else {
						return Some((ancestor_loc, ancestor));
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
			extant.retain(|species| species.death > step);
		};

		let mut cu2 = |ancestor_loc: usize, extant: &mut Vec<Species>| {
			extant.remove(ancestor_loc);
		};

		match model_style {
			2 => cu2(ancestor_loc, extant),
			_ => cu1(extant, step),
		}
	};
	/////////////////////// end "important functions" //////////////////////



	// Everything up to here is just setting up the model, now we run it

	for batch in 0..batches {
		// Set up our data structures
		let mut extant: Vec<Vec<Species>> = Vec::new();

		extant.push(Vec::with_capacity((n_0 as f64 * 1.5).ceil() as usize));
		let doom = doom_timer(x_0, n_0);
		extant[0].push(Species{
			id: 1,
			birth: 0,
			mass: x_0,
			min_mass: x_min,
			death: doom,
			parent: 0,
			niche: 0,
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
				parent: 0,
				niche: 0,
			});
		}


		let start = time::precise_time_ns();

		let mut n_s = 1;
		let mut step = 1;
		let mut recent_ratchet = (x_min, step);  // Used for model_style 5

		// let mut pb = ProgressBar::new(t_max as u64);

		let mut cetacean_flag = true;

		while step <= t_max {
			// pb.inc();

			for nichespace in 0..extant.len() {
				// Downsample our choices from smaller niches so that spec and ext rates match
				if step % (n[0] / n[nichespace]) > 0 {
					continue;
				} else {
					if model_style == 1 {
						// Make sure our extant list is currently accurate
						cleanup(0, &mut extant[nichespace], step);
					}
				}

				let (ancestor_loc, ancestor) = match choose_anc(step, &mut extant[nichespace], recent_ratchet) {
					None => continue,
					Some(s) => s,
				};

				// Sanity check!
				if ancestor.death < step {
					println!("At step {} we picked something that died {}", step, ancestor.death);
				}
				assert!(ancestor.death >= step);

				// Spawn two descendant species
				for _ in 0..2 {
					n_s += 1;

					let mass_d = new_mass(ancestor.mass, ancestor.min_mass);

					// Sanity check!
					assert!(!mass_d.is_nan());

					// See if we've evolved a floor-raising characteristic
					let min_d: f64;
					let niche_d: usize;
					if ratchet {
						// Update ratchet probabilty during initial radiation, if proper model is selected
						let r_eff = if model_style == 4 && n_s < n[nichespace] {
							r_prob_radiation
						} else {
							r_prob
						};

						// #MOMsim
						// if random::<f64>() < r_eff {
						if cetacean_flag && step >= (t_max / 4) && mass_d > 6800.0 && mass_d < 7200.0 {
							cetacean_flag = false;
							// If we get a ratchet, new mass floor is ancestor mass
							min_d = mass_d;
							if model_style == 5 {
								recent_ratchet = (min_d, step + radiation_duration)
							}

							// If we are allowing new nichespaces...
							if random::<f64>() < r_niche_prob {

								// Determine the maximum number of species the new niche can sustain
								// Decreases linearly with respsect to log(x_min)
								let n_new: usize = new_niche_n(min_d);
								// We can't have niches with space for only one species.
								if n_new > min_niche_cap {
									niche_d = extant.len();
									n.push(n_new);
									extant.push(Vec::with_capacity((n_new as f64 * 1.5).ceil() as usize));
									println!("\nCreating new niche dimension, {}, with capacity {} and m_min {}", extant.len(), n_new, min_d);
								} else {
									niche_d = ancestor.niche;
									println!("\nCould have created a new niche dimension, but m_min = {} is too large.", min_d);
								}
							} else {
								niche_d = ancestor.niche;
							}
						} else {
							// Else, the min remains the min of the ancestor
							min_d = ancestor.min_mass;
							niche_d = ancestor.niche;
						}
					} else {
						min_d = ancestor.min_mass;
						niche_d = ancestor.niche;
					}

					let doom = doom_timer(mass_d, n[niche_d]) + step;
					if (doom - step) > 1000000 {
						println!("Currently step = {}.\nSpecies of mass {}, in clade with n = {}, won't die until {}", step, mass_d, n[niche_d], doom);
						for _ in 0..10 {
							println!("Choosing again gives lifespan of {}.", doom_timer(mass_d, n[niche_d]));
						}
						assert!(doom < 1000000);
					}

					extant[niche_d].push(Species{
						id: n_s, 
						birth: step,
						mass: mass_d, 
						min_mass: min_d, 
						death: doom,
						parent: ancestor.id,
						niche: niche_d,
					});

					if write_all {
						all_species.push(Species{
							id: n_s,
							birth: step,
							mass: mass_d,
							min_mass: min_d,
							death: doom,
							parent: ancestor.id,
							niche: niche_d,
						});
					}
				}
				
				// Set the ancestor species to die in cleanup.
				extant[nichespace][ancestor_loc] = Species{
					id: ancestor.id,
					birth: ancestor.birth,
					mass: ancestor.mass,
					min_mass: ancestor.min_mass,
					death: step,
					parent: ancestor.parent,
					niche: ancestor.niche,
				};

				if write_all {
					// Record the fact that this species died here.
					all_species[ancestor.id - 1] = Species{
						id: ancestor.id,
						birth: ancestor.birth,
						mass: ancestor.mass,
						min_mass: ancestor.min_mass,
						death: step,
						parent: ancestor.parent,
						niche: ancestor.niche,
					};
				}

				// If we are told to, check if we're halfway through the sim, and drop a meteor
				if model_style == 3 {
					if step == t_max / 2 {
						// Drop a meteor!
						let mut killed = 0;
						for ii in 0..extant[nichespace].len() {
							let species = extant[nichespace][ii].clone();

							let p_dead = if species.min_mass == x_min {
								// It's part of the seed clade
								0.5 * (1.0 + m_bias)
							} else {
								0.5 * (1.0 - m_bias)
							};

							if random::<f64>() < p_dead {
								killed += 1;

								extant[nichespace][ii] = Species{
									id: species.id,
									birth: species.birth,
									mass: species.mass,
									min_mass: species.min_mass,
									death: step,
									parent: species.parent,
									niche: species.niche,
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
											parent: spec.parent,
											niche: spec.niche,
										};
								}
							}
						}

						println!("A meteor fell and killed off {} out of {} species.", killed, extant[nichespace].len());
					}
				}

				// Clean up extant, depending on the model 
				cleanup(ancestor_loc, &mut extant[nichespace], step);
			}  // End for nichespaces

			step += 1;
		}  // End while < t_max

		if model_style == 2 {
			// Go through the list once and clean it up by removing everything that's dead
			for nichespace in 0..extant.len() {
				extant[nichespace].retain(|species| species.death >= step);
			}
		}

		// End timing
		let end = time::precise_time_ns();
		// pb.finish_print("Complete\n");
		for nichespace in 0..extant.len() {
			println!("Niche {} ended with a capacity of {} and {} extant species.", nichespace, n[nichespace], extant[nichespace].len())
		}
		println!("Ran Model {} for {} species in {} seconds.", model_style, n_s, (end - start) as f64 / 1000000000.0);

		// Print out our final set of extant species
		let path = if batches == 1 {
			format!("extant_m{}_{}_{}_{}.csv", model_style, x_min, x_0, n_0)
		} else {
			format!("extant_m{}_{}_{}_{}_b{}.csv", model_style, x_min, x_0, n_0, batch)
		};
		let mut writer = Writer::from_file(path).unwrap();
		for ns in extant.into_iter() {
			for s in ns.into_iter() {
				writer.encode(s).ok().expect("CSV writer error");
			}
		}

		if write_all {
			let start = time::precise_time_ns();

			// All species from the HashMap
			let path = format!("all_m{}_{}_{}_{}.csv", model_style, x_min, x_0, n_0);
			let mut writer = Writer::from_file(path).unwrap();
			for s in all_species.into_iter() {
				writer.encode(s).ok().expect("CSV writer error");
			}

			let end = time::precise_time_ns();
			println!("Saved all in {} seconds.", (end - start) as f64 / 1000000000.0);
		}
	}

}