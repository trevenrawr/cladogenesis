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

#[derive(Clone, Debug, RustcEncodable)]
pub struct Species {
  id: i64,          // ID, evolution order, "birthdate"
  mass: f64,          // mass (g)
  min_mass: f64,      // mass floor
  death: i64,       // iteration went (or will go) extinct
  parent: i64,      // parent's ID
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
      -a --writeall=[write_all]   'Write out birth order, mass, and death time of all species (default true)'
      -r --ratchet=[ratchet]      'Turn on ratcheting capability (default true)'
      -p --r_prob=[r_prob]        'Probability of a ratchet trait evolving (default 0.0001)' ")
    .get_matches();

  let x_min = value_t!(args.value_of("min"), f64).unwrap_or(1.8);
  let x_0 = value_t!(args.value_of("initial"), f64).unwrap_or(40.0);
  let n = value_t!(args.value_of("nspecies"), usize).unwrap_or(5000);
  let write_all = value_t!(args.value_of("writeall"), bool).unwrap_or(true);
  let ratchet = value_t!(args.value_of("ratchet"), bool).unwrap_or(true);
  let r_prob = value_t!(args.value_of("r_prob"), f64).unwrap_or(0.0001);

  //let r_magnitude: f64 = 0.1;   // x_min increase as result of ratchet (placeholder)

  // Determine how many n_ss to run
  let nu = 1.6;       // mean species lifetime (My)
  let tau = 500.0;     // total simulation time (My)
  let t_max = ((tau / nu) * (n as f64)).ceil() as i64;
  println!("Set to run for {} iterations", t_max);

  // Cope's Rule parameters
  let c1 = 0.33;      // log-lambda intercept
  let c2 = 1.30;      // log-size intercept
  let delta = 0.04;   // systematic bias

  // Monte Carlo distribution parameters
  let sigma = 0.63;   // variance
  let alpha = 0.30;   // power-law tail

  // Extinction parameters
  let beta = 1.0/(n as f64);    // baseline extinction rate
  let rho = 0.025;               // rate of extinction increase

  // Determines how many more rounds until a species goes extinct
  // by drawing from a geometric distribution
  let doom_timer = |mass: f64| {
    let p = (10.0 as f64).powf(rho * mass.log10() + beta.log10());

    let x: f64 = random();
    (x.log10() / (1.0 - p).log10()).ceil() as i64
  };

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
      tt = (r * sigma + mu).exp() * 
        ((random::<f64>() * (1.0 - 1.0 / l1) + 1.0 / l1).powf(alpha)) /
        (random::<f64>().powf(alpha));
    }

    mass_a * tt
  };

  // Set up our data structures
  // extant = HashMap<extinction "date", Vec<(id, mass, mass floor)>>
  // timer = BinaryHeap<usize> keeps the priority queue of death dates to come next
  let mut extant: HashMap<i64, Vec<(i64, f64, f64)>> = HashMap::with_capacity((n as f64 * 1.5).ceil() as usize);
  let mut timer = BinaryHeap::with_capacity((n as f64 * 1.5) as usize);

  let doom = doom_timer(x_0);
  //println!("First doom timer is {}", doom);

  timer.push(-doom);
  {
    let descendant = extant.entry(doom).or_insert(Vec::new());
    descendant.push((1, x_0, x_min));
  }

  let mut all_species: Vec<Species> = Vec::with_capacity(0);
  if write_all {
    all_species = Vec::with_capacity(2 * t_max as usize + 1);
    all_species.push(Species{id: 1, mass: x_0, min_mass: x_min, death: doom, parent: 0});
  }


  // Start timing (for reflection)
  let start = time::precise_time_ns();

  // We start with our total number of species to be 1; the supreme ancestor, above
  let mut n_s = 1;
  while n_s < t_max && !timer.is_empty() {
    // Get the (negative) time from our timer
    let t = -timer.pop().unwrap();
    let r = -t;
    //println!("Now serving time t = {}", t);
    while !timer.is_empty() && &r == timer.peek().unwrap() {
      // Pop any duplicate values
      timer.pop();
    }

    // Pull our ancestors vector from the HashMap
    let ancestors = if let Entry::Occupied(o) = extant.entry(t) {
      o.remove()
    } else {
      panic!("There's no entry in extant for {}", t);
    };

    for (id_a, mass_a, min_a) in ancestors {
      let p_spec = 1.0 - (0.5 / n as f64) * extant.len() as f64;

      // Check to see if this extinction spurred a speciation event
      if random::<f64>() < p_spec {

        // If we speciate, we speciate twice
        for _ in 0..2 {
          // We're now working on creating this new species
          n_s += 1;

          if n_s % 10000 == 0 { println!("Running species {}", n_s); }

          let mass_d = new_mass(mass_a, min_a);
    
          // See if we've evolved a floor-raising characteristic
          let min_d: f64;
          if ratchet {
            if random::<f64>() < r_prob {
              // If we get a ratchet, new mass floor is ancestor mass
              min_d = mass_a;
            } else {
              // Else, the min remains the min of the ancestor
              min_d = min_a;
            }
          } else {
            min_d = min_a;
          }
    
          let doom = doom_timer(mass_d) + t;
          //println!("Species {} will die at {}", n_s, doom);
          
          timer.push(-doom);
          let descendant = extant.entry(doom).or_insert(Vec::new());
          descendant.push((n_s, mass_d, min_d));
    
          if write_all {
            all_species.push(Species{
              id: n_s,
              mass: mass_d,
              min_mass: min_d,
              death: doom,
              parent: id_a
            });
          }
        }
      }
    }
  }

  // End timing
  let end = time::precise_time_ns();

  // Print out our final set of extant species
  let path = format!("extant_{}_{}_{}.csv", x_min, x_0, n);
  let mut writer = Writer::from_file(path).unwrap();
  for (_, s) in extant.into_iter() {
    for ss in s {
      writer.encode(ss).ok().expect("CSV writer error");
    }
  }

  if write_all {
    let path = format!("all_{}_{}_{}.csv", x_min, x_0, n);
    let mut writer = Writer::from_file(path).unwrap();
    for s in all_species.into_iter() {
      writer.encode(s).ok().expect("CSV writer error");
    }
  }

  println!("Ran in {} seconds.", (end - start) / 1000000000);
}
