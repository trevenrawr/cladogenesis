extern crate rand;
extern crate time;
extern crate csv;
extern crate rustc_serialize;
#[macro_use]
extern crate clap;

use rand::random;
use rand::distributions::normal::StandardNormal;
use csv::Writer;

#[derive(Clone, Debug, RustcEncodable)]
pub struct Species {
  id: usize,          // ID, evolution order, "birthdate"
  mass: f64,          // mass (g)
  min_mass: f64,      // mass floor
  death: usize,       // iteration went (or will go) extinct
  parent: usize,      // parent's ID
}


fn main() {
  let args = clap::App::new("Cladogenesis")
    .version("0.1")
    .author("Trevor DiMartino")
    .about("A model of random walk evolution, based on Clauset and Erwin 2008 and modified to include some new capabilities.")
    .args_from_usage("\
      -m --min=[x_min]            'Minimum species mass, g (default 1.8)'
      -i --initial=[x_0]          'Initial species mass, g (default 40.0)'
      -n --nspecies=[n]           'Estimate of number of species at equilibrium (default 5000)'
      -a --writeall=[write_all]   'Write out birth order, mass, and death time of all species (default true)'
      -r --ratchet=[ratchet]      'Turn on ratcheting capability (default false)' ")
    .get_matches();

  let mut x_min = value_t!(args.value_of("min"), f64).unwrap_or(1.8);
  let x_0 = value_t!(args.value_of("initial"), f64).unwrap_or(40.0);
  let n = value_t!(args.value_of("nspecies"), usize).unwrap_or(5000);
  let write_all = value_t!(args.value_of("writeall"), bool).unwrap_or(true);
  let ratchet = value_t!(args.value_of("ratchet"), bool).unwrap_or(false);

  // Ratchet mechanisms
  let r_prob: f64 = 0.0002;   // Probability of a ratchet occurrance (placeholder)
  let r_magnitude: f64 = 0.1;   // 10% increase as result of ratchet (placeholder)

  // Determine how many steps to run
  let nu = 1.6;       // mean species lifetime (My)
  let tau = 60.0;     // total simulation time (My)
  let t_max: usize = ((tau / nu) * (n as f64)).ceil() as usize;

  // Cope's Rule parameters
  let c1 = 0.33;      // log-lambda intercept
  let c2 = 1.30;      // log-size intercept
  let delta = 0.04;   // systematic bias

  // Monte Carlo distribution parameters
  let sigma = 0.63;   // variance
  let alpha = 0.30;   // power-law tail

  // Extinctino parameters
  let beta = 1.0/(n as f64);    // baseline extinction rate
  let rho = 0.025;              // rate of extinction increase

  let doom_timer = |mass: f64| {
    // Determines how many more rounds until a species goes extinct
    // by drawing from a geometric distribution
    let p = (10.0 as f64).powf(rho * mass.log10() + beta.log10());

    let x: f64 = random();
    (x.log10() / (1.0 - p).log10()).ceil() as usize
  };

  // Vec<(id, mass, mass floor, extinction "date")>
  let mut extant: Vec<(usize, f64, f64, usize)> = Vec::with_capacity((n as f64 * 1.5).ceil() as usize);
  let doom = doom_timer(x_0);
  extant.push((1, x_0, x_min, doom));

  let mut all_species: Vec<Species> = Vec::with_capacity(0);
  if write_all {
    all_species = Vec::with_capacity(t_max + 1);
    all_species.push(Species{id: 1, mass: x_0, min_mass: x_min, death: doom, parent: 0});
  }


  // Start timing
  let start = time::precise_time_ns();

  for step in 1..t_max {
    let ancestor: usize = (random::<f64>() * extant.len() as f64).floor() as usize;

    let (id_a, mass_a, min, _) = extant[ancestor];
    x_min = min;

    // See if we've evolved a floor-raising characteristic
    if ratchet {
      if random::<f64>() < r_prob {
        x_min = x_min * (1.0 + r_magnitude);
      }
    }


    // Apply Cope's rule, strengthened for smaller species
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

    let mass_d = mass_a * tt;
    let doom = doom_timer(mass_d) + step;
    extant.push((step, mass_d, x_min, doom));
    if write_all {
      all_species.push(Species{
        id: step + 1,
        mass: mass_d,
        min_mass: x_min,
        death: doom,
        parent: id_a
      });
    }

    // Retain species whose doom timer is later than when we are now
    // ... that is to say, kill off the now-dead ones
    extant.retain(|&(_, _, _, d)| d > step);
  }

  // End timing
  let end = time::precise_time_ns();

  // Print out our final set of extant species
  let path = format!("extant_{}_{}_{}.csv", x_min, x_0, n);
  let mut writer = Writer::from_file(path).unwrap();
  for s in extant.into_iter() {
    writer.encode(s).ok().expect("CSV writer error");
  }

  if write_all {
    // All species from the HashMap
    let path = format!("all_{}_{}_{}.csv", x_min, x_0, n);
    let mut writer = Writer::from_file(path).unwrap();
    for s in all_species.into_iter() {
      writer.encode(s).ok().expect("CSV writer error");
    }
  }

  println!("Ran in {} seconds.", (end - start) / 1000000000);
}
