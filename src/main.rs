
use std::collections::HashMap;

extern crate rand;
extern crate time;
extern crate csv;

use rand::random;
use csv::Writer;

#[derive(Clone, Debug)]
pub struct Species {
  id: i64,            // evolution order, birthdate
  mass: f64,          // mass (g)
  death: i64,         // iteration went (or will go) extinct
}

const N = 5000;  // num species at equil, w/o size bias

fn speciate() {

}

fn doom_timer(mass: f64) -> i64 {
  // Determines how many more rounds until a species goes extinct
  // by drawing from a geometric distribution

  let beta = 1.0/(N as f64);     // baseline extinction rate
  let rho = 0.025;    // rate of extinction increase
  let p = (10.0 as f64).powf(rho * mass.log10() + beta.log10());

  let x: f64 = random();
  (x.log10() / (1.0 - p).log10()).ceil() as i64
}

fn main() {
  // The likely-to-be-toyed-with variables
  let x_min = 1.8;
  let x_0 = 40.0;

  // Determine how many steps to run
  let nu = 1.6;       // mean species lifetime (My)
  let tau = 60.0;     // total simulation time (My)
  let t_max = ((tau / nu) * (N as f64)).round();

  // Cope's Rule parameters
  let c1 = 0.33;      // log-lambda intercept
  let c2 = 1.30;      // log-size intercept
  let delta = 0.04;   // systematic bias

  // Monte Carlo distribution parameters
  let sigma = 0.63;   // variance
  let alpha = 0.30;   // power-law tail

  let mut extant: Vec<f64> = 
    Vec::with_capacity((N as f64 * 1.5).ceil() as usize);
  let mut all_species: Vec<Species> =
    Vec::with_capacity(N * 10);

  // Start timing
  let start = time::precise_time_ns();



  let mut deaths = Vec::with_capacity(10000);
  for ii in 0..10000 {
    deaths.push(doom_timer(x_0));
  }

  let path = format!("deaths{}.csv", x_0);
  let mut writer = Writer::from_file(path).unwrap();
  for d in deaths.into_iter() {
    writer.encode(d).ok().expect("CSV writer error");
  }



  // End timing
  let end = time::precise_time_ns();

  println!("Ran in {} seconds.", (end - start) / (1000000000));
}
