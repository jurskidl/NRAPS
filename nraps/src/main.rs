use std::time::SystemTime;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::collections::HashMap;

pub enum Solver {
    LinAlg,
    Gaussian,
    Jacobian,
    Sor,
}

struct Variables {
    solution: u8,
    testcase: u8,
    analk: u8,
    mattypes: u8,
    energygroups: u8,
    solver: Solver,
    generations: u32,
    histories: u32,
    skip: u16,
    numass: u8,
    numrods: u8,
    roddia: f64,
    rodpitch: f64,
    mpfr: u8,
    mpwr: u8,
    boundl: f64,
    boundr: f64,
}

struct XSData {
    sigtr: Vec<f64>,
    sigis: Vec<f64>,
    sigds: Vec<f64>,
    siga: Vec<f64>,
    sigf: Vec<f64>,
    nut: Vec<f64>,
    chit: Vec<f64>,
}

fn process_input() {
    let file = File::open("../SampleInputFile.txt").expect("Unable to read the file");
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader.lines().map(|x| x.expect("Unable to read line").trim().to_ascii_lowercase()).filter(|x| !x.starts_with("#") && x.contains("=") && x.contains("solution")).collect();

    let vars =
        lines.iter().map(|x| x.split_whitespace().collect()).collect::<HashMap<&str, &str>>();

        println!("{}", vars.get("solution").unwrap().to_owned().parse::<u8>().unwrap());

    // let variables = Variables {
    //     solution: vars.get("solution ").unwrap().to_owned().parse().unwrap(),
    //     testcase: vars.get("testcase ").unwrap().to_owned().parse().unwrap(),
    //     analk: vars.get("analk ").unwrap().to_owned().parse().unwrap(),
    //     mattypes: vars.get("mattypes ").unwrap().to_owned().parse().unwrap(),
    //     energygroups: vars.get("energygroups ").unwrap().to_owned().parse().unwrap(),
    //     solver: match vars.get("solver ").unwrap().to_owned().trim() {
    //         "1" => Solver::Gaussian,
    //         "2" => Solver::Jacobian,
    //         "3" => Solver::Sor,
    //         _ => Solver::LinAlg,
    //     },
    //     generations: vars.get("generations ").unwrap().to_owned().parse().unwrap(),
    //     histories: vars.get("histories ").unwrap().to_owned().parse().unwrap(),
    //     skip: vars.get("skip ").unwrap().to_owned().parse().unwrap(),
    //     numass: vars.get("numass ").unwrap().to_owned().parse().unwrap(),
    //     numrods: vars.get("numrods ").unwrap().to_owned().parse().unwrap(),
    //     roddia: vars.get("roddia ").unwrap().to_owned().parse().unwrap(),
    //     rodpitch: vars.get("rodpitch ").unwrap().to_owned().parse().unwrap(),
    //     mpfr: vars.get("mpfr ").unwrap().to_owned().parse().unwrap(),
    //     mpwr: vars.get("mpwr ").unwrap().to_owned().parse().unwrap(),
    //     boundl: vars.get("boundl ").unwrap().to_owned().parse().unwrap(),
    //     boundr: vars.get("boundr ").unwrap().to_owned().parse().unwrap(),
    // };

    //variables
}

fn main() {
    let now = SystemTime::now();
    for zyn in 0..1000000 {
       process_input();
        print!("{}\n", zyn);
    }

    print!("{}", now.elapsed().unwrap().as_secs());
    //let file_path = "../SampleInputFile.txt";

    //println!("In file {}", file_path);
}