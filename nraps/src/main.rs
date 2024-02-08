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

struct _XSData {
    sigtr: Vec<f64>,
    sigis: Vec<f64>,
    sigds: Vec<f64>,
    siga: Vec<f64>,
    sigf: Vec<f64>,
    nut: Vec<f64>,
    chit: Vec<f64>,
}

fn process_input() -> Variables {
    let file = File::open("../SampleInputFile.txt").expect("Unable to read the file");
    let reader = BufReader::new(file);
    let lines: Vec<_> = reader
        .lines()
        .map(|x| x.expect("Unable to read line").trim().to_ascii_lowercase())
        .filter(|x| !x.starts_with("#") && x.contains("="))
        .collect();

        let vars = lines.into_iter().map(|a| {
            let (key, value) = a.split_once("=").unwrap();
            (key.trim().to_string(), value.to_string())
        }).collect::<HashMap<String,String>>();

    // let vars =
    //     lines.into_iter().map(|x| x.split_whitespace().unwrap()).collect::<HashMap<&str, &str>>();

    // for val in vars.keys() {
    //     println!("{val}");
    // }

    // print!("{}\n",vars.get("solution").unwrap().trim().parse::<u8>().unwrap());
    // print!("{}\n",vars.get("testcase").unwrap());
    // print!("{}\n",vars.get("analk").unwrap());
    // print!("{}\n",vars.get("mattypes").unwrap());
    // print!("{}\n",vars.get("energygroups").unwrap());
    // print!("{}\n",vars.get("generations").unwrap());
    // print!("{}\n",vars.get("histories").unwrap());
    // print!("{}\n",vars.get("skip").unwrap());
    // print!("{}\n",vars.get("numass").unwrap());
    // print!("{}\n",vars.get("numrods").unwrap());
    // print!("{}\n",vars.get("roddia").unwrap());
    // print!("{}\n",vars.get("rodpitch").unwrap());
    // print!("{}\n",vars.get("mpfr").unwrap());
    // print!("{}\n",vars.get("mpwr").unwrap());
    // print!("{}\n",vars.get("boundl").unwrap());
    // print!("{}\n",vars.get("boundr").unwrap());
    

    let variables = Variables {
        solution: vars.get("solution").unwrap().trim().parse().unwrap(),
        testcase: vars.get("testcase").unwrap().trim().parse().unwrap(),
        analk: vars.get("analk").unwrap().trim().parse().unwrap(),
        mattypes: vars.get("mattypes").unwrap().trim().parse().unwrap(),
        energygroups: vars.get("energygroups").unwrap().trim().parse().unwrap(),
        solver: match vars.get("solver").unwrap().trim() {
            "1" => Solver::Gaussian,
            "2" => Solver::Jacobian,
            "3" => Solver::Sor,
            _ => Solver::LinAlg,
        },
        generations: vars.get("generations").unwrap().trim().parse().unwrap(),
        histories: vars.get("histories").unwrap().trim().parse().unwrap(),
        skip: vars.get("skip").unwrap().trim().parse().unwrap(),
        numass: vars.get("numass").unwrap().trim().parse().unwrap(),
        numrods: vars.get("numrods").unwrap().trim().parse().unwrap(),
        roddia: vars.get("roddia").unwrap().trim().parse().unwrap(),
        rodpitch: vars.get("rodpitch").unwrap().trim().parse().unwrap(),
        mpfr: vars.get("mpfr").unwrap().trim().parse().unwrap(),
        mpwr: vars.get("mpwr").unwrap().trim().parse().unwrap(),
        boundl: vars.get("boundl").unwrap().trim().parse().unwrap(),
        boundr: vars.get("boundr").unwrap().trim().parse().unwrap(),
    };

    // println!("{}\n", variables.solution);
    // println!("{}\n", variables.testcase);
    // println!("{}\n", variables.analk);
    // println!("{}\n", variables.mattypes);
    // println!("{}\n", variables.energygroups);
    // println!("{}\n", variables.generations);
    // println!("{}\n", variables.histories);
    // println!("{}\n", variables.skip);
    // println!("{}\n", variables.numass);
    // println!("{}\n", variables.numrods);
    // println!("{}\n", variables.roddia);
    // println!("{}\n", variables.rodpitch);
    // println!("{}\n", variables.mpfr);
    // println!("{}\n", variables.mpwr);
    // println!("{}\n", variables.boundl);
    // println!("{}\n", variables.boundr);
    

    variables
}

fn main() {
    let now = SystemTime::now();
    for zyn in 0..1000000 {
        let variables = process_input();
        print!("{}\n", zyn);
    }

    print!("{}", now.elapsed().unwrap().as_secs());
    //let file_path = "../SampleInputFile.txt";

    //println!("In file {}", file_path);
}