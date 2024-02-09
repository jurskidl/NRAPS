use std::fs::File;
use std::io::{BufReader, BufRead};
use std::collections::HashMap;
// Use these for timing
use std::thread::sleep;
use core::time::Duration;
use std::time::{SystemTime, UNIX_EPOCH};

pub enum Solver {
    LinAlg,
    Gaussian,
    Jacobian,
    Sor,
}

struct Variables {
    solution: u8,
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

fn process_input() -> (Variables, XSData, Vec<u8>) {
    let file = File::open("../SampleInputFile.txt").expect("Unable to read the file");
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader
        .lines()
        .map(|x| x.expect("Unable to read line").trim().to_ascii_lowercase())
        .filter(|x| !x.starts_with("#") && x.contains("="))
        .collect();

    let (vector, hash): (Vec<String>, Vec<String>) = lines
        .into_iter()
        .partition(|x| x.contains("matid"));

    let matid = get_mats(vector);

    let vars = hash.into_iter().map(|a| {
        let (key, value) = a.split_once("=").unwrap();
        (key.trim().to_string(), value.to_string())
    }).collect::<HashMap<String,String>>();

    let variables = Variables {
        solution: vars.get("solution").unwrap().trim().parse().unwrap(),
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

    let xsdata = XSData {
        sigtr: vars.get("sigtr").unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(),
        sigis: vars.get("sigis").unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(),
        sigds: vars.get("sigds").unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(),
        siga: vars.get("siga").unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(),
        sigf: vars.get("sigf").unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(),
        nut: vars.get("nut").unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(),
        chit: vars.get("chit").unwrap().split_whitespace().map(|x| x.parse().unwrap()).collect(),
    };
    
    (variables,xsdata,matid)
}

fn get_mats(vector: Vec<String>) -> Vec<u8> {
    let (var_names, var_values): (Vec<&str>, Vec<&str>) =
        vector.iter().map(|x| x.split_once("=").unwrap()).unzip();

    let mat_pos = var_names
        .iter()
        .enumerate()
        .filter_map(|(index, &ref x)| (x.trim() == "matid").then(|| index))
        .collect::<Vec<usize>>();

    let mut temp: Vec<Vec<u8>> = vec![];
    for index in 0..mat_pos.len() {
        temp.push(
            var_values[mat_pos[index]]
                .split_whitespace()
                .map(|y| y.to_owned().parse::<u8>().unwrap())
                .collect::<Vec<u8>>(),
        )
    }

    // index vector via mat# = matid[numass*((2*numrods)+1)]
    let matid: Vec<u8> = temp.into_iter().flatten().collect::<Vec<u8>>();

    matid
}

fn main() {
    // let (variables, xsdata, matid) = process_input();


    // below is for timing
    let now = SystemTime::now();
    let since_the_epoch = now
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards");

    for zyn in 0..1000000 {
        let (variables, xsdata, matid) = process_input();
        print!("{}\n", zyn);
    }

    print!("{}\n", since_the_epoch.as_millis());

    sleep(Duration::new(30, 0));
}