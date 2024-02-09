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
        .partition(|x| x.contains("matid") || x.contains("sigtr") || x.contains("sigis") || x.contains("sigds") || x.contains("siga") || x.contains("sigf") || x.contains("nut") || x.contains("chit"));

    let (matid, xsdata) = get_data(vector);

    let vars = hash.into_iter().map(|a| {
        let (key, value) = a.split_once("=").unwrap();
        (key.trim().to_string(), value.to_string())
    }).collect::<HashMap<String,String>>();

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
    
    (variables,xsdata,matid)
}

fn get_data(vector: Vec<String>) -> (Vec<u8>, XSData) {
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

    // index vector via mat# = matid[config#*numass*((2*numrods)+1)]
    let matid: Vec<u8> = temp.into_iter().flatten().collect::<Vec<u8>>();

    let variable_names: [&str; 7] = [
        "sigtr",
        "sigis",
        "sigds",
        "siga",
        "sigf",
        "nut",
        "chit",
    ];

    let mut positions: Vec<usize> = vec![];

    for index in 0..variable_names.len() {
        positions.push(
            var_names
                .iter()
                .position(|x| x.trim() == variable_names[index])
                .unwrap(),
        );
    }

    // index into vectors via desired_xs = sigtr[(mat# + mattypes*energygroup) as usize]
    let xsdata = XSData {
        sigtr: var_values[positions[0]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        sigis: var_values[positions[1]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        sigds: var_values[positions[2]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        siga: var_values[positions[3]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        sigf: var_values[positions[4]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        nut: var_values[positions[5]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        chit: var_values[positions[6]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
    };


    (matid, xsdata)
}

fn main() {
    let now = SystemTime::now();
    for zyn in 0..1000000 {
        let (variables, xsdata, matid) = process_input();
        print!("{}\n", zyn);
    }
    
    print!("{}", now.elapsed().unwrap().as_secs());
}