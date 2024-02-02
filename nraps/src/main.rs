use std::{fs, io::BufRead};

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
    let file = fs::read("../SampleInputFile.txt").expect("Unable to read the file");
    let lines: Vec<String> = file
        .lines()
        .map(|x| x.expect("Unable to read line").trim().to_ascii_lowercase())
        .filter(|x| !x.starts_with("#") && x.contains("="))
        .collect();

    let (var_names, var_values): (Vec<&str>, Vec<&str>) =
        lines.iter().map(|x| x.split_once("=").unwrap()).unzip();

    let var_names: Vec<String> = var_names.into_iter().map(|x| x.to_string()).collect();

    let positions = positions(var_names);

    let variables = Variables {
        solution: var_values[positions[0]].trim().parse().unwrap(),
        testcase: var_values[positions[1]].trim().parse().unwrap(),
        analk: var_values[positions[2]].trim().parse().unwrap(),
        energygroups: var_values[positions[3]].trim().parse().unwrap(),
        solver: match var_values[positions[4]].trim() {
            "1" => Solver::Gaussian,
            "2" => Solver::Jacobian,
            "3" => Solver::Sor,
            _ => Solver::LinAlg,
        },
        generations: var_values[positions[5]].trim().parse().unwrap(),
        histories: var_values[positions[6]].trim().parse().unwrap(),
        skip: var_values[positions[7]].trim().parse().unwrap(),
        numass: var_values[positions[8]].trim().parse().unwrap(),
        numrods: var_values[positions[9]].trim().parse().unwrap(),
        roddia: var_values[positions[10]].trim().parse().unwrap(),
        rodpitch: var_values[positions[11]].trim().parse().unwrap(),
        mpfr: var_values[positions[12]].trim().parse().unwrap(),
        mpwr: var_values[positions[13]].trim().parse().unwrap(),
        boundl: var_values[positions[14]].trim().parse().unwrap(),
        boundr: var_values[positions[15]].trim().parse().unwrap(),
    };

    // let indices = var_names
    //     .into_iter()
    //     .position(|x| x.contains(variable_names.iter().map(|x| x)).unwrap())
    //     .collect();

    for index in 0..lines.len() {
        println!("{}", lines[index]);
    }
}

fn positions(vector: Vec<String>) -> Vec<usize> {
    let variable_names: [&str; 23] = [
        "solution",
        "testcase",
        "analk",
        "energygroups",
        "solver",
        "generations",
        "histories",
        "skip",
        "numass",
        "numrods",
        "roddia",
        "rodpitch",
        "mpfr",
        "mpwr",
        "boundl",
        "boundr",
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
            vector
                .iter()
                .position(|x| x.trim() == variable_names[index])
                .unwrap(),
        );
        println!("{}", positions[index]);
    }

    positions
}

fn main() {
    let file_path = "../SampleInputFile.txt";

    process_input();

    println!("In file {}", file_path);

    //let contents = fs::read_to_string(file_path).expect("Should have been able to read the file");

    //println!("With text:\n{contents}");
}
