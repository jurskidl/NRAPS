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

    let (var_names, _var_values): (Vec<&str>, Vec<&str>) =
        lines.iter().map(|x| x.split_once("=").unwrap()).unzip();

    let var_names: Vec<String> = var_names.into_iter().map(|x| x.to_string()).collect();

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

    // let indices = var_names
    //     .into_iter()
    //     .position(|x| x.contains(variable_names.iter().map(|x| x)).unwrap())
    //     .collect();

    // for index in 0.._var_names.len() {
    //     println!("{}", _var_names[index]);
    // }
}

fn main() {
    let file_path = "../SampleInputFile.txt";

    process_input();

    println!("In file {}", file_path);

    //let contents = fs::read_to_string(file_path).expect("Should have been able to read the file");

    //println!("With text:\n{contents}");
}
