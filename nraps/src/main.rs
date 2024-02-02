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

fn process_input() -> (Variables,XSData, XSData, XSData, XSData) {
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

    let mattypes: u8 = var_values[positions[23]].trim().parse::<u8>().unwrap();

    let (sigtr, sigis, sigds, siga, sigf, nut, chit) = xsdata_processor(var_values, positions);

    let uo2 = XSData {
        sigtr: (0..variables.energygroups).into_iter().map(|x: u8| sigtr[(mattypes*x) as usize]).collect(),
        sigis: (0..variables.energygroups).into_iter().map(|x: u8| sigis[(mattypes*x) as usize]).collect(),
        sigds: (0..variables.energygroups).into_iter().map(|x: u8| sigds[(mattypes*x) as usize]).collect(),
        siga: (0..variables.energygroups).into_iter().map(|x: u8| siga[(mattypes*x) as usize]).collect(),
        sigf: (0..variables.energygroups).into_iter().map(|x: u8| sigf[(mattypes*x) as usize]).collect(),
        nut: (0..variables.energygroups).into_iter().map(|x: u8| nut[(mattypes*x) as usize]).collect(),
        chit: (0..variables.energygroups).into_iter().map(|x: u8| chit[(mattypes*x) as usize]).collect(),
    };

    let mox = XSData {
        sigtr: (0..variables.energygroups).into_iter().map(|x: u8| sigtr[((mattypes*x)+1) as usize]).collect(),
        sigis: (0..variables.energygroups).into_iter().map(|x: u8| sigis[((mattypes*x)+1) as usize]).collect(),
        sigds: (0..variables.energygroups).into_iter().map(|x: u8| sigds[((mattypes*x)+1) as usize]).collect(),
        siga: (0..variables.energygroups).into_iter().map(|x: u8| siga[((mattypes*x)+1) as usize]).collect(),
        sigf: (0..variables.energygroups).into_iter().map(|x: u8| sigf[((mattypes*x)+1) as usize]).collect(),
        nut: (0..variables.energygroups).into_iter().map(|x: u8| nut[((mattypes*x)+1) as usize]).collect(),
        chit: (0..variables.energygroups).into_iter().map(|x: u8| chit[((mattypes*x)+1) as usize]).collect(),
    };

    let h2o = XSData {
        sigtr: (0..variables.energygroups).into_iter().map(|x: u8| sigtr[((mattypes*x)+2) as usize]).collect(),
        sigis: (0..variables.energygroups).into_iter().map(|x: u8| sigis[((mattypes*x)+2) as usize]).collect(),
        sigds: (0..variables.energygroups).into_iter().map(|x: u8| sigds[((mattypes*x)+2) as usize]).collect(),
        siga: (0..variables.energygroups).into_iter().map(|x: u8| siga[((mattypes*x)+2) as usize]).collect(),
        sigf: (0..variables.energygroups).into_iter().map(|x: u8| sigf[((mattypes*x)+2) as usize]).collect(),
        nut: (0..variables.energygroups).into_iter().map(|x: u8| nut[((mattypes*x)+2) as usize]).collect(),
        chit: (0..variables.energygroups).into_iter().map(|x: u8| chit[((mattypes*x)+2) as usize]).collect(),
    };

    let cr = XSData {
        sigtr: (0..variables.energygroups).into_iter().map(|x: u8| sigtr[((mattypes*x)+3) as usize]).collect(),
        sigis: (0..variables.energygroups).into_iter().map(|x: u8| sigis[((mattypes*x)+3) as usize]).collect(),
        sigds: (0..variables.energygroups).into_iter().map(|x: u8| sigds[((mattypes*x)+3) as usize]).collect(),
        siga: (0..variables.energygroups).into_iter().map(|x: u8| siga[((mattypes*x)+3) as usize]).collect(),
        sigf: (0..variables.energygroups).into_iter().map(|x: u8| sigf[((mattypes*x)+3) as usize]).collect(),
        nut: (0..variables.energygroups).into_iter().map(|x: u8| nut[((mattypes*x)+3) as usize]).collect(),
        chit: (0..variables.energygroups).into_iter().map(|x: u8| chit[((mattypes*x)+3) as usize]).collect(),
    };

    // let indices = var_names
    //     .into_iter()
    //     .position(|x| x.contains(variable_names.iter().map(|x| x)).unwrap())
    //     .collect();

    for index in 0..lines.len() {
        println!("{}", lines[index]);
    }

    (variables, uo2, mox, h2o, cr)
}

fn positions(vector: Vec<String>) -> Vec<usize> {
    let variable_names: [&str; 24] = [
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
        "mattypes",
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

fn xsdata_processor(var_values: Vec<&str>, positions: Vec<usize>) -> (Vec<f64>,Vec<f64>,Vec<f64>,Vec<f64>,Vec<f64>,Vec<f64>,Vec<f64>){
    let sigtr = var_values[positions[16]].split_whitespace().map(|x| x.to_owned().parse::<f64>().unwrap()).collect::<Vec<f64>>();
    let sigis = var_values[positions[17]].split_whitespace().map(|x| x.to_owned().parse::<f64>().unwrap()).collect::<Vec<f64>>();
    let sigds = var_values[positions[18]].split_whitespace().map(|x| x.to_owned().parse::<f64>().unwrap()).collect::<Vec<f64>>();
    let siga = var_values[positions[19]].split_whitespace().map(|x| x.to_owned().parse::<f64>().unwrap()).collect::<Vec<f64>>();
    let sigf = var_values[positions[20]].split_whitespace().map(|x| x.to_owned().parse::<f64>().unwrap()).collect::<Vec<f64>>();
    let nut = var_values[positions[21]].split_whitespace().map(|x| x.to_owned().parse::<f64>().unwrap()).collect::<Vec<f64>>();
    let chit = var_values[positions[22]].split_whitespace().map(|x| x.to_owned().parse::<f64>().unwrap()).collect::<Vec<f64>>();

    (sigtr, sigis, sigds, siga, sigf, nut, chit)    
}

fn main() {
    let file_path = "../SampleInputFile.txt";

    process_input();

    println!("In file {}", file_path);

    //let contents = fs::read_to_string(file_path).expect("Should have been able to read the file");

    //println!("With text:\n{contents}");
}
