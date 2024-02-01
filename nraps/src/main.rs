use std::{
    fs::File,
    io::{BufRead, BufReader},
};

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

fn process_input(file_path: &str) -> Result<(Variables,XSData, XSData, XSData, XSData, Vec<u8>), Box<dyn std::error::Error>> {
    let file = File::open(file_path).expect("Unable to open the specified file");
    let reader = BufReader::new(file);

    // Initialize the values necessary for use later
    let mut variables = Variables {
        solution: 0,
        testcase: 0,
        analk: 1,
        energygroups: 1,
        solver: Solver::LinAlg,
        generations: 1,
        histories: 1,
        skip: 1,
        numass: 1,
        numrods: 1,
        roddia: 1.0,
        rodpitch: 1.0,
        mpfr: 1,
        mpwr: 1,
        boundl: 1.0,
        boundr: 1.0
    };
    let mut uo2 = XSData{
        sigtr: Vec::new(),
        sigis: Vec::new(),
        sigds: Vec::new(),
        siga: Vec::new(),
        sigf: Vec::new(),
        nut: Vec::new(),
        chit: Vec::new(),
    };
    let mut mox = XSData{
        sigtr: Vec::new(),
        sigis: Vec::new(),
        sigds: Vec::new(),
        siga: Vec::new(),
        sigf: Vec::new(),
        nut: Vec::new(),
        chit: Vec::new(),
    };
    let mut h2o = XSData{
        sigtr: Vec::new(),
        sigis: Vec::new(),
        sigds: Vec::new(),
        siga: Vec::new(),
        sigf: Vec::new(),
        nut: Vec::new(),
        chit: Vec::new(),
    };
    let mut cr = XSData{
        sigtr: Vec::new(),
        sigis: Vec::new(),
        sigds: Vec::new(),
        siga: Vec::new(),
        sigf: Vec::new(),
        nut: Vec::new(),
        chit: Vec::new(),
    };

    let mut materials: Vec<u8> = vec![];

    for line in reader.lines() {
        let line: String = line
            .unwrap()
            .to_ascii_lowercase();

        if !line.starts_with("#") && line.contains("=") {
            let (vars_name, vars_value) = line
                .split_once("=")
                .unwrap();
            match vars_name.trim() {
                "solution" => variables.solution = vars_value.trim().parse::<u8>().unwrap(),
                "testcase" => variables.testcase = vars_value.trim().parse::<u8>().unwrap(),
                "analk" => variables.analk = vars_value.trim().parse::<u8>().unwrap(),
                "energygroups" => variables.energygroups = vars_value.trim().parse::<u8>().unwrap(),
                "solver" => {
                    variables.solver = match vars_value.trim() {
                        "1" => Solver::Gaussian,
                        "2" => Solver::Jacobian,
                        "3" => Solver::Sor,
                        _ => Solver::LinAlg,
                    }
                }
                "generations" => {
                    variables.generations = vars_value.trim().parse::<u32>().unwrap()
                }
                "histories" => variables.histories = vars_value.trim().parse::<u32>().unwrap(),
                "skip" => variables.skip = vars_value.trim().parse::<u16>().unwrap(),
                "numass" => variables.numass = vars_value.trim().parse::<u8>().unwrap(),
                "numrods" => variables.numrods = vars_value.trim().parse::<u8>().unwrap(),
                "roddia" => variables.roddia = vars_value.trim().parse::<f64>().unwrap(),
                "rodpitch" => variables.rodpitch = vars_value.trim().parse::<f64>().unwrap(),
                "mpfr" => variables.mpfr = vars_value.trim().parse::<u8>().unwrap(),
                "mpwr" => variables.mpwr = vars_value.trim().parse::<u8>().unwrap(),
                "boundl" => variables.boundl = vars_value.trim().parse::<f64>().unwrap(),
                "boundr" => variables.boundr = vars_value.trim().parse::<f64>().unwrap(),
                "sigtr" => {
                    let split_vars = vars_value
                    .split_whitespace()
                    .map(|s| s.to_owned())
                    .collect::<Vec<String>>();
                    for index in 0..=variables.energygroups-1 {
                        uo2.sigtr.push(split_vars[(index*4) as usize].parse::<f64>().unwrap());
                        mox.sigtr.push(split_vars[(1+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.sigtr.push(split_vars[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.sigtr.push(split_vars[(3+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                }
                "sigis" => {
                    let split_vars = vars_value
                    .split_whitespace()
                    .map(|s| s.to_owned())
                    .collect::<Vec<String>>();
                    for index in 0..=variables.energygroups-1 {
                        uo2.sigis.push(split_vars[(index*4) as usize].parse::<f64>().unwrap());
                        mox.sigis.push(split_vars[(1+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.sigis.push(split_vars[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.sigis.push(split_vars[(3+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                }
                "sigds" => {
                    let split_vars = vars_value
                    .split_whitespace()
                    .map(|s| s.to_owned())
                    .collect::<Vec<String>>();
                    for index in 0..=variables.energygroups-1 {
                        uo2.sigds.push(split_vars[(index*4) as usize].parse::<f64>().unwrap());
                        mox.sigds.push(split_vars[(1+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.sigds.push(split_vars[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.sigds.push(split_vars[(3+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                }
                "siga" => {
                    let split_vars = vars_value
                    .split_whitespace()
                    .map(|s| s.to_owned())
                    .collect::<Vec<String>>();
                    for index in 0..=variables.energygroups-1 {
                        uo2.siga.push(split_vars[(index*4) as usize].parse::<f64>().unwrap());
                        mox.siga.push(split_vars[(1+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.siga.push(split_vars[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.siga.push(split_vars[(3+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                }
                "sigf" => {
                    let split_vars = vars_value
                    .split_whitespace()
                    .map(|s| s.to_owned())
                    .collect::<Vec<String>>();
                    for index in 0..=variables.energygroups-1 {
                        uo2.sigf.push(split_vars[(index*4) as usize].parse::<f64>().unwrap());
                        mox.sigf.push(split_vars[(1+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.sigf.push(split_vars[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.sigf.push(split_vars[(3+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                }
                "nut" => {
                    let split_vars = vars_value
                    .split_whitespace()
                    .map(|s| s.to_owned())
                    .collect::<Vec<String>>();
                    for index in 0..=variables.energygroups-1 {
                        uo2.nut.push(split_vars[(index*4) as usize].parse::<f64>().unwrap());
                        mox.nut.push(split_vars[(1+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.nut.push(split_vars[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.nut.push(split_vars[(3+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                }
                "chit" => {
                    let split_vars = vars_value
                    .split_whitespace()
                    .map(|s| s.to_owned())
                    .collect::<Vec<String>>();
                    for index in 0..=variables.energygroups-1 {
                        uo2.chit.push(split_vars[(index*4) as usize].parse::<f64>().unwrap());
                        mox.chit.push(split_vars[(1+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.chit.push(split_vars[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.chit.push(split_vars[(3+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                }
                "matid" => {
                    let split_vars = vars_value
                    .split_whitespace()
                    .map(|s| s.to_owned())
                    .collect::<Vec<String>>();
                for index in 0..split_vars.len() {
                materials.push(split_vars[index].parse::<u8>().unwrap());}
                }
                _ => continue,
            }
        }
    }
    return Ok((variables, uo2, mox, h2o, cr,materials));
}

fn main() {
    let file_path = "../SampleInputFile.txt";

    let (_variables, _uo2, _mox, _h2o, _cr, _materials) = process_input(&file_path).unwrap();

    println!("In file {}", file_path);

    //let contents = fs::read_to_string(file_path).expect("Should have been able to read the file");

    //println!("With text:\n{contents}");
}
