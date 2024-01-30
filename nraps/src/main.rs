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
    sigtr: f64,
    sigis: f64,
    sigds: f64,
    siga: f64,
    sigf: f64,
    nut: f64,
    chit: f64,
}

fn process_input(file_path: &str) -> Result<Variables, Box<dyn std::error::Error>> {
    let file = File::open(file_path).expect("Unable to open the specified file");
    let reader = BufReader::new(file);

    // Initialize the values necessary for use later
    let mut variables = Variables {
        solution: 0,
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
        sigtr: 1.0,
        sigis: 1.0,
        sigds: 1.0,
        siga: 1.0,
        sigf: 1.0,
        nut: 1.0,
        chit: 1.0,
    };
    let mut mox = XSData{
        sigtr: 1.0,
        sigis: 1.0,
        sigds: 1.0,
        siga: 1.0,
        sigf: 1.0,
        nut: 1.0,
        chit: 1.0,
    };
    let mut h2o = XSData{
        sigtr: 1.0,
        sigis: 1.0,
        sigds: 1.0,
        siga: 1.0,
        sigf: 1.0,
        nut: 1.0,
        chit: 1.0,
    };
    let mut cr = XSData{
        sigtr: 1.0,
        sigis: 1.0,
        sigds: 1.0,
        siga: 1.0,
        sigf: 1.0,
        nut: 1.0,
        chit: 1.0,
    };
    let mut case: u8 = 0;
    let mut xsdata_flag: bool = false;

    for line in reader.lines() {
        let line: String = line.unwrap();

        if !line.starts_with("#") && line.contains("=") && xsdata_flag == false {
            let split_line: Vec<String> = line
                .to_ascii_lowercase()
                .split("=")
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();
            match split_line[0].trim() {
                "solution" => variables.solution = split_line[1].trim().parse::<u8>().unwrap(),
                "testcase" => case = split_line[1].trim().parse::<u8>().unwrap(),
                "analk" => variables.analk = split_line[1].trim().parse::<u8>().unwrap(),
                "energygroups" => {
                    variables.energygroups = split_line[1].trim().parse::<u8>().unwrap()
                }
                "solver" => {
                    variables.solver = match split_line[1].trim() {
                        "1" => Solver::Gaussian,
                        "2" => Solver::Jacobian,
                        "3" => Solver::Sor,
                        _ => Solver::LinAlg,
                    }
                }
                "generations" => {
                    variables.generations = split_line[1].trim().parse::<u32>().unwrap()
                }
                "histories" => variables.histories = split_line[1].trim().parse::<u32>().unwrap(),
                "skip" => variables.skip = split_line[1].trim().parse::<u16>().unwrap(),
                "numass" => variables.numass = split_line[1].trim().parse::<u8>().unwrap(),
                "numrods" => variables.numrods = split_line[1].trim().parse::<u8>().unwrap(),
                "roddia" => variables.roddia = split_line[1].trim().parse::<f64>().unwrap(),
                "rodpitch" => variables.rodpitch = split_line[1].trim().parse::<f64>().unwrap(),
                "mpfr" => variables.mpfr = split_line[1].trim().parse::<u8>().unwrap(),
                "mpwr" => variables.mpwr = split_line[1].trim().parse::<u8>().unwrap(),
                "boundl" => variables.boundl = split_line[1].trim().parse::<f64>().unwrap(),
                "boundr" => variables.boundr = split_line[1].trim().parse::<f64>().unwrap(),
                "case" => {
                    if split_line[1].trim().parse::<u8>().unwrap() == case {
                            xsdata_flag = true;
                        } else {
                            xsdata_flag = false;
                        }
                    }
                _ => continue,
            }
        } else if !line.starts_with("#") && line.contains("=") && xsdata_flag == true {
            let split_line = line
                .to_ascii_lowercase()
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();
            match split_line[0].trim() {
                "sigtr" => {
                    uo2.sigtr = split_line[2].parse::<f64>().unwrap();
                    mox.sigtr = split_line[3].parse::<f64>().unwrap();
                    h2o.sigtr = split_line[4].parse::<f64>().unwrap();
                    cr.sigtr = split_line[5].parse::<f64>().unwrap();
                }
                "sigis" => {
                    uo2.sigis = split_line[2].parse::<f64>().unwrap();
                    mox.sigis = split_line[3].parse::<f64>().unwrap();
                    h2o.sigis = split_line[4].parse::<f64>().unwrap();
                    cr.sigis = split_line[5].parse::<f64>().unwrap();
                }
                "sigds" => {
                    uo2.sigds = split_line[2].parse::<f64>().unwrap();
                    mox.sigds = split_line[3].parse::<f64>().unwrap();
                    h2o.sigds = split_line[4].parse::<f64>().unwrap();
                    cr.sigds = split_line[5].parse::<f64>().unwrap();
                }
                "siga" => {
                    uo2.siga = split_line[2].parse::<f64>().unwrap();
                    mox.siga = split_line[3].parse::<f64>().unwrap();
                    h2o.siga = split_line[4].parse::<f64>().unwrap();
                    cr.siga = split_line[5].parse::<f64>().unwrap();
                }
                "sigf" => {
                    uo2.sigf = split_line[2].parse::<f64>().unwrap();
                    mox.sigf = split_line[3].parse::<f64>().unwrap();
                    h2o.sigf = split_line[4].parse::<f64>().unwrap();
                    cr.sigf = split_line[5].parse::<f64>().unwrap();
                }
                "nut" => {
                    uo2.nut = split_line[2].parse::<f64>().unwrap();
                    mox.nut = split_line[3].parse::<f64>().unwrap();
                    h2o.nut = split_line[4].parse::<f64>().unwrap();
                    cr.nut = split_line[5].parse::<f64>().unwrap();
                }
                "chit" => {
                    uo2.chit = split_line[2].parse::<f64>().unwrap();
                    mox.chit = split_line[3].parse::<f64>().unwrap();
                    h2o.chit = split_line[4].parse::<f64>().unwrap();
                    cr.chit = split_line[5].parse::<f64>().unwrap();
                }
                _ => xsdata_flag = false,
            }
        }
        //print!("{}\n", case);
    }
    print!("{}\n",uo2.sigtr);
    return Ok(variables);
}

fn main() {
    let file_path = "../SampleInputFile.txt";

    let _variables = process_input(&file_path).unwrap();

    println!("In file {}", file_path);

    //let contents = fs::read_to_string(file_path).expect("Should have been able to read the file");

    //println!("With text:\n{contents}");
}
