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
    sigtr: Vec<f64>,
    sigis: Vec<f64>,
    sigds: Vec<f64>,
    siga: Vec<f64>,
    sigf: Vec<f64>,
    nut: Vec<f64>,
    chit: Vec<f64>,
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
    let mut case: u8 = 0;
    let mut xsdata_flag: bool = false;

    for line in reader.lines() {
        let line: String = line
            .unwrap()
            .to_ascii_lowercase();

        if !line.starts_with("#") && line.contains("=") && xsdata_flag == false {
            let (vars_name, vars_value) = line
                .split_once("=")
                .unwrap();
            match vars_name.trim() {
                "solution" => variables.solution = vars_value.trim().parse::<u8>().unwrap(),
                "testcase" => case = vars_value.trim().parse::<u8>().unwrap(),
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
                "case" => {
                    if vars_value.trim().parse::<u8>().unwrap() == case {
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
            for index in 0..=variables.energygroups-1 {
                match split_line[0].trim() {
                    "sigtr" => {
                        print!("{}\n",split_line[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        uo2.sigtr.push(split_line[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        mox.sigtr.push(split_line[(3+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.sigtr.push(split_line[(4+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.sigtr.push(split_line[(5+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                    "sigis" => {
                        uo2.sigis.push(split_line[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        mox.sigis.push(split_line[(3+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.sigis.push(split_line[(4+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.sigis.push(split_line[(5+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                    "sigds" => {
                        uo2.sigds.push(split_line[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        mox.sigds.push(split_line[(3+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.sigds.push(split_line[(4+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.sigds.push(split_line[(5+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                    "siga" => {
                        uo2.siga.push(split_line[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        mox.siga.push(split_line[(3+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.siga.push(split_line[(4+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.siga.push(split_line[(5+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                    "sigf" => {
                        uo2.sigf.push(split_line[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        mox.sigf.push(split_line[(3+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.sigf.push(split_line[(4+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.sigf.push(split_line[(5+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                    "nut" => {
                        uo2.nut.push(split_line[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        mox.nut.push(split_line[(3+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.nut.push(split_line[(4+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.nut.push(split_line[(5+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                    "chit" => {
                        uo2.chit.push(split_line[(2+(index*4)) as usize].parse::<f64>().unwrap());
                        mox.chit.push(split_line[(3+(index*4)) as usize].parse::<f64>().unwrap());
                        h2o.chit.push(split_line[(4+(index*4)) as usize].parse::<f64>().unwrap());
                        cr.chit.push(split_line[(5+(index*4)) as usize].parse::<f64>().unwrap());
                    }
                    _ => xsdata_flag = false,
                }
            }
        }
    }
    return Ok(variables);
}

fn main() {
    let file_path = "../SampleInputFile.txt";

    let _variables = process_input(&file_path).unwrap();

    println!("In file {}", file_path);

    //let contents = fs::read_to_string(file_path).expect("Should have been able to read the file");

    //println!("With text:\n{contents}");
}
