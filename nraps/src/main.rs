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
    };
    let mut case = "0";

    for line in reader.lines() {
        let line: String = line.unwrap().to_ascii_lowercase();

        if !line.starts_with("#") && line.contains("=") {
            let split_line: Vec<String> = line
                .split("=")
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();
            match split_line[0].trim() {
                "solution" => variables.solution = split_line[1].trim().parse::<u8>().unwrap(),
                "testcase" => case = split_line[1].trim(),
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
                _ => continue,
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
