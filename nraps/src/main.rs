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
    let file = fs::read("../SampleInputFile.txt").expect("Unable to read the file");
    let lines: Vec<String> = file
        .lines()
        .map(|x| x.expect("Unable to read line").trim().to_ascii_lowercase())
        .filter(|x| !x.starts_with("#") && x.contains("="))
        .collect();

    let (var_names, var_values): (Vec<&str>, Vec<&str>) =
        lines.iter().map(|x| x.split_once("=").unwrap()).unzip();

    let var_names: Vec<String> = var_names.into_iter().map(|x| x.to_string()).collect();

    let (positions, mat_pos) = positions(var_names);

    let variables = Variables {
        solution: var_values[positions[0]].trim().parse().unwrap(),
        testcase: var_values[positions[1]].trim().parse().unwrap(),
        analk: var_values[positions[2]].trim().parse().unwrap(),
        mattypes: var_values[positions[3]].trim().parse::<u8>().unwrap(),
        energygroups: var_values[positions[4]].trim().parse().unwrap(),
        solver: match var_values[positions[5]].trim() {
            "1" => Solver::Gaussian,
            "2" => Solver::Jacobian,
            "3" => Solver::Sor,
            _ => Solver::LinAlg,
        },
        generations: var_values[positions[6]].trim().parse().unwrap(),
        histories: var_values[positions[7]].trim().parse().unwrap(),
        skip: var_values[positions[8]].trim().parse().unwrap(),
        numass: var_values[positions[9]].trim().parse().unwrap(),
        numrods: var_values[positions[10]].trim().parse().unwrap(),
        roddia: var_values[positions[11]].trim().parse().unwrap(),
        rodpitch: var_values[positions[12]].trim().parse().unwrap(),
        mpfr: var_values[positions[13]].trim().parse().unwrap(),
        mpwr: var_values[positions[14]].trim().parse().unwrap(),
        boundl: var_values[positions[15]].trim().parse().unwrap(),
        boundr: var_values[positions[16]].trim().parse().unwrap(),
    };

    // index into vectors via desired_xs = sigtr[(mat# + mattypes*energygroup) as usize]
    let xsdata = XSData {
        sigtr: var_values[positions[17]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        sigis: var_values[positions[18]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        sigds: var_values[positions[19]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        siga: var_values[positions[20]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        sigf: var_values[positions[21]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        nut: var_values[positions[22]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
        chit: var_values[positions[23]]
            .split_whitespace()
            .map(|x| x.to_owned().parse::<f64>().unwrap())
            .collect::<Vec<f64>>(),
    };

    // index vector via mat# = matid[config#*numass*((2*numrods)+1)]
    let matid = get_mats(var_values, mat_pos);

    (variables, xsdata, matid)
}

fn positions(vector: Vec<String>) -> (Vec<usize>, Vec<usize>) {
    let variable_names: [&str; 24] = [
        "solution",
        "testcase",
        "analk",
        "mattypes",
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
        //println!("{}", positions[index]);
    }

    let mat_pos = vector
        .iter()
        .enumerate()
        .filter_map(|(index, &ref x)| (x.trim() == "matid").then(|| index))
        .collect::<Vec<usize>>();

    (positions, mat_pos)
}

fn get_mats(vector: Vec<&str>, mat_pos: Vec<usize>) -> Vec<u8> {
    let mut temp: Vec<Vec<u8>> = vec![];
    for index in 0..mat_pos.len() {
        temp.push(
            vector[mat_pos[index]]
                .split_whitespace()
                .map(|y| y.to_owned().parse::<u8>().unwrap())
                .collect::<Vec<u8>>(),
        )
    }

    let matid: Vec<u8> = temp.into_iter().flatten().collect::<Vec<u8>>();
    matid
}

fn main() {
    let file_path = "../SampleInputFile.txt";

    process_input();

    println!("In file {}", file_path);

    //let contents = fs::read_to_string(file_path).expect("Should have been able to read the file");

    //println!("With text:\n{contents}");
}
