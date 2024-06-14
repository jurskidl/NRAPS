use std::fs::File;

use memmap2::MmapOptions;

use crate::{DeltaX, Solver, Variables, XSData};

const EQUALS: u8 = 61;
const NEWLINE: u8 = 10;
const POUND: u8 = 35;

fn skip_line(mut pos: usize, end: usize, buffer: &[u8]) -> usize {
    while buffer[pos] != NEWLINE && pos < end {
        pos += 1;
    }
    pos + 1
}

fn get_index(key: &str) -> usize {
    return match key {
        "analk" => 0,
        "mattypes" => 1,
        "energygroups" => 2,
        "generations" => 3,
        "histories" => 4,
        "skip" => 5,
        "numass" => 6,
        "numrods" => 7,
        "roddia" => 8,
        "rodpitch" => 9,
        "mpfr" => 10,
        "mpwr" => 11,
        "boundl" => 12,
        "boundr" => 13,
        "sigt" =>  14,
        "sigs" => 15,
        "mu" => 16,
        "siga" => 17,
        "sigf" => 18,
        "nut" => 19,
        "chit" => 20,
        "scat" => 21,
        "matid" => 22,
        "solution" => 23,
        "solver" => 24,
        _ => 25,
    }
}


fn scan_ascii_chunk(buffer: &[u8]) -> [String; 26] {
    let end = buffer.len();

    let mut temp: [String; 26] = Default::default();

    let mut pos: usize = 0;
    let mut line_start: usize = 0;
    let mut name_end: usize = 0;
    let mut val_start: usize = 0;

    while pos < end {
        match buffer[pos] {
            POUND => {
                pos = skip_line(pos, end, buffer);
                line_start = pos;
            }
            EQUALS => {
                name_end = pos - 1;
                val_start = pos + 1;
            }
            NEWLINE => {
                if name_end > line_start {
                    let key = String::from_utf8_lossy(&buffer[line_start..name_end])
                        .trim()
                        .to_string()
                        .to_ascii_lowercase();
                    let value = String::from_utf8_lossy(&buffer[val_start..pos])
                        .trim()
                        .to_string();
                    temp[get_index(key.as_str())] += &(" ".to_owned() + &value);
                } else {}
                line_start = pos + 1;
            }
            _ => {}
        }
        pos += 1;
    }
    temp
}

pub fn process_input() -> (Variables, XSData, Vec<u8>, DeltaX, u8, Solver) {
    let file = File::open("./TestCaseC.txt").expect("Unable to read the file");
    let mapped_file = unsafe { MmapOptions::new().map(&file).unwrap() };
    let start: usize = 0;
    let end: usize = mapped_file.len();
    let temp = scan_ascii_chunk(&mapped_file[start..end]);

    let variables = Variables {
        analk: temp[0].trim().parse().unwrap(),
        mattypes: temp[1].trim().parse().unwrap(),
        energygroups: temp[2].trim().parse().unwrap(),
        generations: temp[3].trim().parse().unwrap(),
        histories: temp[4].trim().parse().unwrap(),
        skip: temp[5].trim().parse().unwrap(),
        numass: temp[6].trim().parse().unwrap(),
        numrods: temp[7].trim().parse().unwrap(),
        roddia: temp[8].trim().parse().unwrap(),
        rodpitch: temp[9].trim().parse::<f64>().unwrap()
            - temp[8].trim().parse::<f64>().unwrap(),
        mpfr: temp[10].trim().parse().unwrap(),
        mpwr: temp[11].trim().parse().unwrap(),
        boundl: temp[12].trim().parse().unwrap(),
        boundr: temp[13].trim().parse().unwrap(),
    };

    let deltax = DeltaX {
        fuel: variables.roddia / variables.mpfr as f64,
        water: variables.rodpitch / variables.mpwr as f64,
    };

    // index into vectors via desired_xs = sigtr[(mat# + (energygroup*mattypes) as usize]
    let mut xsdata = XSData {
        sigt: temp[14]
            .split_whitespace()
            .map(|x| x.parse::<f64>().unwrap())
            .collect(),
        sigs: temp[15]
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        mu: temp[16]
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        siga: temp[17]
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        sigf: temp[18]
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        nut: temp[19]
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        chit: temp[20]
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        // Index via [(mattype * energygroups) + ((energygroups * starting_energy) + final_energy)]
        scat_matrix: temp[21]
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        inv_sigtr: Vec::new(),
    };

    for index in 0..xsdata.sigt.len() {
        xsdata.inv_sigtr.push((xsdata.sigt[index] - xsdata.mu[index] * xsdata.sigs[index]).powi(-1));
    }

    let matid: Vec<u8> = temp[22]
        .split_ascii_whitespace()
        .map(|x| x.parse::<u8>().unwrap())
        .collect();
    (
        variables,
        xsdata,
        matid,
        deltax,
        temp[23].trim().parse().unwrap(),
        match temp[24].trim() {
            "1" => Solver::Gaussian,
            "2" => Solver::Jacobian,
            "3" => Solver::SR,
            _ => Solver::LinAlg,
        },
    )
}