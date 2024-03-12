use memmap2::MmapOptions;
use std::collections::HashMap;
use std::fs::File;
use std::iter::repeat;
// For Multithreading
use std::thread;
// Use these for timing
use std::time::SystemTime;

pub const NUM_VARS: usize = 24;
pub const EQUALS: u8 = 61;
pub const SPACE: u8 = 32;
pub const NEWLINE: u8 = 10;
pub const POUND: u8 = 35;

pub enum Solver {
    LinAlg,
    Gaussian,
    Jacobian,
    Sor,
}

struct Variables {
    solution: u8,
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
    mpfr: usize,
    mpwr: usize,
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
// If multithreading
// fn next_end_line(mut end: usize, buffer: &[u8]) -> usize {
//     while buffer[end] != NEWLINE && end < buffer.len() {
//         end += 1;
//     }
//     end + 1
// }

fn skip_line(mut pos: usize, end: usize, buffer: &[u8]) -> usize {
    while buffer[pos] != NEWLINE && pos < end {
        pos += 1;
    }
    pos + 1
}

fn scan_ascii_chunk(start: usize, end: usize, buffer: &[u8]) -> HashMap<String, String> {
    let mut hash: HashMap<String, String> = HashMap::with_capacity(NUM_VARS);

    let mut pos = start;
    let mut line_start = start;
    let mut name_end = start;
    let mut val_start = start;

    // If multithreading
    // if end != buffer.len() && buffer[end] != NEWLINE {
    //     end = next_end_line(end, buffer);
    // }

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
                    hash.entry(key)
                        .and_modify(|existing| *existing = existing.to_owned() + " " + &value) // I don't know why this works but the gods have smiled upon me
                        .or_insert(value);
                } else {
                }
                line_start = pos + 1;
            }
            _ => {}
        }
        pos += 1;
    }
    hash
}

fn process_input() -> (Variables, XSData, Vec<u8>) {
    let file = File::open("../SampleInputFile.txt").expect("Unable to read the file");
    let mapped_file = unsafe { MmapOptions::new().map(&file).unwrap() };
    let start: usize = 0;
    let end: usize = mapped_file.len();
    let hash: HashMap<String, String> = scan_ascii_chunk(start, end, &&mapped_file);

    // For Multithreading
    // let size: usize = mapped_file.len();
    // let threads: usize = thread::available_parallelism().unwrap().get();
    // let chunk_length = size / threads;
    // let starting_points: Vec<usize> = (0..threads).map(|x| x * chunk_length).collect();
    // let mut ending_points: Vec<usize> = Vec::from_iter(starting_points[1..threads].iter().cloned());
    // ending_points.push(size);

    // let mut hash: HashMap<String, String> = HashMap::with_capacity(NUM_VARS);
    // std::thread::scope(|scope| {
    //     let mut handles = Vec::with_capacity(threads);
    //     for thread in 0..threads {
    //         let start = starting_points[thread];
    //         let end = ending_points[thread];
    //         let buffer = &mapped_file;
    //         let handle = scope.spawn(move || scan_ascii_chunk(start, end, &buffer));
    //         handles.push(handle);
    //     }

    //     // Aggregate the results
    //     for handle in handles {
    //         let chunk_result = handle.join().unwrap();
    //         for (key, value) in chunk_result {
    //             hash.entry(key.trim().to_string())
    //                 .and_modify(|existing| *existing = existing.to_owned() + " " + &value)
    //                 .or_insert(value);
    //         }
    //     }
    // });

    let variables = Variables {
        solution: hash.get("solution").unwrap().trim().parse().unwrap(),
        analk: hash.get("analk").unwrap().trim().parse().unwrap(),
        mattypes: hash.get("mattypes").unwrap().trim().parse().unwrap(),
        energygroups: hash.get("energygroups").unwrap().trim().parse().unwrap(),
        solver: match hash.get("solver").unwrap().trim() {
            "1" => Solver::Gaussian,
            "2" => Solver::Jacobian,
            "3" => Solver::Sor,
            _ => Solver::LinAlg,
        },
        generations: hash.get("generations").unwrap().trim().parse().unwrap(),
        histories: hash.get("histories").unwrap().trim().parse().unwrap(),
        skip: hash.get("skip").unwrap().trim().parse().unwrap(),
        numass: hash.get("numass").unwrap().trim().parse().unwrap(),
        numrods: hash.get("numrods").unwrap().trim().parse().unwrap(),
        roddia: hash.get("roddia").unwrap().trim().parse().unwrap(),
        rodpitch: hash.get("rodpitch").unwrap().trim().parse().unwrap(),
        mpfr: hash.get("mpfr").unwrap().trim().parse().unwrap(),
        mpwr: hash.get("mpwr").unwrap().trim().parse().unwrap(),
        boundl: hash.get("boundl").unwrap().trim().parse().unwrap(),
        boundr: hash.get("boundr").unwrap().trim().parse().unwrap(),
    };

    let xsdata = XSData {
        sigtr: hash
            .get("sigtr")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        sigis: hash
            .get("sigis")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        sigds: hash
            .get("sigds")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        siga: hash
            .get("siga")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        sigf: hash
            .get("sigf")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        nut: hash
            .get("nut")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        chit: hash
            .get("chit")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
    };
    let matid: Vec<u8> = hash
        .get("matid")
        .unwrap()
        .split_ascii_whitespace()
        .map(|x| x.parse::<u8>().unwrap())
        .collect();
    (variables, xsdata, matid)
}

// fn mesh_gen(matid: Vec<u8>, mpfr: usize, mpwr: usize) -> Vec<u8> {
//     matid
//         .into_iter()
//         .flat_map(|x| {
//             if x == 0 || x == 1 {
//                 repeat(x).take(mpfr as usize)
//             } else {
//                 repeat(x).take(mpwr as usize)
//             }
//         })
//         .collect()
// }
fn main() {
    // let (variables, xsdata, matid) = process_input();

    // let meshid = mesh_gen(matid, variables.mpfr, variables.mpwr);

    // below is for timing
    let mut now = SystemTime::now();

    for zyn in 0..1000000 {
        if zyn % 10000 == 0 {
            print!(
                "Average time over those 10000 runs was {} microseconds \n",
                now.elapsed().unwrap().as_micros() / 10000
            );
            now = SystemTime::now();
        }
        let (variables, xsdata, matid) = process_input();
    }
}
