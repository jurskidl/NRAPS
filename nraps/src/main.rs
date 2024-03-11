use memmap2::MmapOptions;
use std::collections::HashMap;
use std::fs::File;
use std::iter::repeat;
use std::{result, thread};
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

fn next_end_line(mut end: usize, buffer: &[u8]) -> usize {
    while buffer[end] != NEWLINE && end < buffer.len() {
        end += 1;
    }
    end + 1
}

fn skip_line(mut pos: usize, end: usize, buffer: &[u8]) -> usize {
    while buffer[pos] != NEWLINE && pos < end {
        pos += 1;
    }
    pos
}

fn scan_ascii_chunk(start: usize, mut end: usize, buffer: &[u8]) -> HashMap<String, String> {
    let mut hash: HashMap<String, String> = HashMap::new();
    let mut pos = start;
    let mut line_start = start;
    let mut name_end = start;
    let mut val_start = start;
    if end != buffer.len() && buffer[end] != NEWLINE {
        end = next_end_line(end, buffer);
    }
    //println!("end was {old_end} is now {end} \n");
    while pos < end {
        match buffer[pos] {
            POUND => {
                pos = skip_line(pos, end, buffer);
                line_start = pos;
            }
            EQUALS => {
                if buffer[pos - 1] == SPACE {
                    name_end = pos - 1
                } else {
                    name_end = pos
                }

                if buffer[pos + 1] == SPACE {
                    val_start = pos + 1
                } else {
                    val_start = pos
                }
            }
            NEWLINE => {
                if name_end > line_start {
                    let key = String::from_utf8_lossy(&buffer[line_start..name_end])
                        .trim()
                        .to_string()
                        .to_ascii_lowercase();
                    // println!("{key}");
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

    // If I ever manage to figure out how to merge hashmaps this could be useful
    let size = mapped_file.len();
    let threads: usize = thread::available_parallelism().unwrap().get();
    let chunk_length = size / threads;
    let starting_points: Vec<usize> = (0..threads).map(|x| x * chunk_length).collect();
    let mut ending_points: Vec<usize> = Vec::from_iter(starting_points[1..threads].iter().cloned());
    ending_points.push(size);

    // Using a scoped pool to make it easy to share the immutable data from above.
    // Scan each segment to find station names and values.
    let mut result: HashMap<String, String> = HashMap::with_capacity(NUM_VARS);
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(threads);
        for thread in 0..threads {
            let start = starting_points[thread];
            let end = ending_points[thread];
            let buffer = &mapped_file;
            let handle = scope.spawn(move || scan_ascii_chunk(start, end, &buffer));
            handles.push(handle);
        }

        // Aggregate the results
        for handle in handles {
            let chunk_result = handle.join().unwrap();
            for (key, value) in chunk_result {
                result
                    .entry(key.trim().to_string())
                    .and_modify(|existing| *existing = existing.to_owned() + " " + &value)
                    .or_insert(value);
            }
        }
    });

    // let end = mapped_file.len();
    // let hash = scan_ascii_chunk(end, &&mapped_file);

    for (key, value) in &result {
        println!("{}: {}", key, value);
    }

    let variables = Variables {
        solution: result.get("solution").unwrap().trim().parse().unwrap(),
        analk: result.get("analk").unwrap().trim().parse().unwrap(),
        mattypes: result.get("mattypes").unwrap().trim().parse().unwrap(),
        energygroups: result.get("energygroups").unwrap().trim().parse().unwrap(),
        solver: match result.get("solver").unwrap().trim() {
            "1" => Solver::Gaussian,
            "2" => Solver::Jacobian,
            "3" => Solver::Sor,
            _ => Solver::LinAlg,
        },
        generations: result.get("generations").unwrap().trim().parse().unwrap(),
        histories: result.get("histories").unwrap().trim().parse().unwrap(),
        skip: result.get("skip").unwrap().trim().parse().unwrap(),
        numass: result.get("numass").unwrap().trim().parse().unwrap(),
        numrods: result.get("numrods").unwrap().trim().parse().unwrap(),
        roddia: result.get("roddia").unwrap().trim().parse().unwrap(),
        rodpitch: result.get("rodpitch").unwrap().trim().parse().unwrap(),
        mpfr: result.get("mpfr").unwrap().trim().parse().unwrap(),
        mpwr: result.get("mpwr").unwrap().trim().parse().unwrap(),
        boundl: result.get("boundl").unwrap().trim().parse().unwrap(),
        boundr: result.get("boundr").unwrap().trim().parse().unwrap(),
    };

    let xsdata = XSData {
        sigtr: result
            .get("sigtr")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        sigis: result
            .get("sigis")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        sigds: result
            .get("sigds")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        siga: result
            .get("siga")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        sigf: result
            .get("sigf")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        nut: result
            .get("nut")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
        chit: result
            .get("chit")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect(),
    };

    let matid: Vec<u8> = result
        .get("matid")
        .unwrap()
        .split_ascii_whitespace()
        .map(|x| x.parse::<u8>().unwrap())
        .collect();

    (variables, xsdata, matid)
}

fn mesh_gen(matid: Vec<u8>, mpfr: usize, mpwr: usize) -> Vec<u8> {
    matid
        .into_iter()
        .flat_map(|x| {
            if x == 0 || x == 1 {
                repeat(x).take(mpfr as usize)
            } else {
                repeat(x).take(mpwr as usize)
            }
        })
        .collect()
}
fn main() {
    let (variables, xsdata, matid) = process_input();

    let meshid = mesh_gen(matid, variables.mpfr, variables.mpwr);

    // below is for timing
    // let now = SystemTime::now();

    // for zyn in 0..1000000 {
    //     let (variables, xsdata, matid) = process_input();
    //     print!("{}\n", zyn)
    // }

    // print!("{}\n", now.elapsed().unwrap().as_millis());

    // sleep(Duration::new(30, 0));
}
