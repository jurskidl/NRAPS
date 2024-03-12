use memmap2::MmapOptions;
use rand::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::iter::repeat;
// For Multithreading
// use std::thread;
// Use these for timing
// use std::time::SystemTime;

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
    analk: u8,
    mattypes: u8,
    energygroups: u8,
    solver: Solver,
    generations: usize,
    histories: usize,
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

struct DeltaX {
    fuel: f64,
    water: f64,
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

struct Mesh {
    matid: u8,
    delta_x: f64,
    mesh_left: f64,
    mesh_right: f64,
    mat_left: f64,
    mat_right: f64,
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

fn process_input() -> (Variables, XSData, Vec<u8>, DeltaX, u8) {
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

    let deltax = DeltaX {
        fuel: variables.roddia / variables.mpfr as f64,
        water: variables.rodpitch / variables.mpwr as f64,
    };

    // index into vectors via desired_xs = sigtr[(mat# + ((energygroup-1)*mattypes) as usize]
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

    // index vector via mat# = matid[config#*numass*((2*numrods)+1)]
    let matid: Vec<u8> = hash
        .get("matid")
        .unwrap()
        .split_ascii_whitespace()
        .map(|x| x.parse::<u8>().unwrap())
        .collect();
    (
        variables,
        xsdata,
        matid,
        deltax,
        hash.get("solution").unwrap().trim().parse().unwrap(),
    )
}

fn mesh_gen(matid: Vec<u8>, variables: &Variables, deltax: &DeltaX) -> Vec<Mesh> {
    let mut temp: Vec<u8> = matid
        .into_iter()
        .flat_map(|x| {
            if x == 0 || x == 1 {
                repeat(x).take(variables.mpfr as usize)
            } else {
                repeat(x).take(variables.mpwr as usize)
            }
        })
        .collect();

    // We want to remove the extra water meshes from the center of the array
    // take advantage of code bitshifting instead of dividing in this case by dividing an int by 2
    for _index in 0..variables.mpwr {
        temp.remove(temp.len() / 2);
    }

    // More generic version
    // for index1 in 0..numass - 1 {
    //     for index2 in 0..mpwr {
    //         temp.remove((index1 * temp.len()) / numass);
    //     }
    // }

    // MPWR/2 will be rounded down in the case of an odd int
    temp.drain(0..(variables.mpwr / 2));
    temp.truncate(temp.len() - (variables.mpwr / 2));

    let mut mesh: Vec<Mesh> = Vec::with_capacity(temp.len());

    let mut mesh_left: f64 = 0.0;
    let mut mat_left: f64 = 0.0;
    let mut mat_right: f64 = variables.rodpitch * 0.5;
    let mut last_item = temp[0];
    for item in temp {
        if (item == 0 || item == 1) && item != last_item {
            mat_left = mesh_left;
            mat_right = mat_left + variables.mpfr as f64 * deltax.fuel;
            for _x in 0..variables.mpfr {
                let temp = Mesh {
                    matid: item,
                    delta_x: deltax.fuel,
                    mesh_left: mesh_left,
                    mesh_right: mesh_left + deltax.fuel,
                    mat_left: mat_left,
                    mat_right: mat_right,
                };
                mesh.push(temp);
                mesh_left += deltax.fuel;
                last_item = item
            }
        } else if (item == 0 || item == 1) && item == last_item {
            for _x in 0..variables.mpfr {
                let temp = Mesh {
                    matid: item,
                    delta_x: deltax.fuel,
                    mesh_left: mesh_left,
                    mesh_right: mesh_left + deltax.fuel,
                    mat_left: mat_left,
                    mat_right: mat_right,
                };
                mesh.push(temp);
                mesh_left += deltax.fuel;
            }
        } else if item != last_item {
            mat_left = mesh_left;
            mat_right = mat_left + variables.mpwr as f64 * deltax.water;
            for _x in 0..variables.mpwr {
                let temp = Mesh {
                    matid: item,
                    delta_x: deltax.water,
                    mesh_left: mesh_left,
                    mesh_right: mesh_left + deltax.water,
                    mat_left: mat_left,
                    mat_right: mat_right,
                };
                mesh.push(temp);
                mesh_left += deltax.water;
                last_item = item
            }
        } else {
            for _x in 0..variables.mpwr {
                let temp = Mesh {
                    matid: item,
                    delta_x: deltax.water,
                    mesh_left: mesh_left,
                    mesh_right: mesh_left + deltax.water,
                    mat_left: mat_left,
                    mat_right: mat_right,
                };
                mesh.push(temp);
                mesh_left += deltax.water;
            }
        }
    }

    mesh
}

fn spawn_neutron(variables: &Variables, delta_x: &DeltaX) -> (usize, f64, u8) {
    let spawn_string: String = thread_rng()
        .gen_range(1.0..(variables.numass * variables.numrods) as f64)
        .to_string();
    let (temp_int, temp_dec) = spawn_string.split_once(".").unzip();
    let temp2 = (("0.".to_string() + temp_dec.unwrap())
        .parse::<f64>()
        .unwrap()
        / delta_x.fuel)
        .to_string();
    let (temp2_int, temp2_dec) = temp2.split_once(".").unzip();
    let temp_energy: f64 = thread_rng().gen();
    let chi = if temp_energy >= 0.7 { 1 } else { 2 };
    match variables.energygroups {
        4 => {
            return (
                ((variables.mpfr + variables.mpwr) * temp_int.unwrap().parse::<usize>().unwrap())
                    + variables.mpwr / 2
                    + temp2_int.unwrap().parse::<usize>().unwrap(),
                ("0.".to_string() + temp2_dec.unwrap())
                    .parse::<f64>()
                    .unwrap(),
                chi,
            )
        }
        _ => {
            return (
                ((variables.mpfr + variables.mpwr) * temp_int.unwrap().parse::<usize>().unwrap())
                    + variables.mpwr / 2
                    + temp2_int.unwrap().parse::<usize>().unwrap(),
                ("0.".to_string() + temp2_dec.unwrap())
                    .parse::<f64>()
                    .unwrap(),
                chi,
            )
        }
    }
}

fn monte_carlo(variables: &Variables, xsdata: &XSData, delta_x: &DeltaX, meshid: &Vec<Mesh>) {
    println!("ruunning monte carlo code");
    for x in 0..variables.generations {
        println!("Starting Generation {x}");
        let mut next_gen: f64 = 0.0;
        for _y in 0..variables.histories {
            // spawn_sub_mesh is the partial distance through the mesh
            let (mut mesh_index, spawn_sub_mesh, mut chi): (usize, f64, u8) =
                spawn_neutron(&variables, &delta_x);
            let mut mu: f64 = 2.0 * (thread_rng().gen::<f64>()) - 1.0;
            let mut particle_exists: bool = true;

            while particle_exists == true {
                let s: f64 = thread_rng().gen::<f64>().ln()
                    / xsdata.sigtr
                        [(meshid[mesh_index].matid + variables.mattypes * (chi - 1)) as usize];
                let mut delta_s: f64 = mu * s;

                let mut same_material: bool = true;
                let mut start_x: f64 = meshid[mesh_index].mesh_left + spawn_sub_mesh;

                let end_x = start_x + delta_s;
                while same_material == true {
                    if (mu < 0.0 && mesh_index == 0 && end_x < meshid[mesh_index].mesh_left)
                        || (mu >= 0.0
                            && mesh_index == meshid.len() - 1
                            && end_x > meshid[mesh_index].mesh_right)
                    {
                        // println!("Hitting the boundary");
                        if mu < 0.0 {
                            // tally_length[mesh_index] +=
                            //     (meshid[mesh_index].mesh_left - start_x).abs() / mu;
                            mu *= -variables.boundl;
                            delta_s *= -variables.boundl;
                            start_x = meshid[mesh_index].mesh_left;
                        } else {
                            // tally_length[mesh_index] +=
                            //     (meshid[mesh_index].mesh_right - start_x) / mu;
                            mu *= -variables.boundr;
                            delta_s *= -variables.boundr;
                            start_x = meshid[mesh_index].mesh_right;
                        }
                    } else if mu < 0.0 && end_x < meshid[mesh_index].mat_left {
                        // println!("exiting mat left");
                        // tally_length[mesh_index] +=
                        //     (meshid[mesh_index].mesh_left - start_x) / mu;
                        mesh_index -= 1;
                        while meshid[mesh_index].matid != meshid[mesh_index - 1].matid {
                            // tally_length[mesh_index] +=
                            //     (meshid[mesh_index].mesh_left - meshid[meshindex].mesh_right) / mu;
                            if mesh_index == 1 {
                                break;
                            }
                            mesh_index -= 1;
                        }
                        same_material = false;
                    } else if mu < 0.0 && end_x < meshid[mesh_index].mesh_left {
                        // println!("exiting mesh left");
                        // tally_length[mesh_index] +=
                        //     (meshid[mesh_index].mesh_left - start_x) / mu;
                        delta_s -= start_x - meshid[mesh_index].mesh_left;
                        start_x = meshid[mesh_index].mesh_left;
                        mesh_index -= 1;
                    } else if mu < 0.0 {
                        // println!("interacting left");
                        // tally_length[mesh_index] +=
                        //     (meshid[mesh_index].mesh_left - start_x) / mu
                        let xs_index: usize =
                            (meshid[mesh_index].matid + variables.mattypes * (chi - 1)) as usize;
                        let interaction = thread_rng().gen_range(0.0..xsdata.sigtr[xs_index]);
                        if interaction <= xsdata.sigis[xs_index] {
                            // println!("in scatter");
                            mu = 2.0 * (thread_rng().gen::<f64>()) - 1.0;
                        } else if interaction <= xsdata.sigis[xs_index] + xsdata.sigds[xs_index] {
                            // println!("down scatter");
                            chi += 1;
                            mu = 2.0 * (thread_rng().gen::<f64>()) - 1.0;
                        } else if interaction
                            <= xsdata.sigis[xs_index]
                                + xsdata.sigds[xs_index]
                                + xsdata.siga[xs_index]
                        {
                            // println!("absorption");
                            particle_exists = false;
                            same_material = false;
                        } else {
                            // println!("fission");
                            next_gen += xsdata.nut[xs_index];
                            particle_exists = false;
                            same_material = false;
                        }
                    } else if mu >= 0.0 && end_x > meshid[mesh_index].mat_right {
                        // println!(" exiting mat right");
                        // tally_length[mesh_index] +=
                        //     (meshid[mesh_index].mesh_right - start_x) / mu;
                        mesh_index += 1;
                        while meshid[mesh_index - 1].matid == meshid[mesh_index].matid
                            && mesh_index < (meshid.len() - 1)
                        {
                            // tally_length[mesh_index] +=
                            //     (meshid[mesh_index].mesh_right - meshid[meshindex].mesh_left) / mu;
                            mesh_index += 1;
                        }

                        if mesh_index != meshid.len() - 1 {
                            same_material = false
                        }
                    } else if mu >= 0.0 && end_x > meshid[mesh_index].mesh_right {
                        // println!("exiting mesh right");
                        // tally_length[mesh_index] +=
                        //     (meshid[mesh_index].mesh_right - start_x) / mu;
                        delta_s += meshid[mesh_index].mesh_right - start_x;
                        start_x = meshid[mesh_index].mesh_right;
                        mesh_index += 1;
                    } else {
                        // println!("interacting right");
                        // tally_length[mesh_index] +=
                        //     (meshid[mesh_index].mesh_right - start_x) / mu;
                        let xs_index: usize =
                            (meshid[mesh_index].matid + variables.mattypes * (chi - 1)) as usize;
                        let interaction = thread_rng().gen_range(0.0..xsdata.sigtr[xs_index]);
                        if interaction <= xsdata.sigis[xs_index] {
                            // println!("in scatter");
                            mu = 2.0 * (thread_rng().gen::<f64>()) - 1.0;
                        } else if interaction <= xsdata.sigis[xs_index] + xsdata.sigds[xs_index] {
                            // println!("down scatter");
                            chi += 1;
                            mu = 2.0 * (thread_rng().gen::<f64>()) - 1.0;
                        } else if interaction
                            <= xsdata.sigis[xs_index]
                                + xsdata.sigds[xs_index]
                                + xsdata.siga[xs_index]
                        {
                            // println!("absorption");
                            particle_exists = false;
                            same_material = false;
                        } else {
                            // println!("fission!");
                            next_gen += xsdata.nut[xs_index];
                            particle_exists = false;
                            same_material = false;
                        }
                    }
                }
            }
        }

        println!("keff would be {}", next_gen / variables.histories as f64);
    }
}
fn main() {
    let (variables, xsdata, matid, deltax, solution) = process_input();

    // index meshid via [assembly# * rod#]
    let meshid = mesh_gen(matid, &variables, &deltax);

    match solution {
        1 => monte_carlo(&variables, &xsdata, &deltax, &meshid),
        _ => {} // not implemented
    }

    // below is for timing
    // let mut now = SystemTime::now();

    // for zyn in 0..1000000 {
    //     if zyn % 10000 == 0 {
    //         print!(
    //             "Average time over those 10000 runs was {} microseconds \n",
    //             now.elapsed().unwrap().as_micros() / 10000
    //         );
    //         now = SystemTime::now();
    //     }
    //     let (variables, xsdata, matid) = process_input();
    //     let meshid = mesh_gen(matid, variables.mpfr, variables.mpwr);
    // }
}
