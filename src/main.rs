use csv::Writer;
use memmap2::MmapOptions;
use rand::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::iter::repeat;
use std::process::Command;
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
    inv_sigtr: Vec<f64>,
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
}

struct Tally {
    tally_length0: Vec<f64>,
    tally_length1: Vec<f64>,
}

struct SoltuionResults {
    flux0: Vec<f64>,
    flux1: Vec<f64>,
    fission_source: Vec<f64>,
    k: Vec<f64>,
}

// If multithreading the input processing
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
    let file = File::open("./SampleInputFile.txt").expect("Unable to read the file");
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
        inv_sigtr: hash
            .get("sigtr")
            .unwrap()
            .split_whitespace()
            .map(|x| x.parse::<f64>().unwrap().powi(-1))
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
        variables, xsdata, matid, deltax,
        // hash.get("solution").unwrap().trim().parse().unwrap(),
        1,
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
    // take advantage of code bitshifting instead of dividing in this case by dividing uan int by 2
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
    for item in temp {
        if item == 0 || item == 1 {
            let temp = Mesh {
                matid: item,
                delta_x: deltax.fuel,
                mesh_left: mesh_left,
                mesh_right: mesh_left + deltax.fuel,
            };
            mesh.push(temp);
            mesh_left += deltax.fuel;
        } else {
            let temp = Mesh {
                matid: item,
                delta_x: deltax.water,
                mesh_left: mesh_left,
                mesh_right: mesh_left + deltax.water,
            };
            mesh.push(temp);
            mesh_left += deltax.water;
        }
    }

    mesh
}

fn spawn_neutron(variables: &Variables, delta_x: &DeltaX) -> (usize, f64, f64, u8) {
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
    let chi = if temp_energy >= 0.7 { 0 } else { 1 };
    match variables.energygroups {
        4 => {
            return (
                ((variables.mpfr + variables.mpwr) * temp_int.unwrap().parse::<usize>().unwrap())
                    + variables.mpwr / 2
                    + temp2_int.unwrap().parse::<usize>().unwrap()
                    - 1, // I think I need to subtract one since vectors index @0
                ("0.".to_string() + temp2_dec.unwrap())
                    .parse::<f64>()
                    .unwrap(),
                2.0 * (thread_rng().gen::<f64>()) - 1.0,
                chi,
            );
        }
        _ => {
            return (
                ((variables.mpfr + variables.mpwr) * temp_int.unwrap().parse::<usize>().unwrap())
                    + variables.mpwr / 2
                    + temp2_int.unwrap().parse::<usize>().unwrap()
                    - 1,
                ("0.".to_string() + temp2_dec.unwrap())
                    .parse::<f64>()
                    .unwrap(),
                2.0 * (thread_rng().gen::<f64>()) - 1.0,
                0,
            )
        }
    }
}

fn tally_calc(
    y: usize,
    tally0: f64,
    tally1: f64,
    mesh_end: f64,
    start_x: f64,
    mu: f64,
    chi: u8,
) -> (f64, f64) {
    match chi {
        0 => return (tally0 + ((start_x - mesh_end) / mu).abs(), tally1),
        1 => return (tally0, tally1 + ((start_x - mesh_end) / mu).abs()),
        _ => (tally0, tally1),
    }
}

fn hit_boundary(mu: f64, start_x: f64, delta_s: f64, bound: f64, mesh_end: f64) -> (f64, f64, f64) {
    (
        mu * (-bound),
        (delta_s + (start_x - mesh_end)) * (-bound),
        mesh_end,
    )
}

fn cross_mesh(
    mesh_index: usize,
    mu: f64,
    start_x: f64,
    mesh_end: f64,
    delta_s: f64,
) -> (f64, f64, usize) {
    let mesh_index = if mu >= 0.0 {
        mesh_index + 1
    } else {
        mesh_index - 1
    };

    (delta_s + (start_x - mesh_end), mesh_end, mesh_index)
}

fn interaction(y: usize, xsdata: &XSData, xs_index: usize, chi: u8) -> (bool, u8, f64) {
    let interaction: f64 = thread_rng().gen();
    let absorption: f64 = xsdata.siga[xs_index] * xsdata.inv_sigtr[xs_index];
    let fission: f64 = absorption + xsdata.sigf[xs_index] * xsdata.inv_sigtr[xs_index];
    let in_scatter: f64 = fission + (xsdata.sigis[xs_index] * xsdata.inv_sigtr[xs_index]);
    if interaction < absorption {
        return (false, chi, 0.0);
    } else if interaction < fission {
        return (false, chi, 0.0);
    } else if interaction < in_scatter {
        return (true, chi, 2.0 * (thread_rng().gen::<f64>()) - 1.0);
    } else {
        //down scatter
        return (true, chi + 1, 2.0 * (thread_rng().gen::<f64>()) - 1.0);
    }
}

fn particle_travel(
    y: usize,
    mut tally: Tally,
    meshid: &Vec<Mesh>,
    mut mesh_index: usize,
    mut chi: u8,
    mut mu: f64,
    mut start_x: f64,
    mattypes: u8,
    boundr: f64,
    boundl: f64,
    xsdata: &XSData,
) -> (bool, Tally) {
    // determine particle travel length
    let s: f64 = -thread_rng().gen::<f64>().ln()
        * xsdata.inv_sigtr[(meshid[mesh_index].matid + mattypes * (chi)) as usize];
    let mut delta_s: f64 = mu * s;

    let mut same_material = true;
    while same_material == true {
        let end_x = start_x + delta_s;
        let mesh_end = if mu >= 0.0 {
            meshid[mesh_index].mesh_right
        } else {
            meshid[mesh_index].mesh_left
        };

        if mu * (end_x - mesh_end) >= 0.0 && (mesh_index == 0 || mesh_index == meshid.len() - 1) {
            // hit the boundary
            (
                tally.tally_length0[mesh_index],
                tally.tally_length1[mesh_index],
            ) = tally_calc(
                y,
                tally.tally_length0[mesh_index],
                tally.tally_length1[mesh_index],
                mesh_end,
                start_x,
                mu,
                chi,
            );
            let bound = if mu >= 0.0 { boundr } else { boundl };
            (mu, delta_s, start_x) = hit_boundary(mu, start_x, delta_s, bound, mesh_end);
        } else if mu * (end_x - mesh_end) >= 0.0 {
            // cross mesh boundary
            (
                tally.tally_length0[mesh_index],
                tally.tally_length1[mesh_index],
            ) = tally_calc(
                y,
                tally.tally_length0[mesh_index],
                tally.tally_length1[mesh_index],
                mesh_end,
                start_x,
                mu,
                chi,
            );

            let prev_mat = meshid[mesh_index].matid;

            (delta_s, start_x, mesh_index) = cross_mesh(mesh_index, mu, start_x, mesh_end, delta_s);

            if prev_mat != meshid[mesh_index].matid {
                same_material = false;
            }
        } else {
            // interaction
            (
                tally.tally_length0[mesh_index],
                tally.tally_length1[mesh_index],
            ) = tally_calc(
                y,
                tally.tally_length0[mesh_index],
                tally.tally_length1[mesh_index],
                end_x,
                start_x,
                mu,
                chi,
            );

            let xs_index: usize = (meshid[mesh_index].matid + mattypes * (chi)) as usize;
            let (particle_exists, _chi, _mu) = interaction(y, &xsdata, xs_index, chi);
            if particle_exists == false {
                return (particle_exists, tally);
            }
            chi = _chi;
            mu = _mu;
            delta_s = mu
                * -thread_rng().gen::<f64>().ln()
                * xsdata.inv_sigtr[(meshid[mesh_index].matid + mattypes * (chi)) as usize];
        }
    }
    (true, tally)
}

fn particle_lifetime(
    y: usize,
    mut tally: Tally,
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    variables: &Variables,
    delta_x: &DeltaX,
) -> Tally {
    // spawn_sub_mesh is the partial distance through the mesh
    let (mesh_index, spawn_sub_mesh, mu, chi) = spawn_neutron(&variables, &delta_x);
    let start_x: f64 = meshid[mesh_index].mesh_left + (spawn_sub_mesh * delta_x.fuel);

    let mut particle_exists = true;
    while particle_exists == true {
        (particle_exists, tally) = particle_travel(
            y,
            tally,
            meshid,
            mesh_index,
            chi,
            mu,
            start_x,
            variables.mattypes,
            variables.boundr,
            variables.boundl,
            xsdata,
        );
    }
    return tally;
}

fn monte_carlo(
    variables: &Variables,
    xsdata: &XSData,
    delta_x: &DeltaX,
    meshid: &Vec<Mesh>,
    mut k_new: f64,
) -> SoltuionResults {
    // println!("running monte carlo code");
    let mut results = SoltuionResults {
        flux0: vec![0.0; meshid.len()],
        flux1: vec![0.0; meshid.len()],
        fission_source: vec![0.0; meshid.len()],
        k: vec![0.0; variables.generations],
    };

    for x in 0..variables.generations {
        let mut tally = Tally {
            tally_length0: vec![0.0; meshid.len()],
            tally_length1: vec![0.0; meshid.len()],
        };
        let k = k_new;
        k_new = 0.0;
        // println!("Starting Generation {x}");
        for _y in 0..variables.histories {
            tally = particle_lifetime(_y, tally, xsdata, meshid, variables, delta_x)
        }

        let ((flux0, flux1), fission_source): ((Vec<f64>, Vec<f64>), Vec<f64>) =
            (0..tally.tally_length0.len())
                .into_iter()
                .map(|x| {
                    let delta_x = meshid[x].delta_x;
                    let matid = meshid[x].matid;
                    let flux0 = tally.tally_length0[x] / (k * variables.histories as f64 * delta_x);
                    let fission_source0 =
                        xsdata.nut[matid as usize] * xsdata.sigf[matid as usize] * flux0;
                    let flux1 = tally.tally_length1[x] / (k * variables.histories as f64 * delta_x);
                    let fission_source1 = xsdata.nut[(matid + variables.mattypes) as usize]
                        * xsdata.sigf[(matid + variables.mattypes) as usize]
                        * flux1;
                    let fission_source = fission_source1 + fission_source0;
                    k_new += k * delta_x * fission_source;
                    return ((flux0, flux1), fission_source);
                })
                .unzip();

        for index in 0..meshid.len() {
            results.flux0[index] += flux0[index];
            results.flux1[index] += flux1[index];
            results.fission_source[index] += fission_source[index];
        }

        results.k[x] = k;

        // results.flux0 = results
        //     .flux0
        //     .iter()
        //     .zip(&flux0)
        //     .map(|(a, b)| a + b)
        //     .collect();
        // results.flux1 = results
        //     .flux1
        //     .iter()
        //     .zip(&flux1)
        //     .map(|(a, b)| a + b)
        //     .collect();
        // results.fission_source = results
        //     .fission_source
        //     .iter()
        //     .zip(&fission_source)
        //     .map(|(a, b)| a + b)
        //     .collect();

        // println!("k is {}", k);
    }
    results
}

fn plot_solution(results: SoltuionResults) -> Result<(), Box<dyn Error>> {
    let output_k: Vec<String> = results.k.iter().map(|x| x.to_string()).collect();
    let output_flux0: Vec<String> = results.flux0.iter().map(|x| x.to_string()).collect();
    let output_flux1: Vec<String> = results.flux1.iter().map(|x| x.to_string()).collect();
    let output_fission: Vec<String> = results
        .fission_source
        .iter()
        .map(|x| x.to_string())
        .collect();

    let mut wtr_k = Writer::from_path("./k_eff.csv")?;

    wtr_k.write_record(&output_k)?;
    wtr_k.flush()?;

    let mut wtr = Writer::from_path("./interface.csv")?;

    wtr.write_record(&output_flux0)?;
    wtr.write_record(&output_flux1)?;
    wtr.write_record(&output_fission)?;
    wtr.flush()?;

    Command::new("python").arg("plot.py").spawn()?;

    Ok(())
}

fn main() {
    let (variables, xsdata, matid, deltax, solution) = process_input();

    // index meshid via [assembly# * rod#]
    let meshid = mesh_gen(matid, &variables, &deltax);

    let results = match solution {
        1 => monte_carlo(&variables, &xsdata, &deltax, &meshid, 1.0),
        _ => SoltuionResults {
            flux0: Vec::new(),
            flux1: Vec::new(),
            fission_source: Vec::new(),
            k: Vec::new(),
        }, // not implemented
    };

    let _ = plot_solution(results);

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