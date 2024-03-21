use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::iter::repeat;
use std::process::Command;
// For Multithreading
use std::thread;
// Use these for timing
use std::time::SystemTime;

use csv::Writer;
use memmap2::MmapOptions;
use rand::prelude::*;

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
    skip: usize,
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

struct SolutionResults {
    flux: Vec<Vec<f64>>,
    assembly_average: Vec<Vec<f64>>,
    fission_source: Vec<f64>,
    k: Vec<f64>,
    k_fund: Vec<f64>,
}

fn skip_line(mut pos: usize, end: usize, buffer: &[u8]) -> usize {
    while buffer[pos] != NEWLINE && pos < end {
        pos += 1;
    }
    pos + 1
}

fn scan_ascii_chunk(start: usize, end: usize, buffer: &[u8]) -> HashMap<String, String> {
    let mut hash: HashMap<String, String> = HashMap::with_capacity(24);

    let mut pos = start;
    let mut line_start = start;
    let mut name_end = start;
    let mut val_start = start;

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
                } else {}
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
        rodpitch: hash.get("rodpitch").unwrap().trim().parse::<f64>().unwrap()
            - hash.get("roddia").unwrap().trim().parse::<f64>().unwrap(),
        mpfr: hash.get("mpfr").unwrap().trim().parse().unwrap(),
        mpwr: hash.get("mpwr").unwrap().trim().parse().unwrap(),
        boundl: hash.get("boundl").unwrap().trim().parse().unwrap(),
        boundr: hash.get("boundr").unwrap().trim().parse().unwrap(),
    };

    let deltax = DeltaX {
        fuel: variables.roddia / variables.mpfr as f64,
        water: variables.rodpitch / variables.mpwr as f64,
    };

    // index into vectors via desired_xs = sigtr[(mat# + (energygroup*mattypes) as usize]
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

fn mesh_gen(matid: Vec<u8>, variables: &Variables, deltax: &DeltaX) -> (Vec<Mesh>, Vec<usize>) {
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

    // More generic version
    for index1 in 1..variables.numass {
        for _index2 in 0..variables.mpwr {
            temp.remove(((index1 * temp.len() as u8) / variables.numass) as usize);
        }
    }

    // MPWR/2 will be rounded down in the case of an odd int
    temp.drain(0..(variables.mpwr / 2));
    temp.truncate(temp.len() - (variables.mpwr / 2));

    let fuel_indices: Vec<usize> = temp
        .iter()
        .enumerate()
        .filter(|(_, &r)| r == 0 || r == 1)
        .map(|(index, _)| index)
        .collect::<Vec<_>>();

    let mut mesh: Vec<Mesh> = Vec::with_capacity(temp.len());

    let mut mesh_left: f64 = 0.0;
    for item in temp {
        if item == 0 || item == 1 {
            let temp = Mesh {
                matid: item,
                delta_x: deltax.fuel,
                mesh_left,
                mesh_right: mesh_left + deltax.fuel,
            };
            mesh.push(temp);
            mesh_left += deltax.fuel;
        } else {
            let temp = Mesh {
                matid: item,
                delta_x: deltax.water,
                mesh_left,
                mesh_right: mesh_left + deltax.water,
            };
            mesh.push(temp);
            mesh_left += deltax.water;
        }
    }

    (mesh, fuel_indices)
}

fn spawn_neutron(fuel_indices: &Vec<usize>, variables: &Variables, xsdata: &XSData, meshid: &Vec<Mesh>) -> (usize, f64, f64, u8) {
    let index = fuel_indices[thread_rng().gen_range(0..fuel_indices.len())];
    let chi: f64 = thread_rng().gen();
    let mut temp = 0.0;
    let chit: Vec<f64> = xsdata
        .chit
        .iter()
        .skip(meshid[index].matid as usize)
        .step_by(variables.mattypes as usize)
        .map(|x| {
            temp += x;
            return temp;
        })
        .collect();
    let neutron_energy = chit.iter().position(|&x| x > chi).unwrap() as u8;
    (
        index,
        thread_rng().gen(),
        2.0 * (thread_rng().gen::<f64>()) - 1.0,
        neutron_energy,
    )
}

fn hit_boundary(mu: f64, start_x: f64, delta_s: f64, bound: f64, mesh_end: f64) -> (f64, f64, f64) {
    (
        mu * (-bound),
        (delta_s + (start_x - mesh_end)) * (-bound),
        mesh_end,
    )
}

fn cross_mesh(mesh_index: usize, mu: f64, start_x: f64, mesh_end: f64, delta_s: f64) -> (f64, f64, usize) {
    let mesh_index = if mu >= 0.0 {
        mesh_index + 1
    } else {
        mesh_index - 1
    };

    (delta_s + (start_x - mesh_end), mesh_end, mesh_index)
}

fn interaction(xsdata: &XSData, xs_index: usize, neutron_energy: u8) -> (bool, u8, f64) {
    let interaction: f64 = thread_rng().gen();
    let absorption: f64 = xsdata.siga[xs_index] * xsdata.inv_sigtr[xs_index];
    let fission: f64 = xsdata.siga[xs_index] * xsdata.inv_sigtr[xs_index];
    let in_scatter: f64 = fission + (xsdata.sigis[xs_index] * xsdata.inv_sigtr[xs_index]);
    if interaction < absorption {
        return (false, neutron_energy, 0.0);
    } else if interaction < in_scatter {
        return (
            true,
            neutron_energy,
            2.0 * (thread_rng().gen::<f64>()) - 1.0,
        );
    } else {
        //down scatter
        return (
            true,
            neutron_energy + 1,
            2.0 * (thread_rng().gen::<f64>()) - 1.0,
        );
    }
}

fn particle_travel(
    mut tally: Vec<Vec<f64>>,
    meshid: &Vec<Mesh>,
    mut mesh_index: usize,
    mut neutron_energy: u8,
    mut mu: f64,
    mut start_x: f64,
    mattypes: u8,
    boundr: f64,
    boundl: f64,
    xsdata: &XSData,
) -> (bool, Vec<Vec<f64>>, usize, f64, u8, f64) {
    // determine particle travel length
    let s: f64 = -thread_rng().gen::<f64>().ln() * xsdata.inv_sigtr[(meshid[mesh_index].matid + (mattypes * neutron_energy)) as usize];
    // sum this and divide by number of particles, should be mean free path
    let mut delta_s: f64 = mu * s;

    let mut same_material = true;
    while same_material == true {
        let end_x = start_x + delta_s;
        let mesh_end = if mu >= 0.0 {
            meshid[mesh_index].mesh_right
        } else {
            meshid[mesh_index].mesh_left
        };

        if (mu < 0.0 && (mesh_end - start_x) > (end_x - start_x) && mesh_index == 0)
            || (mu >= 0.0
            && (end_x - start_x) > (mesh_end - start_x)
            && mesh_index == meshid.len() - 1)
        {
            // hit the boundary
            tally[neutron_energy as usize][mesh_index] += ((start_x - mesh_end) / mu).abs();
            let bound = if mu >= 0.0 { boundr } else { boundl };
            (mu, delta_s, start_x) = hit_boundary(mu, start_x, delta_s, bound, mesh_end);
        } else if (end_x - start_x).abs() > (mesh_end - start_x).abs() {
            // cross mesh boundary
            tally[neutron_energy as usize][mesh_index] += ((start_x - mesh_end) / mu).abs();

            let prev_mat = meshid[mesh_index].matid;

            (delta_s, start_x, mesh_index) = cross_mesh(mesh_index, mu, start_x, mesh_end, delta_s);

            if prev_mat != meshid[mesh_index].matid {
                same_material = false;
            }
        } else {
            // interaction
            tally[neutron_energy as usize][mesh_index] += ((start_x - end_x) / mu).abs();

            let xs_index: usize = (meshid[mesh_index].matid + (mattypes * neutron_energy)) as usize;
            let (particle_exists, _neutron_energy, _mu) =
                interaction(&xsdata, xs_index, neutron_energy);
            if particle_exists == false {
                return (
                    particle_exists,
                    tally,
                    mesh_index,
                    mu,
                    neutron_energy,
                    start_x,
                );
            }
            start_x = end_x;
            neutron_energy = _neutron_energy;
            mu = _mu;
            delta_s = mu
                * -thread_rng().gen::<f64>().ln()
                * xsdata.inv_sigtr
                [(meshid[mesh_index].matid + (mattypes * neutron_energy)) as usize];
        }
    }
    (true, tally, mesh_index, mu, neutron_energy, start_x)
}

fn particle_lifetime(
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    fuel_indices: &Vec<usize>,
    variables: &Variables,
    delta_x: &DeltaX,
    start: usize,
    end: usize,
) -> Vec<Vec<f64>> {
    let mut tally: Vec<Vec<f64>> = vec![vec![0.0; meshid.len()]; variables.energygroups as usize];

    for _y in start..end {
        // spawn_sub_mesh is the partial distance through the mesh
        let (mut mesh_index, spawn_sub_mesh, mut mu, mut neutron_energy) =
            spawn_neutron(fuel_indices, &variables, &xsdata, &meshid);
        let mut start_x: f64 = meshid[mesh_index].mesh_left + (spawn_sub_mesh * delta_x.fuel);

        let mut particle_exists: bool = true;
        while particle_exists == true {
            (
                particle_exists,
                tally,
                mesh_index,
                mu,
                neutron_energy,
                start_x,
            ) = particle_travel(
                tally,
                meshid,
                mesh_index,
                neutron_energy,
                mu,
                start_x,
                variables.mattypes,
                variables.boundr,
                variables.boundl,
                xsdata,
            );
        }
    }
    return tally;
}

fn monte_carlo(
    threads: usize,
    variables: &Variables,
    xsdata: &XSData,
    delta_x: &DeltaX,
    meshid: &Vec<Mesh>,
    fuel_indices: &Vec<usize>,
    mut k_new: f64,
) -> SolutionResults {
    // println!("running monte carlo code");
    let mut results = SolutionResults {
        flux: vec![vec![0.0; meshid.len()]; variables.energygroups as usize],
        assembly_average: vec![vec![0.0; meshid.len()]; variables.energygroups as usize],
        fission_source: vec![0.0; meshid.len()],
        k: vec![0.0; variables.generations],
        k_fund: vec![0.0; variables.generations],
    };

    let power_level = 3565e6; // J/s
    let energy_fission = 200e6 * 1.602176634e-19; // J

    // A_core = N_rods * pitch^2, where N_rods = 193 assemblies * 264 rods / assembly and pitch is 1.26
    let diameter: f64 = 2.0 * ((50952.0 * variables.rodpitch) / std::f64::consts::PI).sqrt();
    let core_slice: f64 = 365.76 * diameter;
    let conversion: f64 = power_level / (energy_fission * core_slice);
    let fund: f64 = 1.0 / (variables.generations - (variables.skip - 1)) as f64;

    for x in 0..variables.generations {
        // let mut tally: Vec<Vec<f64>> =
        //     vec![vec![0.0; meshid.len()]; variables.energygroups as usize];
        let k = k_new;
        k_new = 0.0;

        let tally = particle_lifetime(
            xsdata,
            meshid,
            fuel_indices,
            variables,
            delta_x,
            0,
            variables.histories,
        );

        // For Multithreading
        // let threads: usize = 14; //thread::available_parallelism().unwrap().get();
        // let threaded_histories = variables.histories / threads;
        // let starting_points: Vec<usize> = (0..threads).map(|x| x * threaded_histories).collect();
        // let mut ending_points: Vec<usize> = Vec::from_iter(starting_points[1..threads].iter().cloned());
        // ending_points.push(variables.histories);
        //
        // thread::scope(|scope| {
        //     let mut tallies = Vec::with_capacity(threads);
        //     for thread in 0..threads {
        //         let start = starting_points[thread];
        //         let end = ending_points[thread];
        //         let tallied: thread::ScopedJoinHandle<'_, Vec<Vec<f64>>> = scope.spawn(move || {
        //             {
        //                 particle_lifetime(xsdata, meshid, fuel_indices, variables, delta_x, start, end)
        //             }
        //         });
        //         tallies.push(tallied);
        //     }
        //
        //     // Aggregate the results
        //     for tallied in tallies {
        //         let thread_result = tallied.join().unwrap();
        //         for energy in 0..variables.energygroups as usize {
        //             for index in 0..thread_result[energy].len() as usize {
        //                 tally[energy][index] += thread_result[energy][index]
        //             }
        //         }
        //     }
        // });

        for energy in 0..variables.energygroups as usize {
            for index in 0..tally[energy].len() {
                let delta_x = meshid[index].delta_x;
                let matid = meshid[index].matid;
                let flux = tally[energy][index] / (k * variables.histories as f64 * delta_x);
                let fission_source = xsdata.nut
                    [(matid + (variables.mattypes * energy as u8)) as usize]
                    * xsdata.sigf[(matid + (variables.mattypes * energy as u8)) as usize]
                    * flux;
                k_new += k * delta_x * fission_source;
                if x >= variables.skip {
                    results.flux[energy][index] += flux * conversion;
                    results.fission_source[index] += fission_source * fund;
                }
            }
        }
        results.k[x] = k_new;
    }

    for energy in 0..variables.energygroups as usize {
        let average_ass_1 = (0..results.flux[energy].len() / 2)
            .into_iter()
            .map(|x| results.flux[energy][x])
            .sum::<f64>()
            / (results.flux[energy].len() / 2) as f64;
        let average_ass_2 = (results.flux[energy].len() / 2..results.flux[energy].len())
            .into_iter()
            .map(|x| results.flux[energy][x])
            .sum::<f64>()
            / (results.flux[energy].len() / 2) as f64;

        for index in 0..results.flux[energy].len() / 2 {
            results.assembly_average[energy][index] = average_ass_1;
        }
        for index in results.flux[energy].len() / 2..results.flux[energy].len() {
            results.assembly_average[energy][index] = average_ass_2;
        }
    }

    results.k_fund[variables.skip] = results.k[variables.skip];

    for generations in (variables.skip + 1)..variables.generations {
        results.k_fund[generations] = (variables.skip..=generations)
            .into_iter()
            .map(|x| results.k[x])
            .sum::<f64>()
            / (generations - (variables.skip - 1)) as f64;
    }

    results
}

fn plot_solution(
    results: SolutionResults,
    energygroups: u8,
    generations: usize,
    number_meshes: usize,
    assembly_length: f64,
) -> Result<(), Box<dyn Error>> {
    let output_k_fund: Vec<String> = results.k_fund.iter().map(|x| x.to_string()).collect();
    let output_k: Vec<String> = results.k.iter().map(|x| x.to_string()).collect();
    let mut output_flux =
        vec![vec!["0.0".to_string(); results.fission_source.len()]; energygroups as usize];
    for energy in 0..energygroups as usize {
        for index in 0..results.fission_source.len() {
            output_flux[energy][index] = results.flux[energy][index].to_string();
        }
    }
    let mut output_average =
        vec![vec!["0.0".to_string(); results.fission_source.len()]; energygroups as usize];
    for energy in 0..energygroups as usize {
        for index in 0..results.fission_source.len() {
            output_average[energy][index] = results.assembly_average[energy][index].to_string();
        }
    }
    let output_fission: Vec<String> = results
        .fission_source
        .iter()
        .map(|x| x.to_string())
        .collect();

    let _handle1 = thread::spawn(move || {
        let mut wtr_vars = Writer::from_path("./vars.csv").expect("failed to open vars.csv");

        wtr_vars
            .write_record([&assembly_length.to_string()])
            .expect("failed to write to vars.csv");
        wtr_vars
            .write_record([&number_meshes.to_string()])
            .expect("failed to write to vars.csv");
        wtr_vars
            .write_record([&generations.to_string()])
            .expect("failed to write to vars.csv");
        wtr_vars.flush().expect("failed to clear buffer");
    });

    let _handle2 = thread::spawn(move || {
        let mut wtr_k = Writer::from_path("./k_eff.csv").expect("failed to open k_eff.csv");

        wtr_k
            .write_record(&output_k)
            .expect("failed to write to k_eff.csv");
        wtr_k
            .write_record(&output_k_fund)
            .expect("failed to write to k_eff.csv");
        wtr_k.flush().expect("failed to clear buffer");
    });

    let _handle3 = thread::spawn(move || {
        let mut wtr = Writer::from_path("./interface.csv").expect("failed to open interface.csv");

        for energy in 0..energygroups as usize {
            wtr.write_record(&output_flux[energy])
                .expect("failed write to interface.csv");
        }
        for energy in 0..energygroups as usize {
            wtr.write_record(&output_average[energy])
                .expect("failed write to interface.csv");
        }
        wtr.write_record(&output_fission)
            .expect("failed write to interface.csv");
        wtr.flush().expect("failed to clear buffer");
    });

    Command::new("python").arg("plot.py").spawn()?;

    Ok(())
}

fn main() {
    // let now = SystemTime::now();

    let (variables, xsdata, matid, deltax, solution) = process_input();

    let (meshid, fuel_indices) = mesh_gen(matid, &variables, &deltax);

    // let results = match solution {
    //     1 => monte_carlo(&variables, &xsdata, &deltax, &meshid, &fuel_indices, 1.0),
    //     _ => SoltuionResults {
    //         flux: Vec::new(),
    //         assembly_average: Vec::new(),
    //         fission_source: Vec::new(),
    //         k: Vec::new(),
    //         k_fund: Vec::new(),
    //     }, // not implemented
    // };
    //
    // let _ = plot_solution(
    //     results,
    //     variables.energygroups,
    //     variables.generations,
    //     meshid.len(),
    //     meshid[meshid.len() - 1].mesh_right,
    // );
    //
    // print!(
    //     "Run was completed in {} milliseconds \n",
    //     now.elapsed().unwrap().as_micros()
    // );

    // below is for timing
    let mut now = SystemTime::now();

    // for threads in 1..17 {
    for zyn in 0..11 {
        if zyn % 10 == 0 {
            print!(
                "Average time for 1 was {} microseconds \n",
                // threads,
                now.elapsed().unwrap().as_micros() / 10
            );
            now = SystemTime::now();
        }
        let results = match solution {
            1 => monte_carlo(
                1,// threads,
                &variables,
                &xsdata,
                &deltax,
                &meshid,
                &fuel_indices,
                1.0,
            ),
            _ => SolutionResults {
                flux: Vec::new(),
                assembly_average: Vec::new(),
                fission_source: Vec::new(),
                k: Vec::new(),
                k_fund: Vec::new(),
            }, // not implemented
        };
    }
    // }
}
