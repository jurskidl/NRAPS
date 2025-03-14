use std::iter::repeat;
// Use these for timing
use std::time::SystemTime;

use crate::discrete::nalgebra_method;
use crate::mc_code::monte_carlo;
use crate::plot_solution::plot_solution;
use crate::process_input::process_input;

mod discrete;
mod mc_code;
mod plot_solution;
mod process_input;

pub enum Solver {
    LinAlg,
    Gaussian,
    Jacobian,
    SR,
}

struct Variables {
    analk: u8,          // 1 byte
    mattypes: u8,       // 1 byte
    energygroups: u8,   // 1 byte
    generations: usize, // 8 bytes
    histories: usize,   // 8 bytes
    skip: usize,        // 8 bytes
    numass: u8,         // 1 byte
    numrods: u8,        // 1 byte
    roddia: f32,        // 4 bytes
    rodpitch: f32,      // 4 bytes
    mpfr: usize,        // 8 bytes
    mpwr: usize,        // 8 bytes
    boundl: f32,        // 4 bytes
    boundr: f32,        // 4 bytes
                        // 61 bytes used -> Allocates 64 bytes
                        // 3 bytes wasted
}

struct DeltaX {
    fuel: f32,
    water: f32,
}

struct XSData {
    sigt: Vec<f32>,
    sigs: Vec<f32>,
    mu: Vec<f32>,
    siga: Vec<f32>,
    sigf: Vec<f32>,
    nut: Vec<f32>,
    chit: Vec<f32>,
    scat_matrix: Vec<f32>,
    inv_sigtr: Vec<f32>,
}

struct CollapsedXsdata {
    sigt: f32,
    sigs: f32,
    siga: f32,
    sigf: f32,
    mu: f32,
    nu: f32,
    chi: f32,
}

struct Mesh {
    matid: u8,      // 1 byte
    delta_x: f32,   // 4 bytes
    mesh_left: f32, // 4 bytes
    mesh_right: f32, // 4 bytes
                    // 13 bytes used -> 16 bytes allocated
                    // 3 bytes wasted
}

struct SolutionResults {
    flux: Vec<Vec<f32>>,
    assembly_average: Vec<Vec<f32>>,
    fission_source: Vec<f32>,
    k: Vec<f32>,
    k_fund: Vec<f32>,
}

fn mesh_gen(matid: Vec<u8>, variables: &Variables, deltax: &DeltaX) -> (Vec<Mesh>, Vec<usize>) {
    let mut temp: Vec<u8> = matid
        .into_iter()
        .flat_map(|x| {
            if x == 0 || x == 1 {
                repeat(x).take(variables.mpfr)
            } else {
                repeat(x).take(variables.mpwr)
            }
        })
        .collect();

    // We want to remove the extra water meshes from the center of the array
    // This assumes the assembly has water rods in the middle
    for index1 in 1..variables.numass {
        for _index2 in 0..variables.mpwr {
            temp.remove((index1 as usize * temp.len()) / variables.numass as usize);
        }
    }

    // This assumes the assembly starts and ends with a water rod
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

    let mut mesh_left: f32 = 0.0;
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

// fn energy_collapse(
//     meshid: &Vec<Mesh>,
//     mattypes: u8,
//     energygroups: u8,
//     numrods: u8,
//     mpfr: usize,
//     flux: Vec<Vec<f32>>,
//     xsdata: &XSData,
// ) -> () {
//     let mut xs_02: Vec<CollapsedXsdata> = Vec::with_capacity(meshid.len());
//     let mut xs_24: Vec<CollapsedXsdata> = Vec::with_capacity(meshid.len());
//     let mut xs_04: Vec<CollapsedXsdata> = Vec::with_capacity(meshid.len());

//     for index in 0..meshid.len() {
//         let (
//             mut sigf_02,
//             mut sigf_24,
//             mut sigf_04,
//             mut siga_02,
//             mut siga_24,
//             mut siga_04,
//             mut sigs_02,
//             mut sigs_24,
//             mut sigs_04,
//             mut mu_02,
//             mut mu_24,
//             mut mu_04,
//             mut nu_02,
//             mut nu_24,
//             mut nu_04,
//             mut chi_02,
//             mut chi_24,
//             mut chi_04,
//         ) = (
//             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//             0.0,
//         );

//         for energy in 0..energygroups / 2 {
//             sigf_02 += xsdata.sigf[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             sigf_24 += xsdata.sigf[(meshid[index].matid + (mattypes * (energy + 2))) as usize]
//                 * flux[energy as usize][index];
//             siga_02 += xsdata.siga[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             siga_24 += xsdata.siga[(meshid[index].matid + (mattypes * (energy + 2))) as usize]
//                 * flux[energy as usize][index];
//             sigs_02 += xsdata.sigs[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             sigs_24 += xsdata.sigs[(meshid[index].matid + (mattypes * (energy + 2))) as usize]
//                 * flux[energy as usize][index];
//             mu_02 += xsdata.mu[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             mu_24 += xsdata.mu[(meshid[index].matid + (mattypes * (energy + 2))) as usize]
//                 * flux[energy as usize][index];
//             nu_02 += xsdata.nut[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             nu_24 += xsdata.nut[(meshid[index].matid + (mattypes * (energy + 2))) as usize]
//                 * flux[energy as usize][index];
//             chi_02 += xsdata.chit[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             chi_24 += xsdata.chit[(meshid[index].matid + (mattypes * (energy + 2))) as usize]
//                 * flux[energy as usize][index];
//         }

//         for energy in 0..energygroups {
//             sigf_04 += xsdata.sigf[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             siga_04 += xsdata.siga[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             sigs_04 += xsdata.sigs[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             mu_04 += xsdata.mu[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             nu_04 += xsdata.nut[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//             chi_04 += xsdata.chit[(meshid[index].matid + (mattypes * energy)) as usize]
//                 * flux[energy as usize][index];
//         }

//         let flux_cell_02: f32 = (0..energygroups / 2)
//             .into_iter()
//             .map(|energy| flux[energy as usize][index])
//             .sum();
//         let flux_cell_24: f32 = (energygroups / 2..energygroups)
//             .into_iter()
//             .map(|energy| flux[energy as usize][index])
//             .sum();
//         let flux_cell_04: f32 = flux_cell_02 + flux_cell_24;

//         sigf_02 = sigf_02 / flux_cell_02;
//         siga_02 = siga_02 / flux_cell_02;
//         sigs_02 = sigs_02 / flux_cell_02;
//         mu_02 = mu_02 / flux_cell_02;
//         nu_02 = nu_02 / flux_cell_02;
//         chi_02 = chi_02 / flux_cell_02;

//         sigf_24 = sigf_24 / flux_cell_24;
//         siga_24 = siga_24 / flux_cell_24;
//         sigs_24 = sigs_24 / flux_cell_24;
//         mu_24 = mu_24 / flux_cell_24;
//         nu_24 = nu_24 / flux_cell_24;
//         chi_24 = chi_24 / flux_cell_24;

//         sigf_04 = sigf_04 / flux_cell_04;
//         siga_04 = siga_04 / flux_cell_04;
//         sigs_04 = sigs_04 / flux_cell_04;
//         mu_04 = mu_04 / flux_cell_04;
//         nu_04 = nu_04 / flux_cell_04;
//         chi_04 = chi_04 / flux_cell_04;

//         xs_02.push(CollapsedXsdata {
//             sigt: sigf_02 + siga_02 + sigs_02,
//             sigs: sigs_02,
//             siga: siga_02,
//             sigf: sigf_02,
//             mu: mu_02,
//             nu: nu_02,
//             chi: chi_02,
//         });

//         xs_24.push(CollapsedXsdata {
//             sigt: sigf_24 + siga_24 + sigs_24,
//             sigs: sigs_24,
//             siga: siga_24,
//             sigf: sigf_24,
//             mu: mu_24,
//             nu: nu_24,
//             chi: chi_24,
//         });

//         xs_04.push(CollapsedXsdata {
//             sigt: sigf_04 + siga_04 + sigs_04,
//             sigs: sigs_04,
//             siga: siga_04,
//             sigf: sigf_04,
//             mu: mu_04,
//             nu: nu_04,
//             chi: chi_04,
//         })
//     }

//     let sigt_02 = xs_02.iter().map(|x| x.sigt).sum::<f32>() / meshid.len() as f32;
//     let sigt_24 = xs_24.iter().map(|x| x.sigt).sum::<f32>() / meshid.len() as f32;
//     let sigt_04 = xs_04.iter().map(|x| x.sigt).sum::<f32>() / meshid.len() as f32;

//     let siga_02 = xs_02.iter().map(|x| x.siga).sum::<f32>() / meshid.len() as f32;
//     let siga_24 = xs_24.iter().map(|x| x.siga).sum::<f32>() / meshid.len() as f32;
//     let siga_04 = xs_04.iter().map(|x| x.siga).sum::<f32>() / meshid.len() as f32;

//     let sigf_02 = xs_02.iter().map(|x| x.sigf).sum::<f32>() / meshid.len() as f32;
//     let sigf_24 = xs_24.iter().map(|x| x.sigf).sum::<f32>() / meshid.len() as f32;
//     let sigf_04 = xs_04.iter().map(|x| x.sigf).sum::<f32>() / meshid.len() as f32;

//     let sigs_02 = xs_02.iter().map(|x| x.sigs).sum::<f32>() / meshid.len() as f32;
//     let sigs_24 = xs_24.iter().map(|x| x.sigs).sum::<f32>() / meshid.len() as f32;
//     let sigs_04 = xs_04.iter().map(|x| x.sigs).sum::<f32>() / meshid.len() as f32;

//     let mu_02 = xs_02.iter().map(|x| x.mu).sum::<f32>() / meshid.len() as f32;
//     let mu_24 = xs_24.iter().map(|x| x.mu).sum::<f32>() / meshid.len() as f32;
//     let mu_04 = xs_04.iter().map(|x| x.mu).sum::<f32>() / meshid.len() as f32;

//     let nu_02 = xs_02.iter().map(|x| x.nu).sum::<f32>() / meshid.len() as f32;
//     let nu_24 = xs_24.iter().map(|x| x.nu).sum::<f32>() / meshid.len() as f32;
//     let nu_04 = xs_04.iter().map(|x| x.nu).sum::<f32>() / meshid.len() as f32;

//     let chi_02 = xs_02.iter().map(|x| x.chi).sum::<f32>() / meshid.len() as f32;
//     let chi_24 = xs_24.iter().map(|x| x.chi).sum::<f32>() / meshid.len() as f32;
//     let chi_04 = xs_04.iter().map(|x| x.chi).sum::<f32>() / meshid.len() as f32;

//     println!(
//         "2 group 1 xs\nsigt: {}\nsigf: {}\nsiga: {}\nsigs: {}\nmu: {}\nnu: {}\nchi: {}\n",
//         sigt_02, sigf_02, siga_02, sigs_02, mu_02, nu_02, chi_02
//     );
//     println!(
//         "2 group 2 xs\nsigt: {}\nsigf: {}\nsiga: {}\nsigs: {}\nmu: {}\nnu: {}\nchi: {}\n",
//         sigt_24, sigf_24, siga_24, sigs_24, mu_24, nu_24, chi_24
//     );
//     println!(
//         "1 group xs\nsigt: {}\nsigf: {}\nsiga: {}\nsigs: {}\nmu: {}\nnu: {}\nchi: {}\n",
//         sigt_04, sigf_04, siga_04, sigs_04, mu_04, nu_04, chi_04
//     );
// }

fn main() {
    // let now = SystemTime::now();

    // let (variables, xsdata, matid, deltax, solution, solver) = process_input();

    // let (meshid, fuel_indices) = mesh_gen(matid, &variables, &deltax);

    // let results = match (solution, solver) {
    //     (1, _) => monte_carlo(&variables, &xsdata, &deltax, &meshid, &fuel_indices, 1.0),
    //     (_, Solver::LinAlg) => nalgebra_method(
    //         &xsdata,
    //         &meshid,
    //         variables.energygroups,
    //         variables.mattypes,
    //         variables.boundl,
    //         variables.boundr,
    //         variables.numass,
    //     ),
    //     (_, _) => SolutionResults {
    //         flux: Vec::new(),
    //         assembly_average: Vec::new(),
    //         fission_source: Vec::new(),
    //         k: Vec::new(),
    //         k_fund: Vec::new(),
    //     }, // not implemented
    // };

    // // energy_collapse(&meshid, variables.mattypes, variables.energygroups, variables.numrods, variables.mpfr, results.flux.clone(), &xsdata);

    // let _ = plot_solution(
    //     results,
    //     variables.energygroups,
    //     variables.generations,
    //     meshid.len(),
    //     meshid[meshid.len() - 1].mesh_right,
    // );

    // print!(
    //     "Run was completed in {} milliseconds \n",
    //     now.elapsed().unwrap().as_millis()
    // );

    // below is for timing
    let mut now = SystemTime::now();

    for zyn in 0..100000 {
        if zyn % 10000 == 0 {
            print!(
                "Average time over those 10000 runs was {} seconds \n",
                now.elapsed().unwrap().as_micros() / 10000
            );
            now = SystemTime::now();
        }
        let _ = process_input();
    }
}
