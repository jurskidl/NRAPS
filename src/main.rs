use std::iter::repeat;
// Use these for timing
use std::time::SystemTime;

use crate::discrete::{jacobi, nalgebra_method, succ_rel};
use crate::mc_code::monte_carlo;
use crate::process_input::process_input;
use crate::plot_solution::plot_solution;

mod discrete;
mod mc_code;
mod process_input;
mod plot_solution;

pub enum Solver {
    LinAlg,
    Gaussian,
    Jacobian,
    SR,
}

struct Variables {
    analk: u8,
    mattypes: u8,
    energygroups: u8,
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

fn main() {
    let now = SystemTime::now();

    let (variables, xsdata, matid, deltax, solution, solver) = process_input();

    let (meshid, fuel_indices) = mesh_gen(matid, &variables, &deltax);

    let results = match (solution, solver) {
        (1, _) => monte_carlo(&variables, &xsdata, &deltax, &meshid, &fuel_indices, 1.0),
        (_, Solver::LinAlg) => nalgebra_method(
            &xsdata,
            &meshid,
            variables.energygroups,
            variables.mattypes,
            variables.boundl,
            variables.boundr,
            variables.numass,
        ),
        (_, Solver::Jacobian) => jacobi(
            &xsdata,
            &meshid,
            variables.energygroups,
            variables.mattypes,
            variables.boundl,
            variables.boundr,
        ),
        (_, Solver::SR) => succ_rel(
            &xsdata,
            &meshid,
            variables.energygroups,
            variables.mattypes,
            variables.boundl,
            variables.boundr,
        ),
        (_, _) => SolutionResults {
            flux: Vec::new(),
            assembly_average: Vec::new(),
            fission_source: Vec::new(),
            k: Vec::new(),
            k_fund: Vec::new(),
        }, // not implemented
    };

    let _ = plot_solution(
        results,
        variables.energygroups,
        variables.generations,
        meshid.len(),
        meshid[meshid.len() - 1].mesh_right,
    );

    print!(
        "Run was completed in {} milliseconds \n",
        now.elapsed().unwrap().as_millis()
    );

    // below is for timing
    // let mut now = SystemTime::now();

    // for zyn in 0..100 {
    //     if zyn % 10 == 0 {
    //         print!(
    //             "Average time over those 10000 runs was {} seconds \n",
    //             now.elapsed().unwrap().as_secs() / 10
    //         );
    //         now = SystemTime::now();
    //     }
    //     let results = match solution {
    //         1 => monte_carlo(&variables, &xsdata, &deltax, &meshid, &fuel_indices, 1.0),
    //         _ => SolutionResults {
    //             flux: Vec::new(),
    //             fission_source: Vec::new(),
    //             k: Vec::new(),
    //         }, // not implemented
    //     };
    // }
}
