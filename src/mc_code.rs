use rand::prelude::*;
use std::thread;

use crate::{DeltaX, Mesh, SolutionResults, Variables, XSData};

#[inline(always)]
fn energy(
    chi: f32,
    index: usize,
    variables: &Variables,
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
) -> u8 {
    // Set values for use in later step. Optimized away anyway
    let skip = meshid[index].matid as usize;
    let step = variables.mattypes as usize;

    // Use pre-allocated array for better vectorization
    let mut cumulative = 0.0;

    let chit: Vec<f32> = xsdata
        .chit
        .iter()
        .skip(skip)
        .step_by(step)
        .map(|value| {
            cumulative += value;
            cumulative
        })
        .collect();
    chit.partition_point(|&x| x < chi).min(chit.len() - 1) as u8
}

#[inline(always)]
fn direction(mu: f32) -> f32 {
    2.0 * (mu) - 1.0
}

// #[inline(always)]
fn spawn_neutron(
    fuel_indices: &Vec<usize>,
    variables: &Variables,
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
) -> (usize, f32, f32, u8) {
    let index = fuel_indices[thread_rng().gen_range(0..fuel_indices.len())];
    (
        index,
        random::<f32>(),
        direction(random::<f32>()),
        energy(random::<f32>(), index, &variables, &xsdata, &meshid),
    )
}

#[inline(always)]
fn hit_boundary(mu: f32, start_x: f32, delta_s: f32, bound: f32, mesh_end: f32) -> (f32, f32, f32) {
    (
        mu * (-bound),
        (delta_s + (start_x - mesh_end)) * (-bound),
        mesh_end,
    )
}

#[inline(always)]
fn cross_mesh(
    mut mesh_index: usize,
    mu: f32,
    start_x: f32,
    mesh_end: f32,
    delta_s: f32,
) -> (f32, f32, usize) {
    mesh_index = if mu >= 0.0 {
        mesh_index + 1
    } else {
        mesh_index - 1
    };

    (delta_s + (start_x - mesh_end), mesh_end, mesh_index)
}

#[inline(always)]
fn scat_mat_calc(
    energygroups: u8,
    matid: u8,
    neutron_energy: u8,
    inv_sigs: f32,
    scat_matrix: &Vec<f32>,
) -> Vec<f32> {
    let base_idx: usize =
        ((energygroups.pow(2) * matid) + (energygroups * neutron_energy)) as usize;

    // let mut scat_mat = vec![0.0; energygroups as usize];
    let mut cumulative: f32 = 0.0;

    let scat_mat = (0..energygroups as usize)
        .into_iter()
        .map(|_energy| {
            cumulative += scat_matrix[base_idx + _energy];
            cumulative * inv_sigs
        })
        .collect();
    // for energy in 0..energygroups {
    //     scat_mat[energy as usize] = (0..energy + 1)
    //         .into_iter()
    //         .map(|_energy| xsdata.scat_matrix[base_idx + _energy as usize])
    //         .sum::<f32>()
    //         * inv_sigs;
    // }

    scat_mat
}

#[inline(always)]
fn interaction(
    interaction: f32,
    scat_mat: Vec<f32>,
    xsdata: &XSData,
    xs_index: usize,
    neutron_energy: u8,
) -> (bool, u8, f32) {
    let absorption: f32 = xsdata.siga[xs_index] / xsdata.sigt[xs_index];
    let mu = 2.0 * (random::<f32>()) - 1.0;

    let scatter_energy = scat_mat
        .partition_point(|&x| x < random::<f32>())
        .min(scat_mat.len() - 1) as u8;

    match interaction {
        x if x < absorption => (false, neutron_energy, 0.0),
        _ => (true, scatter_energy, mu),
    }
}

fn particle_travel(
    mut tally: Vec<Vec<f32>>,
    meshid: &Vec<Mesh>,
    mut mesh_index: usize,
    mut neutron_energy: u8,
    mut mu: f32,
    mut start_x: f32,
    mattypes: u8,
    energygroups: u8,
    boundr: f32,
    boundl: f32,
    xsdata: &XSData,
) -> (bool, Vec<Vec<f32>>, usize, f32, u8, f32) {
    let xs_index = (meshid[mesh_index].matid + (mattypes * neutron_energy)) as usize;
    let mut delta_s: f32 = mu * -random::<f32>().ln() * xsdata.inv_sigtr[xs_index];

    let mut same_material = true;
    while same_material == true {
        let end_x = start_x + delta_s;
        let mesh_end = if mu >= 0.0 {
            meshid[mesh_index].mesh_right
        } else {
            meshid[mesh_index].mesh_left
        };

        if (mu < 0.0 && mesh_end > end_x && mesh_index == 0)
            || (mu >= 0.0 && end_x > mesh_end && mesh_index == meshid.len() - 1)
        {
            // hit the boundary
            tally[neutron_energy as usize][mesh_index] += ((start_x - mesh_end) / mu).abs();
            let bound = if mu >= 0.0 { boundr } else { boundl };

            if bound > 0.0 {
                (mu, delta_s, start_x) = hit_boundary(mu, start_x, delta_s, bound, mesh_end);
            } else {
                return (false, tally, mesh_index, mu, neutron_energy, start_x);
            }
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

            let scat_mat = scat_mat_calc(
                energygroups,
                meshid[mesh_index].matid,
                neutron_energy,
                1.0 / xsdata.sigs[xs_index],
                &xsdata.scat_matrix,
            );

            let (particle_exists, _neutron_energy, _mu) =
                interaction(random::<f32>(), scat_mat, &xsdata, xs_index, neutron_energy);
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
            delta_s = mu * -random::<f32>().ln() * xsdata.inv_sigtr[xs_index];
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
) -> Vec<Vec<f32>> {
    let mut tally: Vec<Vec<f32>> = vec![vec![0.0; meshid.len()]; variables.energygroups as usize];

    for _y in start..=end {
        // spawn_sub_mesh is the partial distance through the mesh
        let (mut mesh_index, spawn_sub_mesh, mut mu, mut neutron_energy) =
            spawn_neutron(fuel_indices, &variables, &xsdata, &meshid);
        let mut start_x: f32 = meshid[mesh_index].mesh_left + (spawn_sub_mesh * delta_x.fuel);

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
                variables.energygroups,
                variables.boundr,
                variables.boundl,
                xsdata,
            );
        }
    }
    tally
}

fn average_assembly(mut results: SolutionResults, numass: u8, energygroups: u8) -> SolutionResults {
    let mesh_assembly = results.flux[0].len() / numass as usize;
    for energy in 0..energygroups as usize {
        for assembly in 1..=numass as usize {
            let average_ass: f32 = (((assembly - 1) * mesh_assembly)..(assembly * mesh_assembly))
                .into_iter()
                .map(|x| results.flux[energy][x])
                .sum::<f32>()
                / mesh_assembly as f32;
            for index in ((assembly - 1) * mesh_assembly)..(assembly * mesh_assembly) {
                results.assembly_average[energy][index] = average_ass;
            }
        }
    }
    results
}

pub fn monte_carlo(
    variables: &Variables,
    xsdata: &XSData,
    delta_x: &DeltaX,
    meshid: &Vec<Mesh>,
    fuel_indices: &Vec<usize>,
    mut k_new: f32,
) -> SolutionResults {
    let mut results = SolutionResults {
        flux: vec![vec![0.0; meshid.len()]; variables.energygroups as usize],
        assembly_average: vec![vec![0.0; meshid.len()]; variables.energygroups as usize],
        fission_source: vec![0.0; meshid.len()],
        k: vec![0.0; variables.generations],
        k_fund: vec![0.0; variables.generations],
    };

    println!("running MC code");

    for x in 0..variables.generations {
        let mut tally: Vec<Vec<f32>> =
            vec![vec![0.0; meshid.len()]; variables.energygroups as usize];
        let k = k_new;
        k_new = 0.0;

        // For Multithreading
        // One thread for the OS to use
        let threads: usize = thread::available_parallelism().unwrap().get() - 1;
        let threaded_histories = variables.histories / threads;
        let starting_points: Vec<usize> = (0..threads).map(|x| x * threaded_histories).collect();
        let mut ending_points: Vec<usize> =
            Vec::from_iter(starting_points[1..threads].iter().cloned());
        ending_points.push(variables.histories);

        thread::scope(|scope| {
            let mut tallies = Vec::with_capacity(threads);
            for thread in 0..threads {
                let start = starting_points[thread];
                let end = ending_points[thread];
                let tallied: thread::ScopedJoinHandle<'_, Vec<Vec<f32>>> = scope.spawn(move || {
                    {
                        particle_lifetime(
                            xsdata,
                            meshid,
                            fuel_indices,
                            variables,
                            delta_x,
                            start,
                            end,
                        )
                    }
                });
                tallies.push(tallied);
            }

            // Aggregate the results
            for tallied in tallies {
                let thread_result = tallied.join().unwrap();
                for energy in 0..variables.energygroups as usize {
                    for index in 0..thread_result[energy].len() {
                        tally[energy][index] += thread_result[energy][index]
                    }
                }
            }
        });
        let fund: f32 = 1.0 / (variables.generations - (variables.skip - 1)) as f32;

        for energy in 0..variables.energygroups as usize {
            for index in 0..tally[energy].len() {
                let delta_x = meshid[index].delta_x;
                let matid = meshid[index].matid;
                let flux = tally[energy][index] / (k * variables.histories as f32 * delta_x);
                let fission_source = xsdata.nut
                    [(matid + (variables.mattypes * energy as u8)) as usize]
                    * xsdata.sigf[(matid + (variables.mattypes * energy as u8)) as usize]
                    * flux;
                k_new += k * delta_x * fission_source;
                if x >= variables.skip {
                    let conversion: f32 = (3565e6 * k * 36.2)
                        / (200e6
                            * 1.602176634e-19
                            * xsdata.nut[(0 + (variables.mattypes * 1)) as usize]
                            * meshid[meshid.len() - 1].mesh_right);
                    results.flux[energy][index] += flux * conversion * fund;
                    results.fission_source[index] += fission_source * fund;
                }
            }
        }
        results.k[x] = k_new;
    }

    results = average_assembly(results, variables.numass, variables.energygroups);

    results.k_fund[variables.skip] = results.k[variables.skip];

    for generations in (variables.skip + 1)..variables.generations {
        results.k_fund[generations] = (variables.skip..=generations)
            .into_iter()
            .map(|x| results.k[x])
            .sum::<f32>()
            / (generations - (variables.skip - 1)) as f32;
    }
    //println!("{}", results.k_fund[results.k_fund.len() - 1]);

    results
}

mod tests {
    extern crate float_cmp;
    #[allow(unused)]
    use crate::mc_code::{cross_mesh, direction, energy, hit_boundary, interaction, scat_mat_calc};
    #[allow(unused)]
    use crate::{DeltaX, Mesh, SolutionResults, Variables, XSData};
    use float_cmp::approx_eq;
    #[allow(unused)]
    use float_cmp::ApproxEq;

    #[test]
    fn test_boundary() {
        let (mu, delta_s, start_x) = hit_boundary(0.5, 0.5, 1.0, 1.0, 1.0);
        assert!(mu.approx_eq(-0.5, (1e-9, 2)));
        assert!(delta_s.approx_eq(-0.5, (1e-9, 2)));
        assert!(start_x.approx_eq(1.0, (1e-9, 2)));

        let (mu, delta_s, start_x) = hit_boundary(-0.5, 0.5, -1.0, 1.0, 0.0);
        assert!(mu.approx_eq(0.5, (1e-9, 2)));
        assert!(delta_s.approx_eq(0.5, (1e-9, 2)));
        assert!(start_x.approx_eq(0.0, (1e-9, 2)));

        let (mu, delta_s, start_x) = hit_boundary(0.5, 0.5, 0.5, 1.0, 1.0);
        assert!(mu.approx_eq(-0.5, (1e-9, 2)));
        assert!(delta_s.approx_eq(0.0, (1e-9, 2)));
        assert!(start_x.approx_eq(1.0, (1e-9, 2)));

        let (mu, delta_s, start_x) = hit_boundary(-0.5, 0.5, -0.5, 1.0, 0.0);
        assert!(mu.approx_eq(0.5, (1e-9, 2)));
        assert!(delta_s.approx_eq(0.0, (1e-9, 2)));
        assert!(start_x.approx_eq(0.0, (1e-9, 2)));
    }

    #[test]
    fn test_cross_mesh() {
        let (delta_s, start_x, mesh_index) = cross_mesh(1, 0.5, 0.5, 1.0, 1.0);
        assert!(delta_s.approx_eq(0.5, (1e-9, 2)));
        assert!(start_x.approx_eq(1.0, (1e-9, 2)));
        assert_eq!(mesh_index, 2);

        let (delta_s, start_x, mesh_index) = cross_mesh(1, -0.5, 1.0, 0.5, -1.0);
        assert!(delta_s.approx_eq(-0.5, (1e-9, 2)));
        assert!(start_x.approx_eq(0.5, (1e-9, 2)));
        assert_eq!(mesh_index, 0);
    }

    #[test]
    fn test_direction() {
        let particle_direction = direction(1.0);
        assert!(particle_direction.approx_eq(1.0, (1e-9, 2)));

        let particle_direction = direction(0.75);
        assert!(particle_direction.approx_eq(0.5, (1e-9, 2)));

        let particle_direction = direction(0.5);
        assert!(particle_direction.approx_eq(0.0, (1e-9, 2)));

        let particle_direction = direction(0.25);
        assert!(particle_direction.approx_eq(-0.5, (1e-9, 2)));

        let particle_direction = direction(0.0);
        assert!(particle_direction.approx_eq(-1.0, (1e-9, 2)));
    }

    // #[test]
    // fn test_energy() {
    //     let test_energy = energy(1.0, &vec![0.0, 0.0, 0.5, 1.0]);
    //     assert_eq!(test_energy, 3);

    //     let test_energy = energy(0.51, &vec![0.0, 0.0, 0.5, 1.0]);
    //     assert_eq!(test_energy, 3);

    //     let test_energy = energy(0.5, &vec![0.0, 0.0, 0.5, 1.0]);
    //     assert_eq!(test_energy, 2);

    //     let test_energy = energy(0.49, &vec![0.0, 0.0, 0.5, 1.0]);
    //     assert_eq!(test_energy, 2);

    //     let test_energy = energy(1E-9, &vec![0.0, 0.0, 0.5, 1.0]);
    //     assert_eq!(test_energy, 2)
    // }

    #[test]
    fn test_scat_mat_calc() {
        let variables = Variables {
            analk: 1,
            mattypes: 4,
            energygroups: 2,
            generations: 2,
            histories: 10000,
            skip: 1,
            numass: 2,
            numrods: 2,
            roddia: 1.0,
            rodpitch: 0.5,
            mpfr: 1,
            mpwr: 2,
            boundl: 1.0,
            boundr: 1.0,
        };

        let xsdata = XSData {
            sigt: vec![0.200, 0.200, 0.200, 0.1, 1.00, 1.20, 1.10, 1.1],
            sigs: vec![0.200, 0.200, 0.200, 0.0, 0.80, 0.80, 1.10, 0.1],
            mu: vec![0.000, 0.000, 0.000, 0.0, 0.00, 0.00, 0.00, 0.0],
            siga: vec![0.000, 0.000, 0.000, 0.1, 0.20, 0.40, 0.00, 1.0],
            sigf: vec![0.000, 0.000, 0.000, 0.0, 0.18, 0.30, 0.00, 0.0],
            nut: vec![0.000, 0.000, 0.000, 0.0, 1.40, 1.50, 0.00, 0.0],
            chit: vec![1.000, 1.000, 0.000, 0.0, 0.00, 0.00, 0.00, 0.0],
            // Index via [(mattype * energygroups) + ((energygroups * starting_energy) + final_energy)]
            scat_matrix: vec![
                0.185, 0.015, 0.000, 0.800, 0.185, 0.015, 0.000, 0.800, 0.170, 0.030, 0.000, 1.100,
                0.000, 0.000, 0.000, 0.100,
            ],
            inv_sigtr: vec![
                5.0,
                5.0,
                5.0,
                10.0,
                1.00,
                0.83333333,
                0.9090909090909,
                0.9090909090909,
            ],
        };

        let scat_mat_check: Vec<Vec<f32>> = vec![
            vec![0.925, 1.0],
            vec![0.925, 1.0],
            vec![0.85, 1.0],
            vec![f32::NAN, f32::NAN], // skip this since it can't be compared
            vec![0.0, 1.0],
            vec![0.0, 1.0],
            vec![0.0, 1.0],
            vec![0.0, 1.0],
        ];

        let neutron_energy = 0;

        for matid in 0..3 {
            let xs_index = (matid + (variables.mattypes * neutron_energy)) as usize;

            let scat_mat = scat_mat_calc(
                variables.energygroups,
                matid,
                neutron_energy,
                1.0 / xsdata.sigs[xs_index],
                &xsdata.scat_matrix,
            );

            println!("{:?}", scat_mat);
            assert_eq!(
                scat_mat,
                scat_mat_check[((neutron_energy * 4) + matid) as usize]
            );
        }

        let neutron_energy = 1;

        for matid in 0..4 {
            let xs_index = (matid + (variables.mattypes * neutron_energy)) as usize;

            let scat_mat = scat_mat_calc(
                variables.energygroups,
                matid,
                neutron_energy,
                1.0 / xsdata.sigs[xs_index],
                &xsdata.scat_matrix,
            );
            assert_eq!(
                scat_mat,
                scat_mat_check[((neutron_energy * 4) + matid) as usize]
            );
        }
    }
}
