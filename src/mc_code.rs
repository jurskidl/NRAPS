use rand::prelude::*;
use std::thread;

use crate::{DeltaX, Mesh, SolutionResults, Variables, XSData};

fn energy(
    chi: f32,
    index: usize,
    variables: &Variables,
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
) -> u8 {
    let mut temp = 0.0;
    let chit: Vec<f32> = xsdata
        .chit
        .iter()
        .skip(meshid[index].matid as usize)
        .step_by(variables.mattypes as usize)
        .map(|x| {
            temp += x;
            return temp;
        })
        .collect();
    chit.iter().position(|&x| x >= chi).unwrap() as u8
}

fn direction(mu: f32) -> f32 {
    2.0 * (mu) - 1.0
}

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

fn hit_boundary(mu: f32, start_x: f32, delta_s: f32, bound: f32, mesh_end: f32) -> (f32, f32, f32) {
    (
        mu * (-bound),
        (delta_s + (start_x - mesh_end)) * (-bound),
        mesh_end,
    )
}

fn cross_mesh(
    mesh_index: usize,
    mu: f32,
    start_x: f32,
    mesh_end: f32,
    delta_s: f32,
) -> (f32, f32, usize) {
    let mesh_index = if mu >= 0.0 {
        mesh_index + 1
    } else {
        mesh_index - 1
    };

    (delta_s + (start_x - mesh_end), mesh_end, mesh_index)
}

fn interaction(
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    mesh_index: usize,
    mattypes: u8,
    energygroups: u8,
    neutron_energy: u8,
) -> (bool, u8, f32) {
    let interaction: f32 = random();
    let absorption: f32 = xsdata.siga
        [(meshid[mesh_index].matid + (mattypes * neutron_energy)) as usize]
        * xsdata.sigt[(meshid[mesh_index].matid + (mattypes * neutron_energy)) as usize].powi(-1);

    let mut scat_mat = vec![0.0; energygroups as usize];

    let scat: f32 = xsdata.sigs[(meshid[mesh_index].matid + (mattypes * neutron_energy)) as usize];

    for energy in 0..energygroups {
        scat_mat[energy as usize] = (0..energy + 1)
            .into_iter()
            .map(|_energy| {
                xsdata.scat_matrix[((energygroups.pow(2) * meshid[mesh_index].matid)
                    + (energygroups * neutron_energy)
                    + _energy) as usize]
            })
            .sum::<f32>()
            / scat;
    }

    let random_scatter = random::<f32>();

    let scatter_energy = scat_mat.iter().position(|&x| x > random_scatter).unwrap() as u8;

    return match interaction {
        x if x < absorption => (false, neutron_energy, 0.0),
        _ => (true, scatter_energy, 2.0 * (random::<f32>()) - 1.0),
    };
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
    let mut delta_s: f32 = mu
        * -random::<f32>().ln()
        * xsdata.inv_sigtr[(meshid[mesh_index].matid + (mattypes * neutron_energy)) as usize];

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

            let (particle_exists, _neutron_energy, _mu) = interaction(
                &xsdata,
                meshid,
                mesh_index,
                mattypes,
                energygroups,
                neutron_energy,
            );
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
                * -random::<f32>().ln()
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
    return tally;
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
        // One threads for the OS to use
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
    println!("{}", results.k_fund[results.k_fund.len() - 1]);

    results
}

mod tests {
    extern crate float_cmp;
    #[allow(unused)]
    use crate::mc_code::{cross_mesh, direction, energy, hit_boundary, interaction};
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
}
