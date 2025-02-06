use crate::{Mesh, SolutionResults, XSData};
use nalgebra::*;

fn matrix_gen(
    n: usize,
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    neutron_energy: usize,
    mattypes: u8,
    boundl: f32,
    boundr: f32,
    energygroups: u8,
) -> Vec<Vec<f32>> {
    // Generate the matrix A
    let mut a = vec![vec![0.0; n]; n];

    // Set values for original insertion
    let d_curr: f32 = 3.0_f32.powi(-1)
        * xsdata.inv_sigtr[(meshid[0].matid + (mattypes * neutron_energy as u8)) as usize]
        * meshid[0].delta_x.powi(-1);
    let d_next: f32 = 3.0_f32.powi(-1)
        * xsdata.inv_sigtr[(meshid[1].matid + (mattypes * neutron_energy as u8)) as usize]
        * meshid[1].delta_x.powi(-1);
    let d_nextcurr: f32 = (2.0 * d_curr * d_next) * (d_curr + d_next).powi(-1);

    let beta_l: f32 = match boundl {
        x if x == 1.0 => 1.0,
        x if x == 0.0 => 0.25,
        _ => {
            (1.0 - 0.25 * ((1.0 - boundl) / (1.0 + boundl) * (1.0 / d_next)))
                / (1.0 + 0.25 * ((1.0 - boundl) / (1.0 + boundl) * (1.0 / d_curr)))
        }
    };

    // Insert 0,0 and n,n since these differ from the pattern
    a[0][0] = 2.0 * d_curr * (1.0 - beta_l)
        + meshid[0].delta_x
            * (xsdata.sigt[(meshid[0].matid + (mattypes * neutron_energy as u8)) as usize]
                - xsdata.scat_matrix[(((energygroups + 1) * neutron_energy as u8)
                    + (energygroups.pow(2) * meshid[0].matid))
                    as usize])
        + d_nextcurr;

    a[0][1] = -d_nextcurr;

    for x in 1..n - 1 {
        let d_curr: f32 = 3.0_f32.powi(-1)
            * xsdata.inv_sigtr[(meshid[x].matid + (mattypes * neutron_energy as u8)) as usize]
            * meshid[x].delta_x.powi(-1);
        let d_prev: f32 = 3.0_f32.powi(-1)
            * xsdata.inv_sigtr[(meshid[x - 1].matid + (mattypes * neutron_energy as u8)) as usize]
            * meshid[x - 1].delta_x.powi(-1);
        let d_next: f32 = 3.0_f32.powi(-1)
            * xsdata.inv_sigtr[(meshid[x + 1].matid + (mattypes * neutron_energy as u8)) as usize]
            * meshid[x + 1].delta_x.powi(-1);

        let d_prevcurr: f32 = (2.0 * d_curr * d_prev) * (d_curr + d_prev).powi(-1);
        let d_nextcurr: f32 = (2.0 * d_curr * d_next) * (d_curr + d_next).powi(-1);

        a[x][x - 1] = -d_prevcurr;
        a[x][x] = d_prevcurr
            + meshid[x].delta_x
                * (xsdata.sigt[(meshid[x].matid + (mattypes * neutron_energy as u8)) as usize]
                    - xsdata.scat_matrix[(((energygroups + 1) * neutron_energy as u8)
                        + (energygroups.pow(2) * meshid[x].matid))
                        as usize])
            + d_nextcurr;
        a[x][x + 1] = -d_nextcurr;
    }

    // Set values for end insertion
    let d_curr: f32 = 3.0_f32.powi(-1)
        * xsdata.inv_sigtr[(meshid[n - 1].matid + (mattypes * neutron_energy as u8)) as usize]
        / meshid[n - 1].delta_x;
    let d_prev: f32 = 3.0_f32.powi(-1)
        * xsdata.inv_sigtr[(meshid[n - 2].matid + (mattypes * neutron_energy as u8)) as usize]
        / meshid[n - 2].delta_x;
    let d_prevcurr: f32 = (2.0 * d_curr * d_prev) / (d_curr + d_prev);

    // I feel like this shouldn't use d_next
    let beta_r = match boundr {
        x if x == 1.0 => 1.0,
        x if x == 0.0 => 0.25,
        _ => {
            (1.0 - 0.25 * ((1.0 - boundr) / (1.0 + boundr) * (1.0 / d_next)))
                / (1.0 + 0.25 * ((1.0 - boundr) / (1.0 + boundr) * (1.0 / d_curr)))
        }
    };

    // [(mattype * energygroups) + ((energygroups * starting_energy) + final_energy)]

    a[n - 1][n - 2] = -d_prevcurr;
    a[n - 1][n - 1] = 2.0 * d_curr * (1.0 - beta_r)
        + meshid[n - 1].delta_x
            * (xsdata.sigt[(meshid[n - 1].matid + (mattypes * neutron_energy as u8)) as usize]
                - xsdata.scat_matrix[(((energygroups + 1) * neutron_energy as u8)
                    + (energygroups.pow(2) * meshid[n - 1].matid))
                    as usize])
        + d_nextcurr;
    a
}

fn q_gen(
    xsdata: &XSData,
    energygroups: u8,
    mattypes: u8,
    flux: &Vec<Vec<f32>>,
    meshid: &Vec<Mesh>,
) -> Vec<Vec<f32>> {
    let mut q: Vec<Vec<f32>> = vec![vec![0.0; meshid.len()]; energygroups as usize];
    for neutron_energy in 0..energygroups {
        for index in 0..meshid.len() {
            q[neutron_energy as usize][index] = (0..energygroups)
                .into_iter()
                .map(|x| {
                    xsdata.nut[(meshid[index].matid + (mattypes * x)) as usize]
                        * xsdata.sigf[(meshid[index].matid + (mattypes * x)) as usize]
                        * flux[x as usize][index]
                })
                .sum::<f32>()
                * meshid[index].delta_x
                * xsdata.chit[(meshid[index].matid + (mattypes * neutron_energy)) as usize];
        }
    }
    q
}

fn scat_calc(
    index: usize,
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    flux: &Vec<Vec<f32>>,
    neutron_energy: usize,
    energygroups: u8,
) -> f32 {
    let mut scat = 0.0;
    for energy in 0..energygroups {
        if energy != neutron_energy as u8 {
            scat += xsdata.scat_matrix[((energygroups.pow(2) * meshid[index].matid)
                + (energygroups * energy)
                + neutron_energy as u8) as usize]
                * flux[energy as usize][index]
                * meshid[index].delta_x;
        } else {
            continue;
        }
    }
    scat
}

fn average_assembly(flux: Vec<Vec<f32>>, numass: u8, energygroups: u8) -> Vec<Vec<f32>> {
    let mesh_assembly = flux[0].len() / numass as usize;
    let mut average = vec![vec![0.0; flux[0].len()]; energygroups as usize];
    for energy in 0..energygroups as usize {
        for assembly in 1..=numass as usize {
            let average_ass: f32 = (((assembly - 1) * mesh_assembly)..(assembly * mesh_assembly))
                .into_iter()
                .map(|x| flux[energy][x])
                .sum::<f32>()
                / mesh_assembly as f32;
            for index in ((assembly - 1) * mesh_assembly)..(assembly * mesh_assembly) {
                average[energy][index] = average_ass;
            }
        }
    }
    average
}

pub fn nalgebra_method(
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    energygroups: u8,
    mattypes: u8,
    boundl: f32,
    boundr: f32,
    numass: u8,
) -> SolutionResults {
    let n: usize = meshid.len();
    let mut flux: Vec<Vec<f32>> = vec![vec![1.0; n]; energygroups as usize];
    let mut q: Vec<Vec<f32>> = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
    let (mut k, mut delta_flux, mut delta_k): (f32, f32, f32) = (1.0, 1.0, 1.0);

    let mut a_inv_matrix: Vec<Vec<Vec<f32>>> = Vec::with_capacity(energygroups as usize * n * n);

    for neutron_energy in 0..energygroups as usize {
        let temp_a = matrix_gen(
            n,
            xsdata,
            meshid,
            neutron_energy,
            mattypes,
            boundl,
            boundr,
            energygroups,
        );
        let a = DMatrix::from_vec(
            temp_a.len(),
            temp_a.len(),
            temp_a.into_iter().flatten().collect::<Vec<f32>>(),
        );
        let a_inv = a.try_inverse().unwrap();
        let temp_inv: Vec<f32> = a_inv.data.as_vec().to_owned();
        let temp_vec: Vec<Vec<f32>> = temp_inv.chunks(n).map(|x| x.to_owned()).collect();
        a_inv_matrix.push(temp_vec);
    }

    while delta_flux >= 1e-5 && delta_k >= 1e-6 {
        let temp_q = q.clone();

        for neutron_energy in 0..energygroups as usize {
            let mut scat = Vec::with_capacity(n);

            scat.push(scat_calc(
                0,
                xsdata,
                meshid,
                &flux,
                neutron_energy,
                energygroups,
            ));

            for index in 1..n - 1 {
                scat.push(scat_calc(
                    index,
                    xsdata,
                    meshid,
                    &flux,
                    neutron_energy,
                    energygroups,
                ));
            }
            scat.push(scat_calc(
                n - 1,
                xsdata,
                meshid,
                &flux,
                neutron_energy,
                energygroups,
            ));

            //Calculate Flux vector now
            let temp_flux: Vec<f32> = q[neutron_energy]
                .iter()
                .zip(scat.iter())
                .map(|(x, y)| (x * k.powi(-1)) + y)
                .collect();
            let temp: f32 = a_inv_matrix[neutron_energy][0]
                .iter()
                .zip(temp_flux.iter())
                .map(|(x, y)| x * y)
                .sum::<f32>();
            delta_flux = (((flux[neutron_energy][0] - temp) / flux[neutron_energy][0]).abs())
                .max(delta_flux);
            flux[neutron_energy][0] = temp;

            for index in 1..n - 1 {
                let temp_flux: Vec<f32> = q[neutron_energy]
                    .iter()
                    .zip(scat.iter())
                    .map(|(x, y)| (x * k.powi(-1)) + y)
                    .collect();
                let temp: f32 = a_inv_matrix[neutron_energy][index]
                    .iter()
                    .zip(temp_flux.iter())
                    .map(|(x, y)| x * y)
                    .sum::<f32>();
                delta_flux = ((flux[neutron_energy][index] - temp) / flux[neutron_energy][0]).abs();
                flux[neutron_energy][index] = temp;
            }

            let temp_flux: Vec<f32> = q[neutron_energy]
                .iter()
                .zip(scat.iter())
                .map(|(x, y)| (x * k.powi(-1)) + y)
                .collect();
            let temp: f32 = a_inv_matrix[neutron_energy][n - 1]
                .iter()
                .zip(temp_flux.iter())
                .map(|(x, y)| x * y)
                .sum::<f32>();
            delta_flux = (((flux[neutron_energy][n - 1] - temp) / flux[neutron_energy][0]).abs())
                .max(delta_flux);
            flux[neutron_energy][n - 1] = temp;
        }

        q = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
        let temp_k = k;
        k = temp_k * (q.iter().flatten().sum::<f32>() / temp_q.iter().flatten().sum::<f32>());
        delta_k = ((k - temp_k) / temp_k).abs();
    }

    let nut_matrix = (0..energygroups)
        .into_iter()
        .map(|energy| {
            (0..meshid.len())
                .into_iter()
                .map(|index| xsdata.nut[(meshid[index].matid + (mattypes * energy as u8)) as usize])
                .collect::<Vec<f32>>()
        })
        .collect::<Vec<Vec<f32>>>();
    let sigf_matrix = (0..energygroups)
        .into_iter()
        .map(|energy| {
            (0..meshid.len())
                .into_iter()
                .map(|index| {
                    xsdata.sigf[(meshid[index].matid + (mattypes * energy as u8)) as usize]
                })
                .collect::<Vec<f32>>()
        })
        .collect::<Vec<Vec<f32>>>();

    let mut temp = vec![vec![1.0; n]; energygroups as usize];
    for energy in 0..energygroups as usize {
        for index in 0..meshid.len() {
            temp[energy][index] =
                flux[energy][index] * nut_matrix[energy][index] * sigf_matrix[energy][index];
        }
    }

    let power_flux = temp.into_iter().flatten().collect::<Vec<f32>>();
    let power_flux = power_flux
        .chunks(energygroups as usize)
        .map(|x| x.iter().sum::<f32>())
        .collect::<Vec<f32>>();
    let power_constant = 3565e6
        / (1.6022e-13
            * 200.0
            * power_flux
                .iter()
                .zip(meshid)
                .map(|(x, y)| x * y.delta_x)
                .sum::<f32>());

    flux = (0..energygroups as usize)
        .into_iter()
        .map(|energy| {
            flux[energy]
                .iter()
                .map(|x| x * power_constant)
                .collect::<Vec<f32>>()
        })
        .collect::<Vec<Vec<f32>>>();

    let temp_flux = flux.clone();

    println!("{:.10}", k);

    SolutionResults {
        flux,
        assembly_average: average_assembly(temp_flux, numass, energygroups),
        fission_source: Vec::new(),
        k: vec![k],
        k_fund: Vec::new(),
    }
}
