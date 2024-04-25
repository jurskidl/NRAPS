use crate::{Mesh, SolutionResults, XSData};
use nalgebra::*;

fn matrix_gen(
    n: usize,
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    neutron_energy: usize,
    mattypes: u8,
    boundl: f64,
    boundr: f64,
) -> Vec<Vec<f64>> {
    // Generate the matrix A
    let mut a = vec![vec![0.0; n]; n];

    // Set values for original insertion
    let d_curr: f64 = 3.0_f64.powi(-1)
        * xsdata.inv_sigtr[(meshid[0].matid + (mattypes * neutron_energy as u8)) as usize]
        * meshid[0].delta_x.powi(-1);
    let d_next: f64 = 3.0_f64.powi(-1)
        * xsdata.inv_sigtr[(meshid[1].matid + (mattypes * neutron_energy as u8)) as usize]
        * meshid[1].delta_x.powi(-1);
    let d_nextcurr: f64 = (2.0 * d_curr * d_next) * (d_curr + d_next).powi(-1);

    let beta_l: f64 = match boundl {
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
        * (xsdata.inv_sigtr[(meshid[0].matid + (mattypes * neutron_energy as u8)) as usize].powi(-1) 
        - xsdata.sigis[(meshid[0].matid + (mattypes * neutron_energy as u8)) as usize])
        + d_nextcurr;

    a[0][1] = -d_nextcurr;

    for x in 1..n - 1 {
        let d_curr: f64 = 3.0_f64.powi(-1)
            * xsdata.inv_sigtr[(meshid[x].matid + (mattypes * neutron_energy as u8)) as usize]
            * meshid[x].delta_x.powi(-1);
        let d_prev: f64 = 3.0_f64.powi(-1)
            * xsdata.inv_sigtr[(meshid[x - 1].matid + (mattypes * neutron_energy as u8)) as usize]
            * meshid[x - 1].delta_x.powi(-1);
        let d_next: f64 = 3.0_f64.powi(-1)
            * xsdata.inv_sigtr[(meshid[x - 1].matid + (mattypes * neutron_energy as u8)) as usize]
            * meshid[x + 1].delta_x.powi(-1);

        let d_prevcurr: f64 = (2.0 * d_curr * d_prev) * (d_curr + d_prev).powi(-1);
        let d_nextcurr: f64 = (2.0 * d_curr * d_next) * (d_curr + d_next).powi(-1);

        a[x][x - 1] = -d_prevcurr;
        a[x][x] = d_prevcurr
            + meshid[x].delta_x
            * (xsdata.inv_sigtr[(meshid[x].matid + (mattypes * neutron_energy as u8)) as usize].powi(-1) 
            - xsdata.sigis[(meshid[x].matid + (mattypes * neutron_energy as u8)) as usize])
            + d_nextcurr;
        a[x][x + 1] = -d_nextcurr;
    }

    // Set values for end insertion
    let d_curr: f64 = 3.0_f64.powi(-1)
        * xsdata.inv_sigtr[(meshid[n - 1].matid + (mattypes * neutron_energy as u8)) as usize]
        / meshid[n - 1].delta_x;
    let d_prev: f64 = 3.0_f64.powi(-1)
        * xsdata.inv_sigtr[(meshid[n - 2].matid + (mattypes * neutron_energy as u8)) as usize]
        / meshid[n - 2].delta_x;
    let d_prevcurr: f64 = (2.0 * d_curr * d_prev) / (d_curr + d_prev);

    // I feel like this shouldn't use d_next
    let beta_r = match boundr {
        x if x == 1.0 => 1.0,
        x if x == 0.0 => 0.25,
        _ => {
            (1.0 - 0.25 * ((1.0 - boundr) / (1.0 + boundr) * (1.0 / d_next)))
                / (1.0 + 0.25 * ((1.0 - boundr) / (1.0 + boundr) * (1.0 / d_curr)))
        }
    };

    a[n - 1][n - 2] = -d_prevcurr;
    a[n - 1][n - 1] = 2.0 * d_curr * (1.0 - beta_r)
        + meshid[n - 1].delta_x 
        * (xsdata.inv_sigtr[(meshid[n - 1].matid + (mattypes * neutron_energy as u8)) as usize].powi(-1) 
        - xsdata.sigis[(meshid[n - 1].matid + (mattypes * neutron_energy as u8)) as usize])
        + d_nextcurr;
    a
}

fn q_gen(xsdata: &XSData, energygroups: u8, mattypes: u8, flux: &Vec<Vec<f64>>, meshid: &Vec<Mesh>) -> Vec<Vec<f64>> {
    let mut q: Vec<Vec<f64>> = vec![vec![0.0; meshid.len()]; energygroups as usize];
    for neutron_energy in 0..energygroups {
        for index in 0..meshid.len(){
            q[neutron_energy as usize][index] = (0..energygroups).into_iter().map(|x| xsdata.nut[(meshid[index].matid + (mattypes * x)) as usize]
                *  xsdata.sigf[(meshid[index].matid + (mattypes * x)) as usize]
                * flux[x as usize][index]).sum::<f64>()
                * meshid[index].delta_x
                * xsdata.chit[(meshid[index].matid + (mattypes * neutron_energy)) as usize];
        }
    }
    q
}

pub fn nalgebra_method (
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    energygroups: u8,
    mattypes: u8,
    boundl: f64,
    boundr: f64,
) -> SolutionResults {
    let n: usize = meshid.len();
    let mut flux: Vec<Vec<f64>> = vec![vec![1.0; n]; energygroups as usize];
    let q: Vec<Vec<f64>> = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
    let (mut k, mut delta_flux, mut delta_k): (f64, f64, f64) = (1.0, 1.0, 1.0);

    for neutron_energy in 0..energygroups as usize {
        let temp_a = matrix_gen(n, xsdata, meshid, neutron_energy, mattypes, boundl, boundr);
        let a = DMatrix::from_vec(temp_a.len(), temp_a.len(), temp_a.into_iter().flatten().collect::<Vec<f64>>());
        let a_inv = a.try_inverse().unwrap();
        let temp_inv = a_inv.to_vector();
        let temp_vec: Vec<_> = temp_inv.chunks(n).collect();
        let temp2 = 2;
    }

    while delta_flux >= 1e-5 && delta_k >= 1e-6 {
        let temp_q = q.clone();

        for neutron_energy in 0..energygroups as usize {

            // let flux = a_inv[neutron_energy].dot(q[neutron_energy][0] * k.powi(-1));
            // flux_gen[gen+1][energy] = np.dot(A_inv[energy],Q_gen[gen][energy]/k_gen[gen] + SS_gen[gen][energy])

            // let scat: f64 = (0..energygroups).into_iter().map(|x| {
            //     match x {
            //         _neutron_energy =>  0.0,
            //         _ => xsdata.sigds[meshid[0].matid as usize + (mattypes as usize * 0)] * flux[x as usize][0] * meshid[0].delta_x,
            //     }
            // }).sum();
            // let temp_flux = k.powi(-1)
            //     * a[0][0].powi(-1)
            //     * (q[neutron_energy][0] - a[0][1] * flux[neutron_energy][1]);
            // delta_flux = ((flux[neutron_energy][0] - temp_flux)/flux[neutron_energy][0]).abs();
            // flux[neutron_energy][0] = temp_flux;

            // for index in 1..n - 1 {
            //     let temp = k.powi(-1)
            //         * a[index][index].powi(-1)
            //         * (q[neutron_energy][index] - a[index - 1][index] * flux[neutron_energy][index - 1]
            //         - a[index][index + 1] * flux[neutron_energy][index + 1]);
            //     delta_flux = ((flux[neutron_energy][index] - temp)/flux[neutron_energy][index]).abs().max(delta_flux);
            //     flux[neutron_energy][index] = temp;
            // }

            // let temp = k.powi(-1) * 
            //     a[n - 1][n - 1].powi(-1) * (q[neutron_energy][n-1] - a[n - 2][n - 1] * flux[neutron_energy][n - 2]);
            // delta_flux = ((flux[neutron_energy][n - 1] - temp)/flux[neutron_energy][n - 1]).abs().max(delta_flux);
            // flux[neutron_energy][n - 1] = temp;
        }
        

        let q: Vec<Vec<f64>> = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
        let temp_k = k;
        k = temp_k * (q.iter().flatten().sum::<f64>()/temp_q.iter().flatten().sum::<f64>());
        delta_k = ((k - temp_k)/temp_k).abs();

        println!("{}", k);
    }

    SolutionResults {
        flux,
        assembly_average: Vec::new(),
        fission_source: Vec::new(),
        k: vec![k],
        k_fund: Vec::new(),
    }
}

pub fn jacobi(
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    energygroups: u8,
    mattypes: u8,
    boundl: f64,
    boundr: f64,
) -> SolutionResults {
    let n: usize = meshid.len();
    let mut flux: Vec<Vec<f64>> = vec![vec![1000.0; n]; energygroups as usize];
    let q: Vec<Vec<f64>> = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
    let (mut k, mut delta_flux, mut delta_k): (f64, f64, f64) = (1.0, 1.0, 1.0);

    while delta_flux >= 1e-5 && delta_k >= 1e-6 {
        for neutron_energy in 0..energygroups as usize {
            let a = matrix_gen(n, xsdata, meshid, neutron_energy, mattypes, boundl, boundr);
        
            let mut temp_flux: Vec<f64> = vec![0.0; n];
            temp_flux[0] = (1.0/k) * (1.0 / a[0][0]) * (q[neutron_energy][0] - a[0][1] * temp_flux[1]);
            delta_flux = (temp_flux[0] - temp_flux[0]).abs();

            for index in 1..n - 1 {
                temp_flux[index] = (1.0/k)
                    * (1.0 / a[index][index])
                    * (q[neutron_energy][index] - a[index - 1][index] * flux[neutron_energy][index - 1]
                    - a[index][index + 1] * flux[neutron_energy][index + 1]);
                delta_flux = (flux[neutron_energy][index] - temp_flux[index]).abs().max(delta_flux);
            }

            temp_flux[n - 1] =( 1.0/k) * 
            (1.0 / a[n - 1][n - 1]) * (q[neutron_energy][n-1] - a[n - 2][n - 1] * flux[neutron_energy][n - 2]);
            delta_flux = (flux[neutron_energy][n - 1] - temp_flux[n - 1]).abs().max(delta_flux);
            flux[neutron_energy] = temp_flux;
        }
        let temp_q = q.clone();
        let q: Vec<Vec<f64>> = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
        let temp_k = k;
        k = temp_k * (q.iter().flatten().sum::<f64>()/temp_q.iter().flatten().sum::<f64>());
        delta_k = ((k - temp_k)/temp_k).abs();

        println!("{}", k);
    }
    SolutionResults {
        flux,
        assembly_average: Vec::new(),
        fission_source: Vec::new(),
        k: vec![k],
        k_fund: Vec::new(),
    }
}

pub fn succ_rel(
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    energygroups: u8,
    mattypes: u8,
    boundl: f64,
    boundr: f64,
) -> SolutionResults {
    let n: usize = meshid.len();
    let mut flux: Vec<Vec<f64>> = vec![vec![1.0; n]; energygroups as usize];
    let q: Vec<Vec<f64>> = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
    let (mut k, mut delta_flux, mut delta_k): (f64, f64, f64) = (1.0, 1.0, 1.0);
    
    while delta_flux >= 1e-5 && delta_k >= 1e-6 {
        let temp_q = q.clone();

        for neutron_energy in 0..energygroups as usize {
            let a = matrix_gen(n, xsdata, meshid, neutron_energy, mattypes, boundl, boundr);



            // let scat: f64 = (0..energygroups).into_iter().map(|x| {
            //     match x {
            //         _neutron_energy =>  0.0,
            //         _ => xsdata.sigds[meshid[0].matid as usize + (mattypes as usize * 0)] * flux[x as usize][0] * meshid[0].delta_x,
            //     }
            // }).sum();
            let temp_flux = k.powi(-1)
                * a[0][0].powi(-1)
                * (q[neutron_energy][0] - a[0][1] * flux[neutron_energy][1]);
            delta_flux = ((flux[neutron_energy][0] - temp_flux)/flux[neutron_energy][0]).abs();
            flux[neutron_energy][0] = temp_flux;

            for index in 1..n - 1 {
                let temp = k.powi(-1)
                    * a[index][index].powi(-1)
                    * (q[neutron_energy][index] - a[index - 1][index] * flux[neutron_energy][index - 1]
                    - a[index][index + 1] * flux[neutron_energy][index + 1]);
                delta_flux = ((flux[neutron_energy][index] - temp)/flux[neutron_energy][index]).abs().max(delta_flux);
                flux[neutron_energy][index] = temp;
            }

            let temp = k.powi(-1) * 
                a[n - 1][n - 1].powi(-1) * (q[neutron_energy][n-1] - a[n - 2][n - 1] * flux[neutron_energy][n - 2]);
            delta_flux = ((flux[neutron_energy][n - 1] - temp)/flux[neutron_energy][n - 1]).abs().max(delta_flux);
            flux[neutron_energy][n - 1] = temp;
        }
        

        let q: Vec<Vec<f64>> = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
        let temp_k = k;
        k = temp_k * (q.iter().flatten().sum::<f64>()/temp_q.iter().flatten().sum::<f64>());
        delta_k = ((k - temp_k)/temp_k).abs();

        println!("{}", k);
    }

    SolutionResults {
        flux,
        assembly_average: Vec::new(),
        fission_source: Vec::new(),
        k: vec![k],
        k_fund: Vec::new(),
    }
}
