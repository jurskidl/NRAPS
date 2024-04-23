use crate::{Mesh, SolutionResults, XSData};

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
        / meshid[0].delta_x;
    let d_next: f64 = 3.0_f64.powi(-1)
        * xsdata.inv_sigtr[(meshid[1].matid + (mattypes * neutron_energy as u8)) as usize]
        / meshid[1].delta_x;
    let d_nextcurr: f64 = (2.0 * d_curr * d_next) / (d_curr + d_next);

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
        + meshid[0].delta_x * xsdata.siga[(meshid[0].matid + (mattypes * neutron_energy as u8)) as usize]
        + d_nextcurr;

    for x in 1..n - 1 {
        let d_curr: f64 = 3.0_f64.powi(-1)
            * xsdata.inv_sigtr[(meshid[x].matid + (mattypes * neutron_energy as u8)) as usize]
            / meshid[x].delta_x;
        let d_prev: f64 = 3.0_f64.powi(-1)
            * xsdata.inv_sigtr[(meshid[x - 1].matid + (mattypes * neutron_energy as u8)) as usize]
            / meshid[x - 1].delta_x;
        let d_next: f64 = 3.0_f64.powi(-1)
            * xsdata.inv_sigtr[(meshid[x - 1].matid + (mattypes * neutron_energy as u8)) as usize]
            / meshid[x + 1].delta_x;

        let d_prevcurr: f64 = (2.0 * d_curr * d_prev) / (d_curr + d_prev);
        let d_nextcurr: f64 = (2.0 * d_curr * d_next) / (d_curr + d_next);

        a[x - 1][x] = -d_prevcurr;
        a[x][x] = d_prevcurr
            + d_nextcurr
            + meshid[x].delta_x
            * xsdata.siga[(meshid[x].matid + (mattypes * neutron_energy as u8)) as usize];
        a[x][x + 1] = -d_nextcurr;
    }

    // Set values for end insertion
    let d_curr: f64 = 3.0_f64.powi(-1)
        * xsdata.inv_sigtr[(meshid[n - 1].matid + (mattypes * neutron_energy as u8)) as usize]
        / meshid[n - 1].delta_x;
    let d_prev: f64 = 3.0_f64.powi(-2)
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

    a[n - 1][n - 1] = 2.0 * d_curr * (1.0 - beta_r)
        + meshid[n - 1].delta_x
        * xsdata.siga[(meshid[n - 1].matid + (mattypes * neutron_energy as u8)) as usize]
        + d_prevcurr;
    a
}

fn q_gen(xsdata: &XSData, energygroups: u8, mattypes: u8, flux: &Vec<Vec<f64>>, meshid: &Vec<Mesh>) -> Vec<Vec<f64>> {
    let mut q: Vec<Vec<f64>> = vec![vec![0.0; meshid.len()]; energygroups as usize];
    for neutron_energy in 0..energygroups as usize {
        for index in 0..meshid.len(){
            let nut = xsdata.nut[meshid[index].matid as usize + (mattypes as usize * neutron_energy)];
            let sigf = xsdata.sigf[meshid[index].matid as usize + (mattypes as usize * neutron_energy)];
            let deltax = meshid[index].delta_x;
            q[neutron_energy][index] = nut
                *  sigf
                * flux[neutron_energy][index] 
                * deltax;
        }
    }
    q
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
    let mut q: Vec<Vec<f64>> = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
    let (mut k, mut delta_flux, mut delta_k): (f64, f64, f64) = (1.0, 1.0, 1.0);

    while delta_flux >= 1e-5 && delta_k >= 1e-6 {
        let temp_q = q.clone();

        for neutron_energy in 0..energygroups as usize {
            let a = matrix_gen(n, xsdata, meshid, neutron_energy, mattypes, boundl, boundr);
        
            let mut temp_flux: Vec<f64> = vec![0.0; n];
            temp_flux[0] = (1.0/k) * (1.0 / a[0][0]) * (q[neutron_energy][0] - a[0][1] * temp_flux[1]);
            delta_flux = (temp_flux[0] - temp_flux[0]).abs();
            q[neutron_energy][0] = xsdata.nut[meshid[0].matid as usize + (mattypes as usize * neutron_energy)] 
                * xsdata.sigf[meshid[0].matid as usize + (mattypes as usize * neutron_energy)] 
                * temp_flux[0] 
                * meshid[0].delta_x;

            for index in 1..n - 1 {
                temp_flux[index] = (1.0/k)
                    * (1.0 / a[index][index])
                    * (q[neutron_energy][index] - a[index - 1][index] * flux[neutron_energy][index - 1]
                    - a[index][index + 1] * flux[neutron_energy][index + 1]);
                delta_flux = (flux[neutron_energy][index] - temp_flux[index]).abs().max(delta_flux);
                q[neutron_energy][index] = xsdata.nut[meshid[index].matid as usize + (mattypes as usize * neutron_energy)] 
                    * xsdata.sigf[meshid[index].matid as usize + (mattypes as usize * neutron_energy)] 
                    * temp_flux[index] 
                    * meshid[index].delta_x;
            }

            temp_flux[n - 1] =( 1.0/k) * 
            (1.0 / a[n - 1][n - 1]) * (q[neutron_energy][n-1] - a[n - 2][n - 1] * flux[neutron_energy][n - 2]);
            delta_flux = (flux[neutron_energy][n - 1] - temp_flux[n - 1]).abs().max(delta_flux);
            q[neutron_energy][n - 1] = xsdata.nut[meshid[n - 1].matid as usize + (mattypes as usize * neutron_energy)] 
                * xsdata.sigf[meshid[n - 1].matid as usize + (mattypes as usize * neutron_energy)] 
                * temp_flux[n - 1] 
                * meshid[n - 1].delta_x;
            flux[neutron_energy] = temp_flux;
        }
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
    let mut flux: Vec<Vec<f64>> = vec![vec![1000.0; n]; energygroups as usize];
    let mut q: Vec<Vec<f64>> = q_gen(xsdata, energygroups, mattypes, &flux, meshid);
    let (mut k, mut delta_flux, mut delta_k): (f64, f64, f64) = (1.0, 1.0, 1.0);
    
    while delta_flux >= 1e-5 && delta_k >= 1e-6 {
        let temp_q = q.clone();

        for neutron_energy in 0..energygroups as usize {
            let a = matrix_gen(n, xsdata, meshid, neutron_energy, mattypes, boundl, boundr);
            let temp_flux = (1.0/k) * (1.0 / a[0][0]) * (q[neutron_energy][0] - a[0][1] * flux[neutron_energy][1]);
            delta_flux = ((flux[neutron_energy][0] - temp_flux)/flux[neutron_energy][0]).abs();
            flux[neutron_energy][0] = temp_flux;
            q[neutron_energy][0] = xsdata.nut[meshid[0].matid as usize + (mattypes as usize * neutron_energy)] 
                * xsdata.sigf[meshid[0].matid as usize + (mattypes as usize * neutron_energy)] 
                * flux[neutron_energy][0] 
                * meshid[0].delta_x;

            for index in 1..n - 1 {
                let temp = (1.0/k)
                    * (1.0 / a[index][index])
                    * (q[neutron_energy][index] - a[index - 1][index] * flux[neutron_energy][index - 1]
                    - a[index][index + 1] * flux[neutron_energy][index + 1]);
                delta_flux = ((flux[neutron_energy][index] - temp)/flux[neutron_energy][index]).abs().max(delta_flux);
                flux[neutron_energy][index] = temp;
                q[neutron_energy][index] = xsdata.nut[meshid[index].matid as usize + (mattypes as usize * neutron_energy)] 
                    * xsdata.sigf[meshid[index].matid as usize + (mattypes as usize * neutron_energy)] 
                    * flux[neutron_energy][index] 
                    * meshid[index].delta_x;
            }

            let temp = (1.0/k) * 
                (1.0 / a[n - 1][n - 1]) * (q[neutron_energy][n-1] - a[n - 2][n - 1] * flux[neutron_energy][n - 2]);
            delta_flux = ((flux[neutron_energy][n - 1] - temp)/flux[neutron_energy][n - 1]).abs().max(delta_flux);
            flux[neutron_energy][n - 1] = temp;
            q[neutron_energy][n - 1] = xsdata.nut[meshid[n - 1].matid as usize + (mattypes as usize * neutron_energy)] 
                * xsdata.sigf[meshid[n - 1].matid as usize + (mattypes as usize * neutron_energy)] 
                * flux[neutron_energy][n - 1] 
                * meshid[n - 1].delta_x;
        }
        
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
