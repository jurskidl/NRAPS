use rand::seq::index;

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

pub fn jacobi(
    xsdata: &XSData,
    meshid: &Vec<Mesh>,
    energygroups: u8,
    mattypes: u8,
    boundl: f64,
    boundr: f64,
) -> SolutionResults {
    let n: usize = meshid.len();
    let mut flux: Vec<Vec<f64>> = vec![vec![1.0; n]; energygroups as usize];
    for neutron_energy in 0..=energygroups as usize {
        let a = matrix_gen(n, xsdata, meshid, neutron_energy, mattypes, boundl, boundr);
        let q = 1.0;
        let mut delta = 1.0;

        while delta >= 1e-5 {
            let mut temp: Vec<f64> = vec![0.0; n];
            temp[0] = (1.0 / a[0][0]) * (q - a[0][1] * flux[neutron_energy][1]);
            delta = (flux[neutron_energy][0] - temp[0]).abs();

            for index in 1..n - 1 {
                temp[index] = (1.0 / a[index][index])
                    * (q - a[index - 1][index] * flux[neutron_energy][index - 1]
                    - a[index][index + 1] * flux[neutron_energy][index + 1]);
                delta = (flux[neutron_energy][index] - temp[index]).abs().max(delta);
            }

            temp[n - 1] =
                (1.0 / a[n - 1][n - 1]) * (q - a[n - 2][n - 1] * flux[neutron_energy][n - 2]);
            delta = (flux[neutron_energy][n - 1] - temp[n - 1]).abs().max(delta);
            flux[neutron_energy] = temp;
        }
    }
    SolutionResults {
        flux,
        assembly_average: Vec::new(),
        fission_source: Vec::new(),
        k: Vec::new(),
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
    let mut q: Vec<Vec<f64>> = vec![vec![1.0; n]; energygroups as usize];
    for neutron_energy in 0..=energygroups as usize {
        let a = matrix_gen(n, xsdata, meshid, neutron_energy, mattypes, boundl, boundr);
        let mut delta = 1.0;

        while delta >= 1e-5 {
            let temp = (1.0 / a[0][0]) * (q[neutron_energy][0] - a[0][1] * flux[neutron_energy][1]);
            delta = (flux[neutron_energy][0] - temp).abs();
            flux[neutron_energy][0] = temp;

            for index in 1..n - 1 {
                let temp = (1.0 / a[index][index])
                    * (q[neutron_energy][index] - a[index - 1][index] * flux[neutron_energy][index - 1]
                    - a[index][index + 1] * flux[neutron_energy][index + 1]);
                delta = (flux[neutron_energy][index] - temp).abs().max(delta);
                flux[neutron_energy][index] = temp;
            }

            let temp =
                (1.0 / a[n - 1][n - 1]) * (q[neutron_energy][n-1] - a[n - 2][n - 1] * flux[neutron_energy][n - 2]);
            delta = (flux[neutron_energy][n - 1] - temp).abs().max(delta);
            flux[neutron_energy][n - 1] = temp;
        }
    }

    SolutionResults {
        flux,
        assembly_average: Vec::new(),
        fission_source: Vec::new(),
        k: Vec::new(),
        k_fund: Vec::new(),
    }
}
