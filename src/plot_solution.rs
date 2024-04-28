use std::error::Error;
use std::process::Command;
use csv::Writer;

use crate::SolutionResults;

pub fn plot_solution(
    results: SolutionResults,
    energygroups: u8,
    generations: usize,
    number_meshes: usize,
    assembly_length: f64,
) -> Result<(), Box<dyn Error>> {
    let output_k_fund: Vec<String> = results.k_fund.iter().map(|x| x.to_string()).collect();
    let output_k: Vec<String> = results.k.iter().map(|x| x.to_string()).collect();
    let mut output_flux =
        vec![vec!["0.0".to_string(); results.flux[0].len()]; energygroups as usize];
    for energy in 0..energygroups as usize {
        for index in 0..results.fission_source.len() {
            output_flux[energy][index] = results.flux[energy][index].to_string();
        }
    }
    let mut output_average =
        vec![vec!["0.0".to_string(); results.assembly_average[0].len()]; energygroups as usize];
    for energy in 0..energygroups as usize {
        for index in 0..results.assembly_average.len() {
            output_average[energy][index] = results.assembly_average[energy][index].to_string();
        }
    }
    let output_fission: Vec<String> = results
        .fission_source
        .iter()
        .map(|x| x.to_string())
        .collect();

    let mut wtr_vars = Writer::from_path("./vars.csv")?;

    wtr_vars.write_record([&assembly_length.to_string()])?;
    wtr_vars.write_record([&number_meshes.to_string()])?;
    wtr_vars.write_record([&generations.to_string()])?;
    wtr_vars.flush()?;

    let mut wtr_k = Writer::from_path("./k_eff.csv")?;

    wtr_k.write_record(&output_k)?;
    wtr_k.write_record(&output_k_fund)?;
    wtr_k.flush()?;

    let mut wtr = Writer::from_path("./interface.csv")?;

    for energy in 0..energygroups as usize {
        wtr.write_record(&output_flux[energy])?;
    }
    for energy in 0..energygroups as usize {
        wtr.write_record(&output_average[energy])?;
    }
    wtr.write_record(&output_fission)?;
    wtr.flush()?;

    Command::new("python").arg("plot.py").spawn()?;

    Ok(())
}