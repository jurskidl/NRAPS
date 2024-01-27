use std::{
    fs::{self, File},
    io::{BufRead, BufReader},
};

pub enum Solver {
    LinAlg,
    Gaussian,
    Jacobian,
    Sor,
}

fn find_in_file(file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(file_path).expect("Unable to open the specified file");
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line.unwrap().to_ascii_lowercase();
        if line.starts_with("#") || line.is_empty() {
            continue;
        } else if line.contains("solution") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "solution").unwrap();
            let solution = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",solution);
        } else if line.contains("testcase") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "testcase").unwrap();
            let testcase = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",testcase);
        } else if line.contains("configs") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "configs").unwrap();
            let configs = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",configs);
        } else if line.contains("config") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "config").unwrap();
            let config = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",config);
        } else if line.contains("analk") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "analk").unwrap();
            let analk = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",analk);
        } else if line.contains("cases") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "cases").unwrap();
            let cases = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",cases);
        } else if line.contains("mattypes") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "mattypes").unwrap();
            let mattypes = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",mattypes);
        } else if line.contains("energygroups") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "energygroups").unwrap();
            let energygroups = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",energygroups);
        } else if line.contains("solver") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "solver").unwrap();
            let solver = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",solver);
        } else if line.contains("generations") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "generations").unwrap();
            let generations = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",generations);
        } else if line.contains("histories") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "histories").unwrap();
            let histories = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",histories);
        } else if line.contains("skip") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "skip").unwrap();
            let skip = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",skip);
        } else if line.contains("numass") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "numass").unwrap();
            let numass = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",numass);
        } else if line.contains("numrods") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "numrods").unwrap();
            let numrods = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",numrods);
        } else if line.contains("roddia") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "roddia").unwrap();
            let roddia = words[index+2].parse::<f64>().unwrap();
            print!("{}\n",roddia);
        } else if line.contains("rodpitch") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "rodpitch").unwrap();
            let rodpitch = words[index+2].parse::<f64>().unwrap();
            print!("{}\n",rodpitch);
        } else if line.contains("mpfr") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "mpfr").unwrap();
            let mpfr = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",mpfr);
        } else if line.contains("mpwr") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "mpwr").unwrap();
            let mpwr = words[index+2].parse::<i32>().unwrap();
            print!("{}\n",mpwr);
        } else if line.contains("boundl") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "boundl").unwrap();
            let boundl = words[index+2].parse::<f64>().unwrap(); // need to change this to save both values to a vector
            print!("{}\n",boundl);
        } else if line.contains("boundr") && line.contains("="){
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            let index = words.iter().position(|x| x == "boundr").unwrap();
            let boundr = words[index+2].parse::<f64>().unwrap(); // need to change this to save both values to a vector
            print!("{}\n",boundr);
        }
    }

    Ok(())
}

fn main() {
    let file_path = "../SampleInputFile.txt";

    find_in_file(file_path);

    println!("In file {}", file_path);

    //let contents = fs::read_to_string(file_path).expect("Should have been able to read the file");

    //println!("With text:\n{contents}");
}
