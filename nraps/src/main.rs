use std::{
    fs::{self, File},
    io::{BufRead, BufReader},
};



fn find_in_file(file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(file_path).expect("Unable to open the specified file");
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line.unwrap().to_ascii_lowercase();
        if line.starts_with("#") || line.is_empty() {
            continue;
        } else {
            let words = line
                .split_whitespace()
                .map(|s| s.to_owned())
                .collect::<Vec<String>>();

            match words[1]{
                &"solution" => let solution = words[3].parse()
            }
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
