use ff_findpath::utils::{prepare_moves, analyze_folding_path};
use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA;

use std::fs;
use std::path::Path;
use std::io;

/// Reads sequence and structures from a three-line test file.
/// Returns a tuple (sequence, s1, s2) on success.
pub fn load_test_file(path: &str) -> Result<(String, String, String), io::Error> {
    let file_path = Path::new(path);

    // 1. Check if file exists (matches `if not file_path.is_file()`)
    if !file_path.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("Test file '{}' does not exist", path),
        ));
    }

    // 2. Read file content
    let content = fs::read_to_string(file_path)?;

    // 3. Process lines: strip whitespace and remove empty lines
    let lines: Vec<&str> = content
        .lines()
        .map(|line| line.trim())
        .filter(|line| !line.is_empty())
        .collect();

    // 4. Validate line count
    if lines.len() < 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Expected at least three non-empty lines in {}", path),
        ));
    }

    // 5. Return &str slices for sequence and structures
    Ok((
        lines[0].to_string(),
        lines[1].to_string(),
        lines[2].to_string(),
    ))
}


fn main()  -> Result<(), Box<dyn std::error::Error>> {
    // Example usage of load_test_file
    let test_file_path = "../../test_data/short.txt";

    // call the sequence, structures from the test file
    let (seq, struct1, struct2) = load_test_file(test_file_path)?; // These are Strings
    let struct1: &str = &struct1; // Convert String to &str
    let struct2: &str = &struct2;
    let pt1 = PairTable::try_from(struct1).unwrap(); 
    let pt2 = PairTable::try_from(struct2).unwrap();

    let moves = prepare_moves(&pt1, &pt2);
    let (steps, stats) = analyze_folding_path(&seq, &struct1, &moves);

    //println!("Test steps:");
    //    for step in steps {
    //        println!("{}", step);
    //    } 
    println!("Stats: {:?}", stats);
      
    Ok(())

}