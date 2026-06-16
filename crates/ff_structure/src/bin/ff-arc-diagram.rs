use ff_structure::arc_diagram::save_arc_diagram;
use ff_structure::PairTable;

use std::fs;
use std::path::Path;
use std::io;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Input file containing sequence and structure (two lines)
    #[arg(short, long, value_name = "INPUT_FILE")]
    input_file: String,

    /// Output SVG file path
    #[arg(short, long, value_name = "OUTPUT_FILENAME", default_value = "arc.svg")]
    output: String,
}

pub fn load_test_file(path: &str) -> Result<(String, String), io::Error> {
    let file_path = Path::new(path);

    if !file_path.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("Test file '{}' does not exist", path),
        ));
    }

    let content = fs::read_to_string(file_path)?;

    let lines: Vec<&str> = content
        .lines()
        .map(|line| line.trim())
        .filter(|line| !line.is_empty())
        .collect();

    if lines.len() < 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Expected at least two non-empty lines in {}", path),
        ));
    }

    Ok((lines[0].to_string(), lines[1].to_string()))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let (sequence, structure) = load_test_file(&args.input_file)?;

    let pt = PairTable::try_from(structure.as_str())?;

    save_arc_diagram(&sequence, &pt, &args.output)?;

    println!("Arc diagram written to '{}'", args.output);
    Ok(())
}