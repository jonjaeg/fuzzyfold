use std::fs::File;
use std::io::{BufRead, BufReader, stdin};

use anyhow::{anyhow, Result};
use clap::{Parser, ValueEnum};

use fuzzyfold::energy::{ViennaRNA, NucleotideVec, PseudoEnergyModel, parse_structure};
use fuzzyfold::energy::parameters::{RNA_TURNER_2004, RNA_MT09, RNA_DP03, RNA_DP09};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum PkChoice {
    /// Turner 2004 NN + dp03 PK (original HotKnots default)
    Dp03,
    /// Turner 2004 NN + dp09 PK (incorrect pairing — use mt09 instead)
    Dp09,
    /// mt09 NN + dp09 PK (correct pairing; dp09 was trained with mt09)
    Mt09,
}

#[derive(Debug, Parser)]
#[command(name = "ff-calc-pseudo")]
#[command(about = "Evaluate pseudoknot free energy using the Dirks-Pierce model.\n\
                   Parameter sets: dp03 (Turner04+dp03), dp09 (Turner04+dp09), \
                   mt09 (mt09+dp09, recommended).")]
struct Cli {
    /// Input file (FASTA-like: optional >header, sequence, structure), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Pseudoknot parameter set
    #[arg(long, value_enum, default_value = "dp03")]
    pk_params: PkChoice,

    /// Temperature in Celsius (ignored for mt09, which is fixed at 37 °C)
    #[arg(long, default_value = "37.0")]
    celsius: f64,
}

fn read_input(path: &str) -> Result<(String, String)> {
    let reader: Box<dyn BufRead> = if path == "-" {
        Box::new(BufReader::new(stdin()))
    } else {
        Box::new(BufReader::new(File::open(path)?))
    };

    let mut sequence: Option<String> = None;
    let mut structure: Option<String> = None;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            continue;
        }
        if sequence.is_none() {
            sequence = Some(line.replace(' ', ""));
        } else {
            structure = Some(line.replace(' ', ""));
            break;
        }
    }

    let seq = sequence.ok_or_else(|| anyhow!("Missing sequence line"))?;
    let st  = structure.ok_or_else(|| anyhow!("Missing structure line"))?;
    Ok((seq, st))
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let (seq_str, st_str) = read_input(&cli.input)?;

    let seq   = NucleotideVec::try_from_rna(&seq_str)?;
    let loops = parse_structure(&st_str)?;

    let model = match cli.pk_params {
        PkChoice::Dp03 => ViennaRNA::from_thermo_params(&RNA_TURNER_2004, cli.celsius)
            .with_pseudoknot_params(RNA_DP03),
        PkChoice::Dp09 => ViennaRNA::from_thermo_params(&RNA_TURNER_2004, cli.celsius)
            .with_pseudoknot_params(RNA_DP09),
        PkChoice::Mt09 => ViennaRNA::from_andrunescu_params(&RNA_MT09)
            .with_pseudoknot_params(RNA_DP09),
    };

    let energy = model.energy_of_pseudoknotted_structure(&seq, &loops)?;
    println!("{:.4}", energy as f64 / 100.0);
    Ok(())
}
