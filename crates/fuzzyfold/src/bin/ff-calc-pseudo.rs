use std::fs::File;
use std::io::{BufRead, BufReader, stdin};

use anyhow::{Result, anyhow};
use clap::{Parser, ValueEnum};

use fuzzyfold::energy::parameters::{RNA_ANDRONESCU_2007, RNA_DP03, RNA_DP09, RNA_MT09, RNA_TURNER_2004};
use fuzzyfold::energy::{NucleotideVec, PseudoEnergyModel, ViennaRNA, parse_structure};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum PkChoice {
    /// MT09 stacking + DP09 PK params — correct pairing (DP09 was trained with MT09)
    Dp09,
    /// Turner 2004 stacking + DP03 PK params — original HotKnots default
    Dp03,
    /// Andronescu 2007 stacking + DP09 PK params — tests cross-model compatibility
    Andronescu07,
}

#[derive(Debug, Parser)]
#[command(name = "ff-calc-pseudo")]
#[command(
    about = "Evaluate pseudoknot free energy using the Dirks-Pierce model.\n\
                   Parameter sets:\n\
                   dp09  MT09 stacking + dp09 PK  (correct pairing, default)\n\
                   dp03  Turner04 stacking + dp03 PK (original HotKnots)"
)]
struct Cli {
    /// Input file (FASTA-like: optional >header, sequence, structure), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Pseudoknot parameter set
    #[arg(long, value_enum, default_value = "dp09")]
    pk_params: PkChoice,

    /// Temperature in Celsius (ignored for dp09, which uses MT09 fixed at 37 °C)
    #[arg(long, default_value = "37.0")]
    celsius: f64,

    /// Print the energy contribution of each loop individually
    #[arg(short = 'v', long)]
    verbose: bool,
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
    let st = structure.ok_or_else(|| anyhow!("Missing structure line"))?;
    Ok((seq, st))
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let (seq_str, st_str) = read_input(&cli.input)?;

    let seq = NucleotideVec::try_from_rna(&seq_str)?;
    let loops = parse_structure(&st_str)?;

    let model = match cli.pk_params {
        PkChoice::Dp09 => {
            ViennaRNA::from_andrunescu_params(&RNA_MT09).with_pseudoknot_params(RNA_DP09)
        }
        PkChoice::Dp03 => ViennaRNA::from_thermo_params(&RNA_TURNER_2004, cli.celsius)
            .with_pseudoknot_params(RNA_DP03),
        PkChoice::Andronescu07 => {
            ViennaRNA::from_andrunescu_params(&RNA_ANDRONESCU_2007).with_pseudoknot_params(RNA_DP09)
        }
    };

    if cli.verbose {
        let mut total = 0i32;
        for (i, lp) in loops.iter().enumerate() {
            let e = model.energy_of_pseudo_loop(&seq, lp)?;
            total += e;
            println!("{:3}  {:8.4}  {lp}", i, e as f64 / 100.0);
        }
        println!("{}", "-".repeat(60));
        println!("sum  {:8.4}", total as f64 / 100.0);
    } else {
        let energy = model.energy_of_pseudoknotted_structure(&seq, &loops)?;
        println!("{:.4}", energy as f64 / 100.0);
    }
    Ok(())
}
