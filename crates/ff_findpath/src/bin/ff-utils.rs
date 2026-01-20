use ff_findpath::*;
use ff_structure::PairTable;
use std::error::Error;


fn main() -> Result<(), Box<dyn Error>> {
    // test example from max python
    let seq1 = "AGCCAUGAGUGUAUAGUGGGCCUAU";
    let struct1 = ".(((..............)))....";
    let struct2 = ".((................))....";
    let struct3 = ".(..................)....";
    let struct4 = ".........................";
    let struct5 = "..(...............)......";
    let struct6 = "..((.............))......";
    let structEnd = "..((((.........))))......";
    let pt1 = PairTable::try_from(struct1).unwrap();
    let pt2 = PairTable::try_from(struct2).unwrap();
    let pt3 = PairTable::try_from(struct3).unwrap();
    let pt4 = PairTable::try_from(struct4).unwrap();
    let pt5 = PairTable::try_from(struct5).unwrap();
    let pt6 = PairTable::try_from(struct6).unwrap();
    let ptEnd = PairTable::try_from(structEnd).unwrap();
       
    println!("Preparing moves from struct1 to struct2");
    
    let diff = compare_structures(&pt1, &ptEnd);
    println!("StructureDifference:\n{}", diff);
    println!("\nStart with the pt1 structure: {}", pt1);
    // analyze_folding_path(seq, start_struct, moves) -> (Vec<PathStep>, Stats)
    let (path_steps, stats) = analyze_folding_path(seq1, struct1, &diff.move_list);
    
    //println!("\nFolding path:");
    //for step in &path_steps {
    //    println!("{}", step);  // Dein schönes Format!
    //}
    
    println!("\nStats: {:?}", stats);
    Ok(())
}  
