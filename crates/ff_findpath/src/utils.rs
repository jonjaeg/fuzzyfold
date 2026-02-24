use std::fmt;
use ff_energy::{NucleotideVec, ViennaRNA, EnergyModel}; 
use ff_structure::{DotBracketVec, PairTable};

/// `Move` corresponds to an elementary move in the transformation from two RNA structures`.
///
/// it represent the insertion or deletion of a base pair at specific indices (i, j) between two structures. 
/// 
/// # Fields:
/// - `i`: The index of the first base in the pair.
/// - `j`: The index of the second base in the pair.
/// - `is_insertion`: A boolean flag indicating whether the move is an insertion (true) or deletion (false).
/// 
/// 
/// For example, if we have a move with i = 2, j = 5, is_insertion = true, it means that we are inserting the base pair (2, 5) into the structure.
/// The corresponding PairTable representation would have pt[2] = Some(5) and pt[5] = Some(2).
/// 
/// # Note:
/// We define the base pair distance between two structures as the number of moves required to transform one structure into the other.
/// This means one move corresponds to the change of two indices in the PairTable representation.
/// 
/// # Example usage:
/// ```rust
/// use ff_findpath::utils::Move;
/// 
/// let mv = Move { i: 2, j: 5, is_insertion: true };
/// ```
/// 
// Derive Default, makes default Move { i: 0, j: 0, is_insertion: false } (del (0,0) at start of every path)
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Move {
    pub i: usize,
    pub j: usize,
    pub is_insertion: bool, // true = insert, false = delete
}


/// StructureDifference captures the difference between two structures in terms of:
/// - move_list: A list of moves (insertions or deletions of base pairs) required to transform one structure into the other.
/// - bp_distance: The base pair distance between the two structures, defined as the number of base pairs that differ between them.
/// NOTE: Each move in the move_list corresponds to a change of two indices in the PairTable representation,
/// hence a bp_distance of one equals one move.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructureDifference {
    pub move_list: Vec<Move>,
    pub bp_distance: u32,
}



// Display implementations for better readability
impl fmt::Display for Move {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let kind = if self.is_insertion { "ins" } else { "del" };
        write!(f, "{}({}, {})", kind, self.i, self.j)
    }
}

impl fmt::Display for StructureDifference {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "bp_distance: {}", self.bp_distance)?;
        // don't print the number of moves because they are equal to bp_distance
        writeln!(f, "moves ({}):", self.move_list.len())?;
        for mv in &self.move_list {
            writeln!(f, "  - {}", mv)?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct Intermediate {
    pub pt: PairTable,
    pub saddle_energy: f64,
    pub current_energy: f64,
    pub remaining_moves: Vec<Move>, // Moves available to be taken
    pub path: Vec<Move>,            // Sequence of moves taken so far
}




/// compare_structures returns the ElemenrtaryMoves required to transform pt1 into pt2.
/// The data structure StructureDifference captures the list of moves, their corresponding hashes, and the base pair distance.
/// TODO: check again if it is right
/// Compares two PairTables and returns the moves required to transform pt1 into pt2.
pub fn compare_structures(pt1: &PairTable, pt2: &PairTable) -> StructureDifference {
    let mut diff = StructureDifference {
        move_list: vec![],
        bp_distance: 0,
    };

    let length = pt1.len();

    // Iterate through all bases
    for i in 0..length {
        // We only care if the pairing status at i is different
        if pt1[i] != pt2[i] {
            
            // 1. Check for Deletion from pt1
            // logic: if pt1 has a pair (i, j) starting at i
            if let Some(j) = pt1[i] {
                // Ensure we only process the pair once (at the opening bracket, where i < j)
                if i < j {
                    let deletion_move = Move {
                        i,
                        j,
                        is_insertion: false,
                    };

                    diff.move_list.push(deletion_move);
                    diff.bp_distance += 1;
                }
            }

            // 2. Check for Insertion into pt2
            // logic: if pt2 has a pair (i, j) starting at i
            if let Some(j) = pt2[i] {
                // Ensure we only process the pair once (at the opening bracket)
                if i < j {
                    let insertion_move = Move {
                        i,
                        j,
                        is_insertion: true,
                    };

                    
                    diff.move_list.push(insertion_move);
                    diff.bp_distance += 1;
                }
            }
        }
    }

    // Python does `move_list.sort(key=lambda x: x[0])`
    // We sort by the index `i` to maintain consistent order
    // TODO: is sorting useful, because if a insertion of a base pair (i, j) is before a deletion that affects the same indices,
    // then the order matters? ==> maybe other function solves this.
    
    //diff.move_list.sort_by_key(|m| m.i);
    // 2. Sort moves directly
    // Primary Key: Deletion (false) before Insertion (true)
    // Secondary Key: Index i (ascending)
    // Tertiary Key: Index j (ascending)
    diff.move_list.sort_by(|a, b| {
        a.is_insertion.cmp(&b.is_insertion)
            .then(a.i.cmp(&b.i))
            .then(a.j.cmp(&b.j))
    });
    diff
}

/// Prepares the list of moves required to transform pt_start into pt_target.
/// 
/// Translates the Python `_prepare_moves` function.
/// 
/// Note: The Python version maps raw tuples to `MoveState`. 
/// In Rust, `compare_structures` already produces `Move` structs, 
/// so we simply extract the `move_list`.
/// Wrapper for `compare_structures` function.
pub fn prepare_moves(pt_start: &PairTable, pt_target: &PairTable) -> Vec<Move> {
    // compare_structures returns StructureDifference { move_list, hash_list, bp_distance }
    // We only need the move_list for the greedy search.
    let diff = compare_structures(pt_start, pt_target);
    
    diff.move_list
}



/// Represents a single step in the folding trajectory for analysis.
#[derive(Debug, Clone)]
pub struct PathStep {
    pub structure: String,      // Dot-bracket representation
    pub move_applied: Option<Move>, // The move that led to this state (None for start)
    pub energy: f64,            // Energy of this structure
    pub step_index: usize,
}

impl fmt::Display for PathStep {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let move_str = match &self.move_applied {
            Some(mv) => format!("{}", mv),
            None => "Start".to_string(),
        };
        write!(
            f,
            "Step {}: Structure: {}, Move: {}, Energy: {:.2}",
            self.step_index, self.structure, move_str, self.energy
        )
    }
}

/// Summary statistics for the entire path.
#[derive(Debug, Clone)]
pub struct PathStats {
    pub saddle_energy: f64,     // Highest energy point (max_en)
    pub barrier_energy: f64,    // max_en - start_en
    pub start_energy: f64,
    pub end_energy: f64,
}


/// Evaluates a sequence of moves starting from s1.
/// Returns the detailed trajectory (steps) and summary stats.
/// NOTE: This function until now, only calculates energies using the ViennaRNA model!
/// BUG: the move sequence is not always valid! because no check is done if the move can be applied!
/// THIS leads to wrong secondary strucutre --> panics in energy calculation!
pub fn analyze_folding_path(
    sequence: &str,
    s1: &str,
    moves: &Vec<Move>,
) -> (Vec<PathStep>, PathStats) {
    // --- Optimization: Init Model & Sequence Once ---
    let model = ViennaRNA::default();
    let seq_vec = NucleotideVec::from_lossy(sequence);
    // ------------------------------------------------

    let mut current_pt = PairTable::try_from(s1).expect("Invalid start structure s1");
    
    // Use the pre-calculated seq_vec and model to calculate start energy
    let start_energy = model.energy_of_structure(&seq_vec, &current_pt) as f64 / 100.0;
    
    let mut current_max_energy = start_energy;
    let mut steps = Vec::with_capacity(moves.len() + 1);

    steps.push(PathStep {
        structure: s1.to_string(),
        move_applied: None,
        energy: start_energy,
        step_index: 0,
    });
    // print sequence
    println!("{}", sequence);
    // print starting structure and energy
    println!("{} \t\t {}", s1, start_energy);

    for (idx, m) in moves.iter().enumerate() {
        //println!("Applying move: {}", m);
        //println!("Current structure before move: {:?}", current_pt);
        // Apply move
        if m.is_insertion {
            current_pt[m.i] = Some(m.j);
            current_pt[m.j] = Some(m.i);
        } else {
            current_pt[m.i] = None;
            current_pt[m.j] = None;
        }
        //println!("Current PairTable after move: {}", current_pt);
        // print current DotBracket representation
        //println!("{}", DotBracketVec::try_from(&current_pt).unwrap());
        

        // Fast energy check using existing model/vec
        let energy = model.energy_of_structure(&seq_vec, &current_pt) as f64 / 100.0;
        //println!("Current energy: {}", energy);
        println!("{} \t {} \t {}", DotBracketVec::try_from(&current_pt).unwrap(), m, energy);

        if energy > current_max_energy {
            current_max_energy = energy;
        }

        let struct_str = format!("{:?}", current_pt); 

        steps.push(PathStep {
            structure: struct_str,
            move_applied: Some(m.clone()),
            energy,
            step_index: idx + 1,
        });
    }

    let end_energy = steps.last().map(|s| s.energy).unwrap_or(start_energy);

    // print final structure and energy
    //let last_dot_bracket_structure = DotBracketVec::try_from(&current_pt).unwrap();
    //println!("{} \t\t{}",&last_dot_bracket_structure , end_energy);

    let stats = PathStats {
        saddle_energy: current_max_energy,
        barrier_energy: current_max_energy - start_energy,
        start_energy,
        end_energy,
    };

    (steps, stats)
}

#[cfg(test)]
mod tests {
    // when importing ff_structure::PairTable into the utils.rs (outside of the tests module), 
    // This warning happens because PairTable is only used inside your tests module.
    // When you run cargo build (or when Rust Analyzer checks your code), 
    // the #[cfg(test)] module is ignored/stripped out. This leaves the import use ff_structure::PairTable; sitting at the top level with no one using it, hence the "unused import" warning.
    
    use ff_structure::PairTable;
    // use structs defined above (parent module), because tests is a submodule of the module where the structs are defined
    use super::*;

    #[test]
    fn test_valid_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(pt[0], Some(5));
        assert_eq!(pt[1], Some(4));
        assert_eq!(pt[2], None);
        assert_eq!(pt[3], None);
        assert_eq!(pt[4], Some(1));
        assert_eq!(pt[5], Some(0));
    }

    #[test]
    fn test_compare_structures_and_return_structure_difference() {
        // expected value
        let expected = StructureDifference {
            bp_distance: 7,
            move_list: vec![
                Move { i: 1, j: 20, is_insertion: false }, // del(1, 20)
                Move { i: 2, j: 19, is_insertion: false }, // del(2, 19)
                Move { i: 3, j: 18, is_insertion: false }, // del(3, 18)
                Move { i: 2, j: 18, is_insertion: true },  // ins(2, 18)
                Move { i: 3, j: 17, is_insertion: true },  // ins(3, 17)
                Move { i: 4, j: 16, is_insertion: true },  // ins(4, 16)
                Move { i: 5, j: 15, is_insertion: true },  // ins(5, 15)
            ]};
            

        // call the function under test (replace this with your real call)
        let pt1 = PairTable::try_from(".(((..............)))....").unwrap();
        let pt2 = PairTable::try_from("..((((.........))))......").unwrap();
        println!("pt1: {}", pt1);
        println!("pt2: {}", pt2);
        
        let got = compare_structures(&pt1, &pt2);
        println!("Prepared moves:");
        for m in &got.move_list  {
        println!("{:?}", m);
    }
        //println!("Got StructureDifference: {:?}", got);

        assert_eq!(got, expected);
    }

    #[test]
    fn test_vienna_energy_model(){
        use ff_energy::ViennaRNA;
        // energy_of_structure is exported via the EnergyModel trait
        use ff_energy::EnergyModel;
        use ff_energy::NucleotideVec;


        let model = ViennaRNA::default();

        let seq1 = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let struct1 = ".(((..............)))....";
        let expected_energy = -2.2;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq1), &PairTable::try_from(struct1).expect("valid")) as f64/100., expected_energy);
        
        let seq = "UCUACUAUUCCGGCUUGACAUAAAUAUCGAGUGCUCGACC";
        let dbr = "...........(.(((((........)))))..)......";
        let exp_energy = -210;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")), exp_energy);
        
    }
    #[test]
    fn test_analyze_folding_path(){
        use ff_structure::PairTable;
        use super::prepare_moves;

        let seq1 = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let struct1 = ".(((..............)))....";
        let struct2 = "..((((.........))))......";
        let pt1 = PairTable::try_from(struct1).unwrap();
        let pt2 = PairTable::try_from(struct2).unwrap();
        


        let moves = prepare_moves(&pt1, &pt2);
        let (steps, stats) = analyze_folding_path(seq1, struct1, &moves);
        //println!("Test steps:");
        //for step in steps {
        //    println!("{}", step);
        //} 
        println!("Stats: {:?}", stats);
    }
}
