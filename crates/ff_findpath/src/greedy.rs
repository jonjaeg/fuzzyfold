use crate::utils::{Move, Intermediate, PathStats, PathStep, compare_structures};
use ff_energy::{NucleotideVec, ViennaRNA, EnergyModel};
use ff_structure::{PairTable, LoopTable, LoopInfo, DotBracketVec};

/// A lightweight struct to hold a potential neighbor of a secondary structure.
pub struct NeighborCandidate {
    pub pt: PairTable,
    pub applied_move: Move,
    pub original_move_index: usize, // Needed to remove from remaining_moves later
}





/// Generates all valid neighbors of the given input structure.
/// 
/// Checks if the move is valid (topology) and applies it to generate a new PairTable, as well as the move that was valid.
/// 
/// NOTE: This function does NOT calculate energy or clone paths. It only checks if the move is valid (e.g. no pseudoknots) and applies it to generate a new PairTable. 
/// The energy calculation and path cloning are done in the apply_move function, which calls this function to get valid neighbors before calculating their energies and updating the search state.
/// 
/// # Input
/// - `current_pt`: The current structure represented as a PairTable. This is the structure from which we want to generate neighbors by applying the remaining moves.
/// - `available_moves`: A list of moves that can potentially be applied to the current structure to generate neighbors. 
/// 
/// # Output
/// It returns a vector of `NeighborCandidate`, where each candidate includes the new PairTable resulting from applying a valid move, the move that was applied, and the original index of that move in the `available_moves` list.
pub fn generate_valid_neighbors(
    current_pt: &PairTable,
    available_moves: &[Move],
) -> Vec<NeighborCandidate> {
    let mut candidates = Vec::new();

    // LoopTable is calculated once per expansion (O(N))
    let loop_table = match LoopTable::try_from(current_pt) {
        Ok(lt) => lt,
        Err(_) => return Vec::new(), // Invalid current structure, no neighbors
    };

    for (idx, candidate_move) in available_moves.iter().enumerate() {
        let i = candidate_move.i;
        let j = candidate_move.j;

        // --- VALIDITY CHECKS (Topology) ---
        if candidate_move.is_insertion {
            // Check for pseudoknots: i and j must be in the same loop context
            match (&loop_table[i], &loop_table[j]) {
                (LoopInfo::Unpaired { l: l1 }, LoopInfo::Unpaired { l: l2 }) => {
                    if l1 != l2 { continue; } // Different loops = crossing = invalid
                },
                _ => continue, // Already paired
            }
        } else {
            // Deletion: The pair must actually exist
            if current_pt[i] != Some(j) { continue; }
        }

        // --- APPLY MOVE ---
        let mut new_pt = current_pt.clone();
        if candidate_move.is_insertion {
            new_pt[i] = Some(j);
            new_pt[j] = Some(i);
        } else {
            new_pt[i] = None;
            new_pt[j] = None;
        }

        candidates.push(NeighborCandidate {
            pt: new_pt,
            applied_move: candidate_move.clone(),
            original_move_index: idx,
        });
    }

    candidates
}

/// Performs one expansion step from the given Intermediate state.
/// 
/// It generates all valid neighbors by calling the `generate_valid_neighbors` function, calculates their energies,
///
/// # Input
/// - `model`: The ViennaRNA energy model to evaluate energies.
/// - `seq_vec`: The RNA sequence as a vector of nucleotides (NucleotideVec).
/// - `intermediate`: The current state of the search, including the current structure (as a PairTable), the saddle energy so far, the current energy, the list of remaining moves, and the path of moves taken to reach this state.
/// - `max_energy`: An optional parameter to set a maximum energy threshold for the search. If provided, any neighbor with energy above this threshold will be filtered out and not returned as a valid neighbor. 
/// This can be used to prune the search space and focus on more promising trajectories, but it is not used in the greedy algorithm (as it explores all valid neighbors and picks the best one).
/// 
/// # Output
/// It returns a vector of `Intermediate` states, each representing a valid neighbor of the input `intermediate`.
/// 
/// # Internal calls:
/// - generate_valid_neighbors to get valid topologies
pub fn apply_move(
    model: &ViennaRNA,
    seq_vec: &NucleotideVec,
    intermediate: &Intermediate,
    max_energy: Option<f64>,
) -> Vec<Intermediate> {
    
    // 1. Get valid topologies
    let candidates = generate_valid_neighbors(&intermediate.pt, &intermediate.remaining_moves);
    let mut results = Vec::new();

    for cand in candidates {
        // 2. Calculate Energy (Expensive)
        let en = model.energy_of_structure(seq_vec, &cand.pt) as f64 / 100.0;
        

        // NOT USED IN GREEDY as we always explore all moves and pick the best
        // 3. Filter (Optimization: Drop before cloning paths if eneregy too high)
        if let Some(max) = max_energy {
            if en >= max { continue; }
        }

        // 4. Update State (Expensive Clones happen ONLY here)
        let new_saddle = f64::max(intermediate.saddle_energy, en);
        
        let mut next_remaining = intermediate.remaining_moves.clone();
        next_remaining.remove(cand.original_move_index);

        let mut next_path = intermediate.path.clone();
        next_path.push(cand.applied_move);

        results.push(Intermediate {
            pt: cand.pt,
            saddle_energy: new_saddle,
            current_energy: en,
            remaining_moves: next_remaining,
            path: next_path,
        });
    }

    results
}

/// Calculates a folding path and an energy barrier between two RNA secondary structures using the greedy heuristics (Morgan-Higgs algorithm).
/// # Input 
/// It takes a 
/// -  `model`: The ViennaRNA energy model to evaluate energies., 
/// -  `sequence`: the RNA sequence as a string,
/// -  `s1`: the starting structure in dot-bracket notation, 
/// -  `s2`: the target structure in dot-bracket notation,

/// 
/// # Output
/// It returns a `Result` which is either:
/// - `Ok((Vec<PathStep>, PathStats))` containing the folding path (as a vector of `PathStep`) and some statistics about the path (as `PathStats`), or 
/// - `Err(String)` containing an error message if the input parameters are invalid or if the search fails (e.g. due to energy constraints or topological issues).
/// 
/// # Internal calls:
/// - `compare_structures()`: This function compares the starting and target structures to generate a list of moves that can transform one structure into the other. This is done once at the beginning to prepare the move lists for both directions.
/// - `apply_move()`: This function generates neighboring structures by applying the allowed moves to the current structure, while also calculating their energies
/// NOTE: the `max_energy` parameter is not used in greedy search, as it is a greedy algorithm that explores all valid neighbors and picks the best one.
pub fn greedy_find_path(
    model: &ViennaRNA,
    sequence: &str,
    s1: &str,
    s2: &str,
) -> Result<(Vec<PathStep>, PathStats), String> {
    
    // 1. Setup & Validation
    let seq_vec = NucleotideVec::from_lossy(sequence);
    let pt_start = PairTable::try_from(s1).map_err(|_| "Invalid start structure s1")?;
    let pt_target = PairTable::try_from(s2).map_err(|_| "Invalid target structure s2")?;

    // 2. Prepare Moves
    let comparison = compare_structures(&pt_start, &pt_target);
    let move_list = comparison.move_list;
    let total_steps = move_list.len();

    // 3. Initialize Energies
    let start_energy = model.energy_of_structure(&seq_vec, &pt_start) as f64 / 100.0;
    
    // 4. Initialize Trajectory
    let mut trajectory = Vec::with_capacity(total_steps + 1);
    
    // Record step 0 (Start)
    trajectory.push(PathStep {
        structure: DotBracketVec::try_from(&pt_start).unwrap().to_string(), 
        move_applied: None,
        energy: start_energy,
        step_index: 0,
    });

    // 5. Initialize Search State
    let mut current = Intermediate {
        pt: pt_start,
        saddle_energy: start_energy,
        current_energy: start_energy,
        remaining_moves: move_list,
        path: Vec::new(),
    };

    // 6. Greedy Loop
    for step in 1..=total_steps {
        
        // A. Generate Neighbors
        let candidates = apply_move(model, &seq_vec, &current, None);

        if candidates.is_empty() {
            return Err(format!("Greedy search stuck at step {}: No valid topology neighbors.", step));
        }

        // B. Select Best Neighbor (Min saddle, then Min current)
        let best_candidate = candidates
            .into_iter()
            .min_by(|a, b| {
                a.saddle_energy.partial_cmp(&b.saddle_energy)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| {
                        a.current_energy.partial_cmp(&b.current_energy)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
            })
            .unwrap();

        // C. Record Step
        let last_move = best_candidate.path.last().cloned();
        
        trajectory.push(PathStep {
            structure: DotBracketVec::try_from(&best_candidate.pt).unwrap().to_string(),
            move_applied: last_move,
            energy: best_candidate.current_energy,
            step_index: step,
        });

        // D. Advance State
        current = best_candidate;
    }

    // 7. Finalize Stats
    let stats = PathStats {
        saddle_energy: current.saddle_energy,
        barrier_energy: current.saddle_energy - start_energy,
        start_energy: start_energy,
        end_energy: current.current_energy,
    };

    Ok((trajectory, stats))
}




// #####################################
//                  TESTS
// #####################################




#[cfg(test)]
mod tests {
    use super::*;
    

    #[test]
    fn test_generate_valid_neighbors_01() {
        let db1 = "..((..))";
        //let db2 = "((......))";

        let pt1 = PairTable::try_from(db1).expect("Invalid structure");
        //let pt2 = PairTable::try_from(db2).expect("Invalid structure");
        let moves = vec![
            Move { i: 0, j: 7, is_insertion: true }, // InValid
            Move { i: 2, j: 7, is_insertion: false }, // Valid
            Move { i: 0, j: 5, is_insertion: true }, // Invalid (crossing)
        ];
        let neighbors = generate_valid_neighbors(&pt1, &moves);
        assert_eq!(neighbors.len(), 1, "Only one valid move should be found");
        assert_eq!(neighbors[0].applied_move.i, 2);
        assert_eq!(neighbors[0].applied_move.j, 7);
        assert_eq!(neighbors[0].applied_move.is_insertion, false);
    }

    #[test]
fn test_generate_valid_neighbors_02() {
    // --- SCENARIO 1: Simple Valid Moves ---
    // Start: "...." (Length 4, all unpaired)
    let pt_empty = PairTable::try_from("....").unwrap();
    
    // Potential moves:
    // 1. Pair (0, 3) -> "(..)"
    // 2. Pair (1, 2) -> ".(.)"
    // Initialize manually
    let moves = vec![
        Move { i: 0, j: 3, is_insertion: true }, 
        Move { i: 1, j: 2, is_insertion: true }, 
    ];

    let neighbors = generate_valid_neighbors(&pt_empty, &moves);

    println!("--- Scenario 1: Empty Structure ---");
    assert_eq!(neighbors.len(), 2, "Both moves should be valid on empty structure");
    // Verify the structures in the neighbors
    // Neighbor 0 should have 0-3 paired
    assert_eq!(neighbors[0].pt[0], Some(3));

}

    #[test]
    fn test_generate_valid_neighbors_03_pseudoknot(){
        // --- SCENARIO 2: Crossing / Pseudoknot Rejection ---
    // Start: "(...)." (Length 6)
    // Pair exists at 0-4. 
    // We try to pair index 3 (inside the loop) with index 5 (outside the loop).
    let pt_crossing = PairTable::try_from("(...).").unwrap();
    
    // This move would create a crossing link (pseudoknot)
    let bad_moves = vec![Move { i: 3, j: 5, is_insertion: true }];
    
    let neighbors_cross = generate_valid_neighbors(&pt_crossing, &bad_moves);
    
    println!("\n--- Scenario 2: Crossing Check ---");
    println!("Trying to pair (3,5) on structure '(...).'");
    println!("No valid neighbors, as expected. Would create pseudoknot.");
    assert_eq!(neighbors_cross.len(), 0, "Crossing move should be rejected automatically");

    }

    #[test]
    fn test_generate_valid_neighbors_printing_debugging(){
        // Does not test anything, just for visual inspection and debugging of the generated neighbors.
        let db1 = ".(((..............)))....";
        let db2 = "..((((.........))))......";

        let pt1 = PairTable::try_from(db1).expect("Invalid structure");
        let pt2 = PairTable::try_from(db2).expect("Invalid structure");

        let comparison = compare_structures(&pt1, &pt2);
        let move_list = comparison.move_list;

        println!("bp_distance ({}) to go from Structure 1 to Structure 2", move_list.len());
        println!("-----------------------------------------------------");

        let neighbors = generate_valid_neighbors(&pt1, &move_list);
        println!("Generated only {} valid neighbors from Structure 1 towards Structure 2", neighbors.len());
        for (idx, n) in neighbors.iter().enumerate() {
            println!("\nNeighbor {}: {:?} --> New Structure: {}", idx+1, n.applied_move, DotBracketVec::try_from(&n.pt).unwrap());
        }       



    }

    #[test]
    fn test_apply_move() {
        let seq1  = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let str1 = ".(((..............)))....";
        let str2 = "..((((.........))))......";
        let model = ViennaRNA::default();


        // 1. Setup & Validation
        let seq_vec = NucleotideVec::from_lossy(seq1);
        let pt_start = PairTable::try_from(str1).unwrap();
        let pt_target = PairTable::try_from(str2).unwrap();

        // 2. Prepare Moves
        let comparison = compare_structures(&pt_start, &pt_target);
        let move_list = comparison.move_list;
        

        // 3. Initialize Energies
        let start_energy = model.energy_of_structure(&seq_vec, &pt_start) as f64 / 100.0;


        // Initialize Intermediate
        let current = Intermediate {
            pt: pt_start,
            saddle_energy: start_energy,
            current_energy: start_energy,
            remaining_moves: move_list,
            path: Vec::new(),
        };
        let expanded = apply_move(&model, &seq_vec, &current, None);
        println!("Expanded into {} new structures.", expanded.len());
        println!("-------------------------------------");
        //println!("Expanded States:");
        for (idx, state) in expanded.iter().enumerate() {
            let db = DotBracketVec::try_from(&state.pt).unwrap();
            println!("State {}: Structure: {}\t Current Energy: {:.2} kcal/mol", idx+1, db, state.current_energy);
        }

        assert!(expanded.len() > 0, "At least one expansion should be possible");


    }




    #[test]
    fn test_greedy_find_path_integration_01() {
        // Simple hairpin closing example
        // Assume ViennaRNA default model
        
        let seq = "GGGGUUUUCCCC";
        let s1 = "............"; //(Unpaired)
        let s2 = "((((....))))"; //(Full Hairpin)
        
        // Initialize Model (This might fail if params aren't found in your env)
        // Ensure you have valid params or a mock model available.
        let model = ViennaRNA::default();
        let result = greedy_find_path(&model, seq, s1, s2);

        match result {
            Ok((trajectory, stats)) => {
                println!("Path found with {} steps", trajectory.len());
                println!("Saddle Energy: {:.2} kcal/mol", stats.saddle_energy);
                println!("Barrier Height: {:.2} kcal/mol", stats.barrier_energy);
                
                // Assertions
                assert_eq!(trajectory.len(), 5, "Should have 5 PathStep entries in the trajectory (including start and end)");
                // Check if the final path actually formed the structure
                // We know it starts empty, so all moves should be insertions
                
                for m in trajectory{
                    println!("{:?}", m)
                }
                }
            
            Err(e) => panic!("Greedy path search failed: {}", e),
        }
    }

    #[test]
    fn test_greedy_find_path_integration_02() {

       // Assume ViennaRNA default model
        
        let seq = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let s1 = ".(((..............)))....";
        let s2 = "..((((.........))))......";
        
        // Initialize Model (This might fail if params aren't found in your env)
        // Ensure you have valid params or a mock model available.
        let model = ViennaRNA::default();
        let result = greedy_find_path(&model, seq, s1, s2);

        match result {
            Ok((trajectory, stats)) => {
                println!("Path found with {} steps", trajectory.len());
                println!("Saddle Energy: {:.2} kcal/mol", stats.saddle_energy);
                println!("Barrier Height: {:.2} kcal/mol", stats.barrier_energy);
                
                // Assertions
                assert_eq!(trajectory.len(), 8, "Should have 8 PathStep entries in the trajectrory (including start and end)");
                // Check if the final path actually formed the structure
                for m in trajectory{
                    println!("{:?}", m)
                }
                
            },
            Err(e) => panic!("Greedy path search failed: {}", e),
        }
    }

    #[test]
    fn test_greedy_find_path_printing_debugging(){
        // Does not test anything, just prints the path and stats for visual inspection and debugging.
        let seq = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let s1 = ".(((..............)))....";
        let s2 = "..((((.........))))......";
        
        // Initialize Model (This might fail if params aren't found in your env)
        // Ensure you have valid params or a mock model available.
        let model = ViennaRNA::default();

        let (steps, stats) = greedy_find_path(&model, seq, s1, s2).unwrap();
        println!("Test steps:");
        for step in steps {
            println!("{} \t {} \t {}", step.structure,step.move_applied.unwrap_or_default(), step.energy );
        } 
        println!("Stats: {}", stats);
    }
}


