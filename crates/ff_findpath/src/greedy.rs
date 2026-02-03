use crate::utils::{Move, Intermediate, PathStats, PathStep, compare_structures};
use ff_energy::{NucleotideVec, ViennaRNA, EnergyModel};
use ff_structure::{PairTable, LoopTable, LoopInfo, DotBracketVec};

/// A lightweight struct to hold a potential neighbor before we commit to it.
pub struct NeighborCandidate {
    pub pt: PairTable,
    pub applied_move: Move,
    pub original_move_index: usize, // Needed to remove from remaining_moves later
}

// #[derive(Clone, Debug)]
// pub struct Intermediate {
// pub pt: PairTable,
// pub saddle_energy: f64,
// pub current_energy: f64,
// pub remaining_moves: Vec<Move>,
// pub path: Vec<Move>,
// }





/// generate_valid_neighbors checks which moves are valid from the current PairTable,
/// applies them, and returns a list of valid neighbor candidates without calculating energy.
/// Checks topology and generates valid PairTables. 
/// Does NOT calculate energy or clone paths yet.
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

/// greedy expand_state performs one expansion step from the given Intermediate state.
/// that means it generates all valid neighbors, calculates their energies,
/// filters them based on max_energy, and constructs new Intermediate states.
/// Calls the generate_valid_neighbors, calculates energy, filters, and builds Intermediate.
/// max_energy is an optional threshold to prune high-energy states early. (idea for findpath. with greedy this is not used!)
pub fn expand_state(
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


/// Finds a folding path from s1 to s2 using a greedy strategy.
///
/// At each step, it generates all valid moves that lead toward the target,
/// evaluates their energy, and picks the one that minimizes the saddle energy 
/// (and breaks ties with current energy).
pub fn _greedy_find_path(
    model: &ViennaRNA,
    sequence: &str, // Passed as &str for convenience, converted inside if needed
    s1: &str,
    s2: &str,
) -> Result<(Vec<Move>, f64, f64), String> {
    
    // 1. Setup Data Structures
    let seq_vec = &NucleotideVec::from_lossy(sequence); // Assuming NucleotideVec has From<&str>
    let pt_start = PairTable::try_from(s1).map_err(|_| "Invalid start structure s1")?;
    let pt_target = PairTable::try_from(s2).map_err(|_| "Invalid target structure s2")?;

    // 2. Prepare Moves (Compare start vs target)
    let comparison = compare_structures(&pt_start, &pt_target);
    let move_list = comparison.move_list;

    // Edge Case: If structures are identical, return empty path with 0 barrier
    if move_list.is_empty() {
        // Calculate energy of the static structure
        let start_energy = model.energy_of_structure(seq_vec, &pt_start) as f64 / 100.0;
        return Ok((Vec::new(), start_energy, 0.0));
    }

    // 3. Initialize Starting State
    let start_energy = model.energy_of_structure(&seq_vec, &pt_start) as f64 / 100.0;
    
    let mut current = Intermediate {
        pt: pt_start.clone(),
        saddle_energy: start_energy,
        current_energy: start_energy,
        remaining_moves: move_list.clone(),
        path: Vec::new(),
    };

    let total_steps = move_list.len();

    // 4. Greedy Loop
    // We iterate exactly 'total_steps' times because we must apply every move
    // in the list to reach the target structure.
    for _step in 1..=total_steps {
        
        // Generate and evaluate neighbors
        // We pass 'None' for max_energy as the Python script passes None in the loop
        let candidates = expand_state(model, &seq_vec, &current, None);

        if candidates.is_empty() {
            return Err("Greedy search got stuck: No valid moves available to progress toward target.".to_string());
        }

        // Greedy Selection: Minimize (saddle_energy, current_energy)
        // unwrap is safe because we checked !candidates.is_empty()
        let best_candidate = candidates
            .into_iter()
            .min_by(|a, b| {
                a.saddle_energy
                    .partial_cmp(&b.saddle_energy)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| {
                        a.current_energy
                            .partial_cmp(&b.current_energy)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
            })
            .unwrap();

        // Advance to next state
        current = best_candidate;
    }

    // 5. Finalize Results
    // The barrier energy is (Highest Point - Start Energy)
    let barrier_energy = current.saddle_energy - start_energy;

    Ok((current.path, current.saddle_energy, barrier_energy))
}


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
        let candidates = expand_state(model, &seq_vec, &current, None);

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
fn test_generate_simple_neighbors_02() {
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
    assert_eq!(neighbors_cross.len(), 0, "Crossing move should be rejected automatically");
}

    #[test]
    fn test_valid_move_generation() {
        let seq1  = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let str1 = ".(((..............)))....";
        let str2 = "..((((.........))))......";

        // Corrected Slicing for Balance (Length 12)
        // 1..6  => "GCCAU" (Indices 1,2,3,4,5). str1: "(((.." (3 open)
        // 13..20 => "AGUGGGC" (Indices 13..19). str1: ")))...." (3 closed)
        // Result: "(((..)))...." -> Balanced.
        let _subseq = format!("{}{}", &seq1[1..6], &seq1[13..20]);
        let substr1 = format!("{}{}", &str1[1..6], &str1[13..20]);
        let substr2 = format!("{}{}", &str2[1..6], &str2[13..20]); // Needs careful checking if balanced

        // Setup structures
        let pt1 = PairTable::try_from(substr1.as_str()).expect("substr1 invalid");
        let pt2 = PairTable::try_from(substr2.as_str()).expect("substr2 invalid");

        // Calculate differences (Moves to go from pt1 -> pt2)
        let diff = compare_structures(&pt1, &pt2);
        
        // --- TEST THE MOVE GENERATOR ---
        // We do NOT need 'model' or 'seq_vec' here!
        let candidates = generate_valid_neighbors(&pt1, &diff.move_list);

        println!("Found {} valid neighbors", candidates.len());
        
        for c in candidates {
            println!("Move: {:?}, New Structure Valid? {}", c.applied_move, true);
            // Optional: assert that c.pt is what you expect
        }
    }





    #[test]
    fn test_greedy_find_path_integration_01() {
        // Simple hairpin closing example
        // Seq: GGGGUUUUCCCC
        // S1:  ............ (Unpaired)
        // S2:  ((((....)))) (Full Hairpin)
        
       // Assume ViennaRNA default model
        
        let seq = "GGGGUUUUCCCC";
        let s1 = "............";
        let s2 = "((((....))))";
        
        // Initialize Model (This might fail if params aren't found in your env)
        // Ensure you have valid params or a mock model available.
        let model = ViennaRNA::default();
        let result = _greedy_find_path(&model, seq, s1, s2);

        match result {
            Ok((path, saddle, barrier)) => {
                println!("Path found with {} steps", path.len());
                println!("Saddle Energy: {:.2} kcal/mol", saddle);
                println!("Barrier Height: {:.2} kcal/mol", barrier);
                
                // Assertions
                assert_eq!(path.len(), 4, "Should take 4 base pair insertions to close the stem");
                // Check if the final path actually formed the structure
                // We know it starts empty, so all moves should be insertions
                for m in path {
                    assert!(m.is_insertion, "All moves should be insertions for this test case");
                }
            },
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
        let result = _greedy_find_path(&model, seq, s1, s2);

        match result {
            Ok((path, saddle, barrier)) => {
                println!("Path found with {} steps", path.len());
                println!("Saddle Energy: {:.2} kcal/mol", saddle);
                println!("Barrier Height: {:.2} kcal/mol", barrier);
                
                // Assertions
                assert_eq!(path.len(), 7, "Should take 7 base pair insertions to close the stem");
                // Check if the final path actually formed the structure
                // We know it starts empty, so all moves should be insertions
                for m in path{
                    println!("{:?}", m)
                }
                
            },
            Err(e) => panic!("Greedy path search failed: {}", e),
        }
    }

    #[test]
    fn test_greedy_find_path_printing(){
        let seq = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let s1 = ".(((..............)))....";
        let s2 = "..((((.........))))......";
        
        // Initialize Model (This might fail if params aren't found in your env)
        // Ensure you have valid params or a mock model available.
        let model = ViennaRNA::default();

        let (steps, stats) = greedy_find_path(&model, seq, s1, s2).unwrap();
        println!("Test steps:");
        for step in steps {
            println!("{} \t {} \t {} kcal/mol", step.structure,step.move_applied.unwrap_or_default(), step.energy );
        } 
        println!("Stats: {:?}", stats);
    }
}


