use ff_structure::PairTable;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::fmt;



/// Move represent the insertion or deletion of a base pair (i, j) between two structures. 
/// If is_insertion is true, it represents the insertion of the base pair (i, j) into the structure.
/// If is_insertion is false, it represents the deletion of the base pair (i, j) from the structure.
/// Move corresponds to an elementary move in the transformation from two consecutive intermediate structures S_t to S_t+1.
/// For example, if we have a move with i = 2, j = 5, is_insertion = true, it means that we are inserting the base pair (2, 5) into the structure.
/// The corresponding PairTable representation would have pt[2] = Some(5) and pt[5] = Some(2).
/// We define the base pair distance between two structures as the number of moves required to transform one structure into the other.
/// This means one move corresponds to the change of two indices in the PairTable representation.
#[derive(Debug, Clone)]
pub struct Move {
    pub i: usize,
    pub j: usize,
    pub is_insertion: bool, // true = insert, false = delete
}


/// StructureDifference captures the difference between two structures in terms of:
/// - move_list: A list of moves (insertions or deletions of base pairs) required to transform one structure into the other.
/// - hash_list: A list of unique hashes corresponding to each move in the move_list. These hashes are computed based on the indices of the base pairs involved in the moves.
/// - bp_distance: The base pair distance between the two structures, defined as the number of base pairs that differ between them.
/// NOTE: Each move in the move_list corresponds to a change of two indices in the PairTable representation,
/// hence a bp_distance of one equals one move.
#[derive(Debug, Clone)]
pub struct StructureDifference {
    pub move_list: Vec<Move>, 
    pub hash_list: Vec<u64>,
    pub bp_distance: u32,
}

// I hash a base pair (i, j) into a u64 value
// do i need also the boolean flag is_insertion?
// because i can not differentiate between insertion and deletion only with (i, j)
// --> i check in the function that generates the StructureDifference between two structures

/// Helper to hash indices.
/// To match Python's logic:
/// - Deletions use negative values (-i, -j)
/// - Insertions use positive values (i, j)
fn int_hash_64(i: i64, j: i64) -> u64 {
    let mut hasher = DefaultHasher::new();
    (i, j).hash(&mut hasher);
    hasher.finish()
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
        // don't print the number of moves and hashes because they are equal to bp_distance
        //writeln!(f, "moves ({}):", self.move_list.len())?;
        for mv in &self.move_list {
            writeln!(f, "  - {}", mv)?;
        }
        writeln!(f, "hashes ({}):", self.hash_list.len())?;
        for h in &self.hash_list {
            writeln!(f, "  - {}", h)?;
        }
        Ok(())
    }
}

pub struct Intermediate {
    pub pt: PairTable,
    pub saddle_energy: f64,
    pub current_energy: f64,
    // We track which moves from the global list are still valid
    pub remaining_moves: Vec<Move>,
}

/// compare_structures returns the ElemenrtaryMoves required to transform pt1 into pt2.
/// The data structure StructureDifference captures the list of moves, their corresponding hashes, and the base pair distance.
/// TODO: check again if it is right
/// Compares two PairTables and returns the moves required to transform pt1 into pt2.
pub fn compare_structures(pt1: &PairTable, pt2: &PairTable) -> StructureDifference {
    let mut diff = StructureDifference {
        move_list: vec![],
        hash_list: vec![],
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

                    // Python behavior: hash (-i, -j) for deletions
                    let h = int_hash_64(-(i as i64), -(j as i64));
                    
                    diff.move_list.push(deletion_move);
                    diff.hash_list.push(h);
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

                    // Python behavior: hash (i, j) for insertions
                    let h = int_hash_64(i as i64, j as i64);

                    diff.move_list.push(insertion_move);
                    diff.hash_list.push(h);
                    diff.bp_distance += 1;
                }
            }
        }
    }

    // Python does `move_list.sort(key=lambda x: x[0])`
    // We sort by the index `i` to maintain consistent order
    // TODO: is sorting useful, because if a insertion of a base pair (i, j) is before a deletion that affects the same indices,
    // then the order matters? ==> maybe other function solves this.
    diff.move_list.sort_by_key(|m| m.i);

    diff
}


#[cfg(test)]
mod tests {
    // when importing ff_structure::PairTable into the utils.rs (outside of the tests module), 
    // This warning happens because PairTable is only used inside your tests module.
    // When you run cargo build (or when Rust Analyzer checks your code), 
    // the #[cfg(test)] module is ignored/stripped out. This leaves the import use ff_structure::PairTable; sitting at the top level with no one using it, hence the "unused import" warning.
    
    use ff_structure::PairTable;

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
    fn test_compare_structures() {
        let pt1 = PairTable::try_from(".(((..............)))....").unwrap();
        let pt2 = PairTable::try_from("..((((.........))))......").unwrap();
        println!("pt1: {}", pt1);
        println!("pt2: {}", pt2);
        let diff = super::compare_structures(&pt1, &pt2);
        assert_eq!(diff.bp_distance, 7);
        //assert_eq!(diff.move_list)
        println!("{}", diff);
    }
}