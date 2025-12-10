use ff_structure::PairTable;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;



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

fn int_hash_64(i: isize, j: isize) -> u64 {
    let mut hasher = DefaultHasher::new();
    (i, j).hash(&mut hasher);
    hasher.finish()
}


pub struct Intermediate {
    pub pt: PairTable,
    pub saddle_energy: f64,
    pub current_energy: f64,
    // We track which moves from the global list are still valid
    pub remaining_moves: Vec<usize>, 
}



fn compare_structures(pt1: &PairTable, pt2: &PairTable) -> StructureDifference {
    let mut diff = StructureDifference {
        move_list: vec![],
        hash_list: vec![],
        bp_distance: 0,
    };

    let length = pt1.len();

    for i in 0..length {
        if pt1[i] != pt2[i] {
            if i < pt1[i].expect("PairTable is empty! (has None value") {
                diff.move_list.push((-i, -pt1[i]));
                diff.hash_list.push();
                diff.bp_distance += 1;
            }
        }
    }
    
    
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


}

