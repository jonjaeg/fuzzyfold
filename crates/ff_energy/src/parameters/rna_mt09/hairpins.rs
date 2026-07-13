use crate::Base;
use crate::parameters::parameterset::LoopEntry;

const fn convert<const N: usize>(bytes: &[u8; N]) -> [Base; N] {
    let mut out = [Base::A; N];
    let mut i = 0;
    while i < N {
        out[i] = match bytes[i] {
            b'A' => Base::A,
            b'C' => Base::C,
            b'G' => Base::G,
            b'U' => Base::U,
            _ => panic!("Invalid base"),
        };
        i += 1;
    }
    out
}

macro_rules! loop_entry {
    ($seq:literal, $val:expr) => {
        LoopEntry { seq: &convert($seq), val: $val }
    };
}

pub static TRILOOPS: &[LoopEntry] = &[];

pub static TETRALOOPS: &[LoopEntry] = &[
	loop_entry!(b"GGGGAC", -130),
	loop_entry!(b"GGUGAC", -62),
	loop_entry!(b"CGAAAG", -166),
	loop_entry!(b"GGAGAC", -168),
	loop_entry!(b"CGCAAG", -158),
	loop_entry!(b"GGAAAC", -184),
	loop_entry!(b"CGGAAG", -196),
	loop_entry!(b"CUUCGG", -233),
	loop_entry!(b"CGUGAG", -186),
	loop_entry!(b"CGAAGG", -139),
	loop_entry!(b"CUACGG", -142),
	loop_entry!(b"GGCAAC", -216),
	loop_entry!(b"CGCGAG", 20),
	loop_entry!(b"UGAGAG", -135),
	loop_entry!(b"CGAGAG", -102),
	loop_entry!(b"AGAAAU", -133),
	loop_entry!(b"CGUAAG", -97),
	loop_entry!(b"CUAACG", -117),
	loop_entry!(b"UGAAAG", -80),
	loop_entry!(b"GGAAGC", -67),
	loop_entry!(b"GGGAAC", -84),
	loop_entry!(b"UGAAAA", -52),
	loop_entry!(b"AGCAAU", -69),
	loop_entry!(b"AGUAAU", -68),
	loop_entry!(b"CGGGAG", -18),
	loop_entry!(b"AGUGAU", -63),
	loop_entry!(b"GGCGAC", -33),
	loop_entry!(b"GGGAGC", 6),
];

pub static HEXALOOPS: &[LoopEntry] = &[];

