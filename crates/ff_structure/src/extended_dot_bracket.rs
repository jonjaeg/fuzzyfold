//! Extended-dot-bracket notation.
//! 
//! This keeps the canonical dot-bracket implementation seperate from the extended-dot-bracket.
//! Especially for Pseudoknots, we need more characters for nested structures.
//! This implementation uses two seperate enums:
//! - ExtendedDotBracket
//! - BracketKind

use crate::DotBracket;
use crate::StructureError;



/// Represents a single character in extended dot-bracket notation,
/// covering simple and pseudoknot-encoding bracket types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ExtendedDotBracket {
    Unpaired,           // '.'
    Break,              // '+' or '&'
    Open(BracketKind),
    Close(BracketKind),
}

/// The bracket type, ordered by conventional pseudoknot nesting level.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BracketKind {
    Round,    // ( )
    Square,   // [ ]
    Curly,    // { }
    Angle,    // < >
    UpperA,   // A a
    UpperB,   // B b
    UpperC,   // C c
    UpperD,   // D d
}

impl BracketKind {
    pub fn open_char(self) -> char {
        match self {
            BracketKind::Round  => '(',
            BracketKind::Square => '[',
            BracketKind::Curly  => '{',
            BracketKind::Angle  => '<',
            BracketKind::UpperA => 'A',
            BracketKind::UpperB => 'B',
            BracketKind::UpperC => 'C',
            BracketKind::UpperD => 'D',
        }
    }

    pub fn close_char(self) -> char {
        match self {
            BracketKind::Round  => ')',
            BracketKind::Square => ']',
            BracketKind::Curly  => '}',
            BracketKind::Angle  => '>',
            BracketKind::UpperA => 'a',
            BracketKind::UpperB => 'b',
            BracketKind::UpperC => 'c',
            BracketKind::UpperD => 'd',
        }
    }

    /// All kinds in their conventional pseudoknot-nesting order.
    pub fn all() -> &'static [BracketKind] {
        &[
            BracketKind::Round,
            BracketKind::Square,
            BracketKind::Curly,
            BracketKind::Angle,
            BracketKind::UpperA,
            BracketKind::UpperB,
            BracketKind::UpperC,
            BracketKind::UpperD,
        ]
    }
}

impl TryFrom<char> for ExtendedDotBracket {
    type Error = StructureError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            '.' => return Ok(ExtendedDotBracket::Unpaired),
            '+' | '&' => return Ok(ExtendedDotBracket::Break),
            _ => {}
        }
        for &kind in BracketKind::all() {
            if c == kind.open_char()  { return Ok(ExtendedDotBracket::Open(kind));  }
            if c == kind.close_char() { return Ok(ExtendedDotBracket::Close(kind)); }
        }
        Err(StructureError::InvalidToken(
            format!("character '{}'", c),
            "extended dot-bracket".to_string(),
            0,
        ))
    }
}

impl From<ExtendedDotBracket> for char {
    fn from(edb: ExtendedDotBracket) -> Self {
        match edb {
            ExtendedDotBracket::Unpaired     => '.',
            ExtendedDotBracket::Break        => '+',
            ExtendedDotBracket::Open(kind)   => kind.open_char(),
            ExtendedDotBracket::Close(kind)  => kind.close_char(),
        }
    }
}

/// Convenience: lossless downcast when only Round brackets are present.
impl TryFrom<ExtendedDotBracket> for DotBracket {
    type Error = StructureError;

    fn try_from(edb: ExtendedDotBracket) -> Result<Self, Self::Error> {
        match edb {
            ExtendedDotBracket::Unpaired              => Ok(DotBracket::Unpaired),
            ExtendedDotBracket::Break                 => Ok(DotBracket::Break),
            ExtendedDotBracket::Open(BracketKind::Round)  => Ok(DotBracket::Open),
            ExtendedDotBracket::Close(BracketKind::Round) => Ok(DotBracket::Close),
            other => Err(StructureError::InvalidToken(
                format!("character '{}'", char::from(other)),
                "simple dot-bracket".to_string(),
                0,
            )),
        }
    }
}