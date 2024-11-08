use std::fmt::Debug;

use crate::errors::{AStarError, Result};
use crate::parse::{Mode, Record};

type WaveFrontVec<'a> = Vec<Option<WaveFrontTensor<'a>>>;

const SCHEME: ScoringScheme = ScoringScheme {
    mismatch: 4,
    gap_opening: 2,
    gap_extension: 6,
};

pub fn wfa_align<'a>(seq1: &Record, seq2: &Record, mode: Mode) -> Result<'a, ()> {
    let mut ocean = match mode {
        Mode::Global => Ocean::global(),
        _ => return Err(AStarError::AlignmentError("not implemented")),
    };
    ocean.expand(&seq1.seq, &seq2.seq)?;
    Ok(())
}
#[derive(Default, Debug, Clone, PartialEq, PartialOrd, Ord, Eq)]
enum State {
    #[default]
    M,
    D,
    I,
}

struct ScoringScheme {
    mismatch: i32,
    gap_opening: i32,
    gap_extension: i32,
}

#[derive(Default, PartialEq, Eq, PartialOrd, Ord, Clone)]
struct WaveFrontElement<'a> {
    offset: i32,
    // backtrace info
    parents: Vec<&'a WaveFrontElement<'a>>,
    state: State,
}

impl<'a> WaveFrontElement<'a> {
    fn x(&self, diag: i32) -> usize {
        (self.offset - diag.min(0)) as usize
    }
    fn y(&self, diag: i32) -> usize {
        (self.offset + diag.max(0)) as usize
    }
}

impl<'a> Debug for WaveFrontElement<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Ok(())
    }
}

#[derive(Default, Debug)]
struct WaveFront<'a> {
    hi: i32,
    lo: i32,
    elements: Vec<WaveFrontElement<'a>>,
}

impl<'a> WaveFront<'a> {
    fn expand(&mut self, seq1: &[u8], seq2: &[u8]) {
        for (i, element) in self.elements.iter_mut().enumerate() {
            while seq1[element.y(self.lo + i as i32)] == seq2[element.x(self.lo + i as i32)] {
                element.offset += 1;
            }
        }
    }
}

#[derive(Default, Debug)]
struct WaveFrontTensor<'a> {
    i: Option<WaveFront<'a>>,
    d: Option<WaveFront<'a>>,
    m: Option<WaveFront<'a>>,
}

impl<'a> WaveFrontTensor<'a> {}

#[derive(Debug)]
enum Ocean<'a> {
    Global(WaveFrontVec<'a>),
    Local(WaveFrontVec<'a>),
    SemiGlobal(WaveFrontVec<'a>),
}

impl<'a> Ocean<'a> {
    fn global() -> Self {
        Self::Global(vec![Some(WaveFrontTensor {
            i: None,
            d: None,
            m: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: vec![WaveFrontElement {
                    offset: 0,
                    state: State::M,
                    parents: Vec::default(),
                }],
            }),
        })])
    }

    fn expand(&mut self, seq1: &[u8], seq2: &[u8]) -> Result<'a, ()> {
        match self {
            Ocean::Global(ref mut wfs) => {
                let s = wfs.len() as i32;
                if (s < SCHEME.mismatch || wfs[(s - SCHEME.mismatch) as usize].is_none())
                    && (s < SCHEME.gap_extension
                        || wfs[(s - SCHEME.gap_extension) as usize].is_none())
                    && (s < SCHEME.gap_opening - SCHEME.gap_extension
                        || wfs[(s - SCHEME.gap_opening - SCHEME.gap_extension) as usize].is_none())
                {
                    wfs.push(None);
                    return Ok(());
                }
            }
            Ocean::SemiGlobal(_) => return Err(AStarError::AlignmentError("not implemented")),
            Ocean::Local(_) => return Err(AStarError::AlignmentError("not implemented")),
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_init() {
        //todo!()
    }
    #[test]
    fn test_expand() {
        //todo!()
    }
}
