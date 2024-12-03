use std::fmt::Debug;
use std::ops::Index;

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
    mismatch: usize,
    gap_opening: usize,
    gap_extension: usize,
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

    fn new() -> Option<Self> {
        Some(Self::default())
    }
}

impl<'a> Debug for WaveFrontElement<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "element, {:#?}, {:#?}", self.state, self.offset)?;
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

    fn with_hi_lo(hi: i32, lo: i32) -> Self {
        Self {
            hi,
            lo,
            elements: Vec::default(),
        }
    }

    fn new(
        state: State,
        other_gap_open: Option<&WaveFront>,
        other_gap_extend: Option<&WaveFront>,
        other_mismatch: Option<&WaveFront>,
    ) -> Option<Self> {
        Some(Self::default())
    }

    fn get_element(&self, idx: i32) -> Option<&WaveFrontElement<'a>> {
        self.elements.get((idx - self.lo) as usize)
    }

    fn get_offset(&self, idx: i32) -> Option<&i32> {
        if let Some(wf) = self.elements.get((idx - self.lo) as usize) {
            Some(&wf.offset)
        } else {
            None
        }
    }
}

impl<'a> Index<i32> for WaveFront<'a> {
    type Output = WaveFrontElement<'a>;
    fn index(&self, index: i32) -> &Self::Output {
        &self.elements[(index - self.lo) as usize]
    }
}

#[derive(Default, Debug)]
struct WaveFrontTensor<'a> {
    i: Option<WaveFront<'a>>,
    d: Option<WaveFront<'a>>,
    m: Option<WaveFront<'a>>,
}

impl<'a> WaveFrontTensor<'a> {
    fn new(
        other_gap_open: Option<&WaveFrontTensor>,   // s - o - e
        other_gap_extend: Option<&WaveFrontTensor>, // s - e
        other_mismatch: Option<&WaveFrontTensor>,   // s - x
    ) -> Option<Self> {
        let hi: i32 = [
            other_gap_open.and_then(|r| r.m.as_ref().map(|r| r.hi)),
            other_mismatch.and_then(|r| r.m.as_ref().map(|r| r.hi)),
            other_gap_extend.and_then(|r| r.i.as_ref().map(|r| r.hi)),
            other_gap_extend.and_then(|r| r.d.as_ref().map(|r| r.hi)),
        ]
        .iter()
        .filter_map(|&opt| opt)
        .max()?
            + 1;

        let lo: i32 = [
            other_gap_open.and_then(|r| r.m.as_ref().map(|r| r.lo)),
            other_mismatch.and_then(|r| r.m.as_ref().map(|r| r.lo)),
            other_gap_extend.and_then(|r| r.i.as_ref().map(|r| r.lo)),
            other_gap_extend.and_then(|r| r.d.as_ref().map(|r| r.lo)),
        ]
        .iter()
        .filter_map(|&opt| opt)
        .min()?
            - 1;
        println!("lo: {}, hi: {}", lo, hi);

        // TODO: calculate lo and hi dynamically, as different wavefronts might have different boundaries?. also use I and D for M

        let (mut i, mut d, mut m) = (
            WaveFront::with_hi_lo(hi, lo),
            WaveFront::with_hi_lo(hi, lo),
            WaveFront::with_hi_lo(hi, lo),
        );

        // iterate over all elemetns in teh corresponding wavefrotns and append to i, d, m. need to do this safely (handel none values)

        for idx in lo..=hi {
            i.elements.push(WaveFrontElement {
                offset: *[
                    other_gap_open
                        .and_then(|r| r.m.as_ref().map(|r| r.get_offset(idx - 1).unwrap_or(&-1))),
                    other_gap_extend
                        .and_then(|r| r.i.as_ref().map(|r| r.get_offset(idx - 1).unwrap_or(&-1))),
                ]
                .iter()
                .filter_map(|&opt| opt)
                .max()
                .unwrap_or(&-1),
                parents: Vec::default(),
                state: State::I,
            });
            d.elements.push(WaveFrontElement {
                offset: *[
                    other_gap_open
                        .and_then(|r| r.m.as_ref().map(|r| r.get_offset(idx + 1).unwrap_or(&-1))),
                    other_gap_extend
                        .and_then(|r| r.d.as_ref().map(|r| r.get_offset(idx + 1).unwrap_or(&-1))),
                ]
                .iter()
                .filter_map(|&opt| opt)
                .max()
                .unwrap_or(&-1),
                parents: Vec::default(),
                state: State::D,
            });
            m.elements.push(WaveFrontElement {
                offset: *[
                    other_mismatch
                        .and_then(|r| r.m.as_ref().map(|r| r.get_offset(idx + 1).unwrap_or(&-1))),
                    i.get_offset(idx),
                    d.get_offset(idx),
                ]
                .iter()
                .filter_map(|&opt| opt)
                .max()
                .unwrap_or(&-1),
                parents: Vec::default(),
                state: State::M,
            });
        }

        println!("{:#?}", m);
        println!("{:#?}", d);
        println!("{:#?}", i);

        m.elements.retain(|e| e.offset != -1);
        i.elements.retain(|e| e.offset != -1);
        d.elements.retain(|e| e.offset != -1);

        Some(Self {
            i: if !i.elements.is_empty() {
                Some(i)
            } else {
                None
            },
            d: if !d.elements.is_empty() {
                Some(d)
            } else {
                None
            },
            m: if !m.elements.is_empty() {
                Some(m)
            } else {
                None
            },
        })
    }
}

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
                let s = wfs.len();
                wfs.push(WaveFrontTensor::new(
                    wfs.get((s - SCHEME.gap_opening - SCHEME.gap_extension) as usize) // assuming this does not go below 0. need to check
                        .and_then(|r| r.as_ref()),
                    wfs.get((s - SCHEME.gap_extension) as usize)
                        .and_then(|r| r.as_ref()),
                    wfs.get((s - SCHEME.mismatch) as usize)
                        .and_then(|r| r.as_ref()),
                ));
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
    #[test]
    fn test_wavefront_tensor_new_all_none() {
        let tensor = WaveFrontTensor::new(None, None, None);
        assert!(
            tensor.is_none(),
            "Tensor should be None when all inputs are None"
        );
    }

    #[test]
    fn test_wavefront_tensor_new_single_input() {
        // Mock input tensor with known hi and lo
        let mock_tensor = WaveFrontTensor {
            i: Some(WaveFront {
                hi: -1,
                lo: 2,
                elements: vec![
                    WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::I,
                    };
                    4
                ],
            }),
            d: Some(WaveFront {
                hi: -2,
                lo: 3,
                elements: vec![WaveFrontElement {
                    offset: 1,
                    parents: Vec::default(),
                    state: State::D,
                }],
            }),
            m: Some(WaveFront {
                lo: -2,
                hi: 3,
                elements: vec![
                    WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::I
                    };
                    6
                ],
            }),
        };

        let result = WaveFrontTensor::new(Some(&mock_tensor), None, None).unwrap();

        // Assert hi and lo
        assert_eq!(result.m.as_ref().unwrap().lo, -3, "lo should decrease by 1");
        assert_eq!(result.m.as_ref().unwrap().hi, 4, "hi should increase by 1");

        // Assert the range is correctly initialized
        let expected_size = 8; // From -3 to 4
        assert_eq!(
            result.m.as_ref().unwrap().elements.len(),
            expected_size,
            "Wavefront I should cover the correct range"
        );
    }
    #[test]
    fn test_wavefront_tensor_new_full_inputs() {
        return;
        //TODO
        // Mock input tensors
        let gap_open = WaveFrontTensor {
            m: Some(WaveFront {
                lo: -5,
                hi: 2,
                elements: vec![
                    WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::M
                    };
                    8
                ],
            }),
            i: Some(WaveFront {
                lo: -5,
                hi: 2,
                elements: vec![
                    WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::I
                    };
                    8
                ],
            }),
            d: Some(WaveFront {
                lo: -5,
                hi: 2,
                elements: vec![
                    WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::D
                    };
                    8
                ],
            }),
        };
        let gap_extend = WaveFrontTensor {
            i: Some(WaveFront {
                lo: -3,
                hi: 6,
                elements: vec![],
            }),
            d: None,
            m: None,
        };

        let result = WaveFrontTensor::new(Some(&gap_open), Some(&gap_extend), None).unwrap();

        // Assert hi and lo
        assert_eq!(
            result.i.as_ref().unwrap().lo,
            -6,
            "lo should be minimum across inputs minus 1"
        );
        assert_eq!(
            result.i.as_ref().unwrap().hi,
            7,
            "hi should be maximum across inputs plus 1"
        );

        // Assert the range
        let expected_size = 14; // From -6 to 7
        assert_eq!(
            result.i.as_ref().unwrap().elements.len(),
            expected_size,
            "Wavefront I should cover the correct range"
        );
    }
    #[test]
    fn test_wavefront_tensor_new_mixed_inputs() {
        return;
        //TODO
        let gap_open = WaveFrontTensor {
            m: Some(WaveFront {
                lo: -1,
                hi: 1,
                elements: vec![],
            }),
            i: None,
            d: None,
        };

        let result = WaveFrontTensor::new(Some(&gap_open), None, None).unwrap();

        // Assert hi and lo
        assert_eq!(result.i.as_ref().unwrap().lo, -2, "lo should decrease by 1");
        assert_eq!(result.i.as_ref().unwrap().hi, 2, "hi should increase by 1");

        // Check `m` is properly initialized
        assert!(result.d.is_none(), "D state should remain None");
    }
    #[test]
    fn test_wavefront_tensor_new_negative_range() {
        return;
        //TODO
        let mismatch = WaveFrontTensor {
            m: Some(WaveFront {
                lo: -10,
                hi: -5,
                elements: vec![],
            }),
            i: None,
            d: None,
        };

        let result = WaveFrontTensor::new(None, None, Some(&mismatch)).unwrap();

        // Assert hi and lo
        assert_eq!(
            result.m.as_ref().unwrap().lo,
            -11,
            "lo should decrease by 1"
        );
        assert_eq!(result.m.as_ref().unwrap().hi, -4, "hi should increase by 1");

        // Assert range is handled properly
        let expected_size = 8; // From -11 to -4
        assert_eq!(
            result.m.as_ref().unwrap().elements.len(),
            expected_size,
            "Wavefront M should cover the correct range"
        );
    }
}
