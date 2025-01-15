use std::fmt::Debug;
use std::ops::Index;

use crate::errors::{AlignerError, Result};
use crate::parse::{Mode, Record};

//TODO: rewrite to use less options (optimally). alternatively use Vec<Option<WaveFrontElement>>> in current implementation

type WaveFrontVec<'a> = Vec<Option<WaveFrontTensor<'a>>>;

const SCHEME: ScoringScheme = ScoringScheme {
    mismatch: 4,
    gap_opening: 2,
    gap_extension: 6,
};

pub fn wfa_align<'a>(seq1: &Record, seq2: &Record, mode: Mode) -> Result<'a, ()> {
    let mut ocean = match mode {
        Mode::Global => Ocean::global(),
        _ => return Err(AlignerError::AlignmentError("not implemented")),
    };
    while !ocean.is_converged(&seq1.seq, &seq2.seq) {
        ocean.expand(&seq1.seq, &seq2.seq)?;
    }
    let s = if let Ocean::Global(w) = &ocean {
        w.len()
    } else {
        0
    };
    println!("converged with score {}: {:#?}", s, ocean);
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

#[derive(Debug, Default, PartialEq, Eq)]
struct HiLoTracker {
    cur_hi: i32,
    cur_lo: i32,
    lo_set: bool,
}

impl HiLoTracker {
    fn with_lo_hi(lo: i32, hi: i32) -> Self {
        Self {
            cur_hi: hi,
            cur_lo: lo,
            lo_set: false,
        }
    }
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
    /*
    fn new() -> Option<Self> {
        Some(Self::default())
    }
    */
}

impl<'a> Debug for WaveFrontElement<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "element, {:#?}, {:#?}", self.state, self.offset)?;
        Ok(())
    }
}

#[derive(Default, Debug, PartialEq, Eq)]
struct WaveFront<'a> {
    // are hi and lo inclusive?? yes
    hi: i32,
    lo: i32,
    elements: Vec<Option<WaveFrontElement<'a>>>,
}

impl<'a> WaveFront<'a> {
    fn expand(&mut self, seq1: &[u8], seq2: &[u8]) {
        for (i, element) in self.elements.iter_mut().enumerate() {
            if let Some(ref mut e) = element {
                while seq1[e.y(self.lo + i as i32)] == seq2[e.x(self.lo + i as i32)] {
                    e.offset += 1;
                }
            } else {
                continue;
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
    /*
    fn new(
        state: State,
        other_gap_open: Option<&WaveFront>,
        other_gap_extend: Option<&WaveFront>,
        other_mismatch: Option<&WaveFront>,
    ) -> Option<Self> {
        Some(Self::default())
    }
    */

    fn get_element(&self, idx: i32) -> Option<&WaveFrontElement<'a>> {
        self.elements
            .get((idx - self.lo) as usize)
            .and_then(|e| e.as_ref())
    }

    fn get_offset(&self, idx: i32) -> Option<&i32> {
        if let Some(wf) = self.get_element(idx) {
            Some(&wf.offset)
        } else {
            None
        }
    }

    fn is_converged(&self, seq1: &[u8], seq2: &[u8]) -> bool {
        self.elements
            .iter()
            .enumerate()
            .filter_map(|(i, element)| {
                element
                    .as_ref()
                    .map(|e| (e.x(self.lo + i as i32), e.y(self.lo + i as i32)))
            })
            .any(|(x, y)| x == seq1.len() && y == seq2.len())
    }
}

impl<'a> Index<i32> for WaveFront<'a> {
    type Output = Option<WaveFrontElement<'a>>;
    fn index(&self, index: i32) -> &Self::Output {
        &self.elements[(index - self.lo) as usize]
    }
}

#[derive(Default, Debug, Eq, PartialEq)]
struct WaveFrontTensor<'a> {
    i: Option<WaveFront<'a>>,
    d: Option<WaveFront<'a>>,
    m: Option<WaveFront<'a>>,
}

impl<'a> WaveFrontTensor<'a> {
    fn expand(&mut self, seq1: &[u8], seq2: &[u8]) {
        if let Some(ref mut wf) = self.m {
            wf.expand(seq1, seq2);
        }
    }

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
        // the above should be correct
        // TODO: calculate lo and hi dynamically, as different wavefronts might have different boundaries?(should not be necessary, as lo, hi = min, max -> others are at most hi, lo). also use I and D for M

        let (mut i, mut d, mut m) = (
            WaveFront::with_hi_lo(hi, lo),
            WaveFront::with_hi_lo(hi, lo),
            WaveFront::with_hi_lo(hi, lo),
        );

        let (mut i_tracker, mut d_tracker, mut m_tracker) = (
            HiLoTracker::with_lo_hi(lo, hi),
            HiLoTracker::with_lo_hi(lo, hi),
            HiLoTracker::with_lo_hi(lo, hi),
        );

        // iterate over all elemetns in teh corresponding wavefrotns and append to i, d, m. need to do this safely (handel none values)

        for idx in lo..=hi {
            //TODO make sure to push NONE values also. Additionally fux lo, hi calcs for I, D, M

            let offset: Option<i32> = [
                other_gap_open
                    .and_then(|r| r.m.as_ref().map(|r| r.get_offset(idx + 1)))
                    .and_then(|o| o),
                other_gap_extend
                    .and_then(|r| r.d.as_ref().map(|r| r.get_offset(idx + 1)))
                    .and_then(|o| o),
            ]
            .into_iter()
            .max()
            .flatten()
            .copied();
            if let Some(offset) = offset {
                d.elements.push(Some(WaveFrontElement {
                    offset,
                    parents: Vec::default(),
                    state: State::D,
                }));

                d_tracker.cur_hi = idx;
                if !d_tracker.lo_set {
                    d_tracker.cur_lo = idx;
                    d_tracker.lo_set = true;
                }
            } else {
                d.elements.push(None);
            }

            let offset: Option<i32> = [
                other_gap_open
                    .and_then(|r| r.m.as_ref().map(|r| r.get_offset(idx - 1)))
                    .flatten(),
                other_gap_extend
                    .and_then(|r| r.i.as_ref().map(|r| r.get_offset(idx - 1)))
                    .flatten(),
            ]
            .into_iter()
            .max()
            .flatten()
            .copied();
            if let Some(offset) = offset {
                i.elements.push(Some(WaveFrontElement {
                    offset: offset + 1,
                    parents: Vec::default(),
                    state: State::I,
                }));

                i_tracker.cur_hi = idx;
                if !i_tracker.lo_set {
                    i_tracker.cur_lo = idx;
                    i_tracker.lo_set = true;
                }
            } else {
                i.elements.push(None);
            }
            let offset: Option<i32> = [
                other_mismatch
                    .and_then(|r| r.m.as_ref().map(|r| r.get_offset(idx).map(|o| o + 1)))
                    .flatten(),
                i.get_offset(idx).copied(),
                d.get_offset(idx).copied(),
            ]
            .into_iter()
            .max()
            .flatten();
            if let Some(offset) = offset {
                m.elements.push(Some(WaveFrontElement {
                    offset,
                    parents: Vec::default(),
                    state: State::M,
                }));

                m_tracker.cur_hi = idx;
                if !m_tracker.lo_set {
                    m_tracker.cur_lo = idx;
                    m_tracker.lo_set = true;
                }
            } else if m_tracker.lo_set {
                m.elements.push(None);
            }
        }

        (i.lo, i.hi) = (i_tracker.cur_lo, i_tracker.cur_hi);
        (d.lo, d.hi) = (d_tracker.cur_lo, d_tracker.cur_hi);
        (m.lo, m.hi) = (m_tracker.cur_lo, m_tracker.cur_hi);

        i.elements.rotate_left(lo.abs_diff(i.lo) as usize);
        i.elements.truncate(i.hi.abs_diff(i.lo) as usize + 1);
        d.elements.rotate_left(lo.abs_diff(d.lo) as usize);
        d.elements.truncate(d.hi.abs_diff(d.lo) as usize + 1);
        m.elements.truncate(m.hi.abs_diff(m.lo) as usize + 1);
        /*
        println!("{:#?}", m);
        println!("{:#?}", d);
        println!("{:#?}", i);
        */
        Some(Self {
            i: if i_tracker.lo_set { Some(i) } else { None },
            d: if d_tracker.lo_set { Some(d) } else { None },
            m: if m_tracker.lo_set { Some(m) } else { None },
        })
    }

    fn is_converged(&self, seq1: &[u8], seq2: &[u8]) -> bool {
        if let Some(i) = &self.i {
            if i.is_converged(seq1, seq2) {
                return true;
            }
        }
        if let Some(d) = &self.d {
            if d.is_converged(seq1, seq2) {
                return true;
            }
        }
        if let Some(m) = &self.m {
            if m.is_converged(seq1, seq2) {
                return true;
            }
        }
        false
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
                elements: vec![Some(WaveFrontElement {
                    offset: 0,
                    state: State::M,
                    parents: Vec::default(),
                })],
            }),
        })])
    }

    fn expand(&mut self, seq1: &[u8], seq2: &[u8]) -> Result<'a, ()> {
        match self {
            Ocean::Global(ref mut wfs) => {
                let s = wfs.len() as i32;
                wfs.push(WaveFrontTensor::new(
                    wfs.get((s - SCHEME.gap_opening - SCHEME.gap_extension) as usize) // assuming this does not go below 0. need to check. shoul not matter unless wf.len() ~ usize::MAX
                        .and_then(|r| r.as_ref()),
                    wfs.get((s - SCHEME.gap_extension) as usize)
                        .and_then(|r| r.as_ref()),
                    wfs.get((s - SCHEME.mismatch) as usize)
                        .and_then(|r| r.as_ref()),
                ));
                if let Some(Some(ref mut wf)) = wfs.get_mut(s as usize) {
                    wf.expand(seq1, seq2);
                }
            }
            Ocean::SemiGlobal(_) => return Err(AlignerError::AlignmentError("not implemented")),
            Ocean::Local(_) => return Err(AlignerError::AlignmentError("not implemented")),
        }
        Ok(())
    }

    fn is_converged(&self, seq1: &[u8], seq2: &[u8]) -> bool {
        if let Ocean::Global(v) = self {
            if let Some(Some(t)) = v.last() {
                return t.is_converged(seq1, seq2);
            }
        }
        false
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_init() {
        //TODO
    }
    #[test]
    fn test_expand() {
        //TODO
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
    fn recurrance_eq() {
        let m_tensor_full = WaveFrontTensor {
            i: Some(WaveFront {
                hi: -1,
                lo: 2,
                elements: vec![
                    Some(WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::I,
                    });
                    4
                ],
            }),
            d: Some(WaveFront {
                hi: -2,
                lo: 3,
                elements: vec![Some(WaveFrontElement {
                    offset: 1,
                    parents: Vec::default(),
                    state: State::D,
                })],
            }),
            m: Some(WaveFront {
                lo: -2,
                hi: 3,
                elements: vec![
                    Some(WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::I
                    });
                    6
                ],
            }),
        };
        let m_tensor_simple = WaveFrontTensor {
            i: None,
            d: None,
            m: Some(WaveFront {
                lo: -2,
                hi: 3,
                elements: vec![
                    Some(WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::I
                    });
                    6
                ],
            }),
        };

        let m_tensor_simple_gap = WaveFrontTensor {
            i: Some(WaveFront {
                hi: -1,
                lo: 2,
                elements: vec![
                    Some(WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::I,
                    });
                    4
                ],
            }),
            d: Some(WaveFront {
                hi: -2,
                lo: 3,
                elements: vec![Some(WaveFrontElement {
                    offset: 1,
                    parents: Vec::default(),
                    state: State::D,
                })],
            }),
            m: None,
        };

        assert_eq!(
            WaveFrontTensor::new(Some(&m_tensor_simple), None, None),
            WaveFrontTensor::new(Some(&m_tensor_full), None, None),
            "only M relevant in s - o - e"
        );
        assert_eq!(
            WaveFrontTensor::new(None, None, Some(&m_tensor_simple)),
            WaveFrontTensor::new(None, None, Some(&m_tensor_full)),
            "only M relevant in s - x"
        );
        assert_eq!(
            WaveFrontTensor::new(None, Some(&m_tensor_simple_gap), None),
            WaveFrontTensor::new(None, Some(&m_tensor_full), None),
            "only I, D relevant in s - e"
        );
    }

    #[test]
    fn test_initial() {
        // currently this test shows false hi and lo values in i and d, as well as the need for placefolders in WaveFront
        let initial = WaveFrontTensor {
            i: None,
            d: None,
            m: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: vec![Some(WaveFrontElement {
                    offset: 0,
                    parents: Vec::default(),
                    state: State::M,
                })],
            }),
        };

        let true_res_o = WaveFrontTensor {
            i: Some(WaveFront {
                hi: 1,
                lo: 1,
                elements: vec![Some(WaveFrontElement {
                    offset: 1,
                    parents: Vec::default(),
                    state: State::I,
                })],
            }),
            d: Some(WaveFront {
                hi: -1,
                lo: -1,
                elements: vec![Some(WaveFrontElement {
                    offset: 0,
                    parents: Vec::default(),
                    state: State::D,
                })],
            }),
            m: Some(WaveFront {
                hi: 1,
                lo: -1,
                elements: vec![
                    Some(WaveFrontElement {
                        offset: 0,
                        parents: Vec::default(),
                        state: State::M,
                    }),
                    None,
                    Some(WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::M,
                    }),
                ],
            }),
        };
        let true_res_m = WaveFrontTensor {
            i: None,
            d: None,
            m: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: vec![Some(WaveFrontElement {
                    offset: 1,
                    parents: Vec::default(),
                    state: State::M,
                })],
            }),
        };

        assert_eq!(
            WaveFrontTensor::new(Some(&initial), None, None),
            Some(true_res_o)
        );
        assert_eq!(
            WaveFrontTensor::new(None, None, Some(&initial)),
            Some(true_res_m)
        );
    }

    #[test]
    fn test_full() {
        //TODO
        let _initial_o_e = WaveFrontTensor {
            i: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
            d: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
            m: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
        };
        let _initial_e = WaveFrontTensor {
            i: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
            d: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
            m: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
        };
        let _initial_m = WaveFrontTensor {
            i: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
            d: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
            m: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
        };

        let _res = WaveFrontTensor {
            i: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
            d: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
            m: Some(WaveFront {
                hi: 0,
                lo: 0,
                elements: Vec::default(),
            }),
        };
        /*
        assert_eq!(
            Some(res),
            WaveFrontTensor::new(Some(&initial_o_e), Some(&initial_e), Some(&initial_m))
        );
        */
    }

    #[test]
    fn test_iteration() {
        //TODO
        let mut initial = Ocean::global();
        let query = b"AAAATTTTCCCC";
        let db = b"AAAATCTCC";
        assert!(initial.expand(query, db).is_ok());
        //assert_eq!(initial, ...)
        assert!(initial.expand(query, db).is_ok());
        //assert_eq!(initial, ...)
        assert!(initial.expand(query, db).is_ok());
        //assert_eq!(initial, ...)
        assert!(initial.expand(query, db).is_ok());
        //assert_eq!(initial, ...)
        assert!(initial.expand(query, db).is_ok());
        //assert_eq!(initial, ...)
        assert!(initial.expand(query, db).is_ok());
        //assert_eq!(initial, ...)
    }

    #[test]
    fn test_converge() {
        let initial = Ocean::global();
        let query = b"AACATCAY";
        let db = b"ATAGTAG";
        assert!(!initial.is_converged(query, db))
    }
}
