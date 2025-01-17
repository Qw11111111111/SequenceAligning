use std::char;
use std::fmt::{Debug, Display};
use std::marker::PhantomData;
use std::ops::Index;

use crate::errors::{AlignerError, Result};
use crate::parse::{Mode, Record};

//TODO: rewrite to use less options (optimally). alternatively use Vec<Option<WaveFrontElement>>> in current implementation
//TODO: rewrite this entire shitfest

type WaveFrontVec<'a> = Vec<Option<WaveFrontTensor<'a>>>;

const MINLENGTH: u32 = 5;
const MAXDIFF: u32 = 20;

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
    while ocean.is_converged(&seq1.seq, &seq2.seq).is_none() {
        ocean.expand(&seq1.seq, &seq2.seq)?;
    }
    let s = if let Ocean::Global(w) = &ocean {
        w.len()
    } else {
        0
    };
    println!("converged with score {}: ", s);
    let t = ocean.traceback(&seq1.seq, &seq2.seq);
    println!("{}", t[0]);
    println!("{:#?}", t[0]);

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
    parents: Vec<State>,
    state: State,
    _marker: PhantomData<&'a u8>,
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
    fn get_distance(&self, seq1: &[u8], seq2: &[u8], diag: i32) -> i32 {
        //TODO: check if seqs are ordered correctly
        let left_v = seq1.len() as i32 - self.offset - diag;
        let left_h = seq2.len() as i32 - self.offset;
        left_v.max(left_h)
    }
}

impl<'a> Debug for WaveFrontElement<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Element {{")?;
        writeln!(
            f,
            "\tstate: {:#?}\n\toffset: {:#?}",
            self.state, self.offset
        )?;
        writeln!(f, "\tparents: {:#?}", self.parents)?;
        writeln!(f, "}}")?;
        Ok(())
    }
}

#[derive(Default, Debug, PartialEq, Eq, Clone)]
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
                while (e.y(self.lo + i as i32) < seq1.len() && e.x(self.lo + i as i32) < seq2.len())
                    && (seq1[e.y(self.lo + i as i32)] == seq2[e.x(self.lo + i as i32)])
                {
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
    fn get_state(&self, idx: i32) -> Option<&State> {
        if let Some(wf) = self.get_element(idx) {
            Some(&wf.state)
        } else {
            None
        }
    }

    fn is_converged(&self, seq1: &[u8], seq2: &[u8]) -> Option<&WaveFrontElement> {
        self.elements
            .iter()
            .enumerate()
            .filter_map(|(i, element)| {
                element
                    .as_ref()
                    .map(|e| (e, e.x(self.lo + i as i32), e.y(self.lo + i as i32)))
            })
            .find(|(_e, x, y)| *x == seq2.len() - 1 && *y == seq1.len() - 1)
            .map(|(e, _, _)| e)
    }
}

impl<'a> Index<i32> for WaveFront<'a> {
    type Output = Option<WaveFrontElement<'a>>;
    fn index(&self, index: i32) -> &Self::Output {
        &self.elements[(index - self.lo) as usize]
    }
}

fn get_parents<'a>(offset: i32, wfs: Vec<&'a WaveFrontElement>) -> Vec<State> {
    let mut parents = Vec::new();
    for wf in &wfs {
        if wf.offset == offset {
            parents.push(wf.state.clone());
        }
    }
    parents
}

#[derive(Default, Debug, Eq, PartialEq, Clone)]
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
                    parents: get_parents(
                        offset,
                        [
                            other_gap_open
                                .and_then(|r| r.m.as_ref().map(|r| r.get_element(idx + 1))),
                            other_gap_extend
                                .and_then(|r| r.d.as_ref().map(|r| r.get_element(idx + 1))),
                        ]
                        .iter()
                        .flatten()
                        .filter_map(|i| *i)
                        .collect(),
                    ),
                    offset,
                    state: State::D,
                    _marker: PhantomData,
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
                    parents: get_parents(
                        offset,
                        [
                            other_gap_open
                                .and_then(|r| r.m.as_ref().map(|r| r.get_element(idx - 1))),
                            other_gap_extend
                                .and_then(|r| r.i.as_ref().map(|r| r.get_element(idx - 1))),
                        ]
                        .iter()
                        .flatten()
                        .filter_map(|i| *i)
                        .collect(),
                    ),
                    state: State::I,
                    _marker: PhantomData,
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
                    parents: get_parents(
                        offset,
                        [
                            other_mismatch.and_then(|r| {
                                r.m.as_ref().map(|r| {
                                    r.get_element(idx).map(|f| WaveFrontElement {
                                        offset: f.offset + 1,
                                        parents: Vec::default(),
                                        state: State::M,
                                        _marker: PhantomData,
                                    })
                                })
                            }),
                            Some(i.get_element(idx).cloned()), //.as_ref(),
                            Some(d.get_element(idx).cloned()),
                        ]
                        .iter()
                        .flatten()
                        .filter_map(|i| i.as_ref())
                        .collect(),
                    ),
                    state: State::M,
                    _marker: PhantomData,
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

    fn is_converged(&self, seq1: &[u8], seq2: &[u8]) -> Option<&WaveFrontElement> {
        if let Some(i) = &self.i {
            if let Some(i) = i.is_converged(seq1, seq2) {
                return Some(i);
            }
        }
        if let Some(d) = &self.d {
            if let Some(d) = d.is_converged(seq1, seq2) {
                return Some(d);
            }
        }
        if let Some(m) = &self.m {
            if let Some(m) = m.is_converged(seq1, seq2) {
                return Some(m);
            }
        }
        None
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
                    _marker: PhantomData,
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
                self.trim(seq1, seq2);
            }
            Ocean::SemiGlobal(_) => return Err(AlignerError::AlignmentError("not implemented")),
            Ocean::Local(_) => return Err(AlignerError::AlignmentError("not implemented")),
        }
        Ok(())
    }

    fn trim(&mut self, seq1: &[u8], seq2: &[u8]) {
        //TODO fix: what happens if m.hi > i/d.hi, but m.lo > i.hi/d.hi (currently crash)
        let wfs = if let Self::Global(ref mut wfs) = self {
            wfs
        } else {
            return;
        };
        let current = if let Some(Some(wf)) = wfs.last_mut() {
            wf
        } else {
            return;
        };
        let m = if let Some(ref mut m) = current.m {
            m
        } else {
            return;
        };
        if m.lo.abs_diff(m.hi) <= MINLENGTH {
            return;
        }

        let mut min_d = 0;
        for diag in m.lo..=m.hi {
            if let Some(offset) = m.get_element(diag) {
                let dist = offset.get_distance(seq1, seq2, diag);
                min_d = min_d.min(dist);
            }
        }
        let mut next_d = m
            .elements
            .first()
            .expect("first element is ensured to be Some")
            .as_ref()
            .unwrap()
            .get_distance(seq1, seq2, m.lo);
        while m.lo < m.hi && next_d.abs_diff(min_d) > MAXDIFF {
            //println!("wtf {:#?}", m);
            m.lo += 1;
            m.elements.remove(0);
            while m.get_element(m.lo).is_none() {
                if m.lo == m.hi {
                    break;
                }
                //println!("wtf2 {:#?}", m);
                m.lo += 1;
                m.elements.remove(0);
            }
            next_d = m
                .elements
                .first()
                .expect("first element is ensured to be Some")
                .as_ref()
                .unwrap()
                .get_distance(seq1, seq2, m.lo);
        }
        let mut next_d = m
            .elements
            .last()
            .expect("first element is ensured to be Some")
            .as_ref()
            .unwrap()
            .get_distance(seq1, seq2, m.hi);
        while m.hi > m.lo && next_d.abs_diff(min_d) > MAXDIFF {
            m.hi -= 1;
            m.elements.pop();
            while m.get_element(m.hi).is_none() {
                if m.lo == m.hi {
                    break;
                }
                m.hi -= 1;
                m.elements.pop();
            }
            next_d = m
                .elements
                .last()
                .expect("first element is ensured to be Some")
                .as_ref()
                .unwrap()
                .get_distance(seq1, seq2, m.hi);
        }
        // println!("dhucwoierfgiyu4rewgcfyuoervtgv");

        // println!("{:#?}", current.i);
        // println!("{:#?}", current.d);
        if let Some(ref mut i) = current.i {
            // println!("{:#?}", i.elements);
            let t = if i.lo < m.lo {
                i.elements.rotate_left(i.lo.abs_diff(m.lo) as usize);
                i.lo.abs_diff(m.lo) as usize
                    + if i.hi > m.hi {
                        i.hi.abs_diff(m.hi) as usize
                    } else {
                        0
                    }
            } else if i.hi > m.hi {
                i.hi.abs_diff(m.hi) as usize
            } else {
                0
            };
            // println!("{}, {}, {}, {}", i.lo, i.hi, m.lo, m.hi);
            // println!("{}", t);
            // println!("{:#?}", i.elements);
            // println!("{}", i.elements.len() - t);
            i.elements.truncate(i.elements.len() - t);
            // println!("huhu");
            i.hi = i.hi.min(m.hi);
            i.lo = i.lo.max(m.lo);
        }
        if let Some(ref mut d) = current.d {
            let t = if d.lo < m.lo {
                // println!("heufhiwhufiwhfw");
                // println!("d: {:#?}", d);
                // println!("{}, {}", m.lo, m.hi);
                d.elements.rotate_left(d.lo.abs_diff(m.lo) as usize);
                // println!("well yes i am dumb");
                d.lo.abs_diff(m.lo) as usize
                    + if d.hi > m.hi {
                        d.hi.abs_diff(m.hi) as usize
                    } else {
                        0
                    }
            } else if d.hi > m.hi {
                d.hi.abs_diff(m.hi) as usize
            } else {
                0
            };
            // println!("{}, {}, {}, {}", d.lo, d.hi, m.lo, m.hi);
            // println!("{}", t);
            // println!("{:#?}", d.elements);
            d.elements.truncate(d.elements.len() - t);
            d.hi = d.hi.min(m.hi);
            d.lo = d.lo.max(m.lo);
        }
    }

    fn is_converged(&self, seq1: &[u8], seq2: &[u8]) -> Option<&WaveFrontElement> {
        if let Ocean::Global(v) = self {
            if let Some(Some(t)) = v.last() {
                return t.is_converged(seq1, seq2);
            }
        }
        None
    }

    fn traceback(&self, seq1: &[u8], seq2: &[u8]) -> Vec<Alignment> {
        let mut diag = seq1.len() as i32 - seq2.len() as i32;
        let mut states: Vec<Alignment> = vec![Alignment {
            seq1: Vec::new(),
            seq2: Vec::new(),
        }];
        let mut last = if let Some(last) = self.is_converged(seq1, seq2) {
            last
        } else {
            return Vec::default();
        };
        let l = if let Ocean::Global(l) = self {
            l.len()
        } else {
            0
        };
        println!("huhu, diag: {}\n{:#?}\nscore: {}", diag, last, l);
        self.rec_tr(diag, seq1, seq2, states, last.clone(), l)
    }

    fn rec_tr(
        &self,
        diag: i32,
        seq1: &[u8],
        seq2: &[u8],
        mut current: Vec<Alignment>,
        next_e: WaveFrontElement,
        current_score: usize,
    ) -> Vec<Alignment> {
        if diag == 0 && next_e.offset == 0 {
            println!("ret");
            return current;
        }

        for next_score_d in [
            SCHEME.mismatch as usize,
            SCHEME.gap_extension as usize,
            (SCHEME.gap_opening + SCHEME.gap_extension) as usize,
        ] {
            if next_score_d > current_score {
                println!("well shit");
                continue;
            }
            let next_score = current_score - next_score_d;
            println!("yeah, score: {}", next_score);
            match self {
                Self::Global(wf_tensors) => {
                    if let Some(tensor) = wf_tensors.get(next_score) {
                        // potential parents of next_e in tensor. need to find parent point. return only one alignment for now.
                        // TODO: find all alignmwnts, instead of returning after first one.
                        let m = SCHEME.mismatch as usize;
                        let e = SCHEME.gap_extension as usize;
                        if next_score_d == m {
                            if next_e.state != State::M && next_e.parents.contains(&State::M) {
                                if let Some(wf) = tensor
                                    .as_ref()
                                    .and_then(|t| t.m.as_ref().and_then(|m| m.get_element(diag)))
                                {
                                    println!("mismatch");
                                    current[0]
                                        .seq1
                                        .extend(seq1[wf.y(diag)..next_e.y(diag)].iter().rev());
                                    current[0]
                                        .seq2
                                        .extend(seq2[wf.x(diag)..next_e.x(diag)].iter().rev());
                                    return self.rec_tr(
                                        diag,
                                        seq1,
                                        seq2,
                                        current,
                                        wf.clone(),
                                        next_score,
                                    );
                                }
                            }
                        } else if next_score_d == e {
                            if next_e.parents.contains(&State::D) {
                                println!("extend");
                                if let Some(wf) = tensor.as_ref().and_then(|t| {
                                    t.d.as_ref().and_then(|d| d.get_element(diag - 1))
                                }) {
                                    current[0]
                                        .seq1
                                        .extend(seq1[wf.y(diag)..next_e.y(diag)].iter().rev());
                                    current[0].seq2.push(b'-');
                                    current[0]
                                        .seq2
                                        .extend(seq2[wf.x(diag)..next_e.x(diag)].iter().rev());
                                    return self.rec_tr(
                                        diag - 1,
                                        seq1,
                                        seq2,
                                        current,
                                        wf.clone(),
                                        next_score,
                                    );
                                }
                            }
                            if let Some(wf) = tensor
                                .as_ref()
                                .and_then(|t| t.i.as_ref().and_then(|d| d.get_element(diag + 1)))
                            {
                                println!("extend");
                                current[0].seq1.push(b'-');
                                current[0]
                                    .seq1
                                    .extend(seq1[wf.y(diag)..next_e.y(diag)].iter().rev());

                                current[0]
                                    .seq2
                                    .extend(seq2[wf.x(diag)..next_e.x(diag)].iter().rev());
                                return self.rec_tr(
                                    diag + 1,
                                    seq1,
                                    seq2,
                                    current,
                                    wf.clone(),
                                    next_score,
                                );
                            }
                        } else if next_e.parents.contains(&State::M) {
                            println!("open");
                            match next_e.state {
                                State::D => {
                                    if let Some(wf) = tensor.as_ref().and_then(|t| {
                                        t.d.as_ref().and_then(|d| d.get_element(diag - 1))
                                    }) {
                                        current[0]
                                            .seq1
                                            .extend(seq1[wf.y(diag)..next_e.y(diag)].iter().rev());

                                        current[0].seq2.push(b'-');
                                        current[0]
                                            .seq2
                                            .extend(seq2[wf.x(diag)..next_e.x(diag)].iter().rev());
                                        return self.rec_tr(
                                            diag - 1,
                                            seq1,
                                            seq2,
                                            current,
                                            wf.clone(),
                                            next_score,
                                        );
                                    }
                                }
                                State::I => {
                                    if let Some(wf) = tensor.as_ref().and_then(|t| {
                                        t.i.as_ref().and_then(|d| d.get_element(diag + 1))
                                    }) {
                                        current[0].seq1.push(b'-');
                                        current[0]
                                            .seq1
                                            .extend(seq1[wf.y(diag)..next_e.y(diag)].iter().rev());

                                        current[0]
                                            .seq2
                                            .extend(seq2[wf.x(diag)..next_e.x(diag)].iter().rev());
                                        return self.rec_tr(
                                            diag + 1,
                                            seq1,
                                            seq2,
                                            current,
                                            wf.clone(),
                                            next_score,
                                        );
                                    }
                                }
                                State::M => {
                                    if let Some(wf) = tensor.as_ref().and_then(|t| {
                                        t.i.as_ref().and_then(|d| d.get_element(diag + 1))
                                    }) {
                                        current[0].seq1.push(b'-');
                                        current[0]
                                            .seq1
                                            .extend(seq1[wf.y(diag)..next_e.y(diag)].iter().rev());

                                        current[0]
                                            .seq2
                                            .extend(seq2[wf.x(diag)..next_e.x(diag)].iter().rev());
                                        return self.rec_tr(
                                            diag + 1,
                                            seq1,
                                            seq2,
                                            current,
                                            wf.clone(),
                                            next_score,
                                        );
                                    }
                                    if let Some(wf) = tensor.as_ref().and_then(|t| {
                                        t.d.as_ref().and_then(|d| d.get_element(diag - 1))
                                    }) {
                                        current[0]
                                            .seq1
                                            .extend(seq1[wf.y(diag)..next_e.y(diag)].iter().rev());

                                        current[0].seq1.push(b'-');
                                        current[0]
                                            .seq2
                                            .extend(seq2[wf.x(diag)..next_e.x(diag)].iter().rev());
                                        return self.rec_tr(
                                            diag - 1,
                                            seq1,
                                            seq2,
                                            current,
                                            wf.clone(),
                                            next_score,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
                Self::Local(_) => {}
                Self::SemiGlobal(_) => {}
            }
        }
        println!("huh");
        current
    }

    fn rec_tr_2(
        &self,
        seq1: &[u8],
        seq2: &[u8],
        current_diag: i32,
        current_score: usize,
        mut current_element: WaveFrontElement,
        alignments: Vec<Alignment>,
    ) -> Vec<Alignment> {
        if current_diag == 0 && current_element.offset == 0 {
            return alignments;
        }
        match current_element.state {
            State::M => {
                for parent in current_element.parents {
                    match parent {
                        State::M => {
                            //msimatch
                            // diag = diag, offset = offset -1
                            if current_score < SCHEME.mismatch as usize {
                                return Vec::new();
                            }
                            let wfs = if let Self::Global(wfs) = &self {
                                wfs
                            } else {
                                return Vec::new();
                            };
                            if let Some(prev_element) = wfs
                                .get(current_score - SCHEME.mismatch as usize)
                                .and_then(|wf| {
                                    wf.as_ref().and_then(|wf| {
                                        wf.m.as_ref().and_then(|m| m.get_element(current_diag))
                                    })
                                })
                            {
                                let d_offset = prev_element.offset.abs_diff(current_element.offset);
                                // generate all possible alignemnts
                            }
                        }
                        State::D => {
                            // could be any, delegates to D, I
                        }
                        State::I => {
                            // could be any, delegates to D, I
                        }
                    }
                }
            }
            State::D => {
                for parent in current_element.parents {
                    match parent {
                        State::M => {
                            // gap open
                            // diag = diag - 1, offset = offset
                            let prev_diag = current_diag - 1;
                        }
                        State::D => {
                            // gap extend
                            // diag = diag - 1, offset = offset
                            let prev_diag = current_diag - 1;
                        }
                        State::I => {
                            // impossible
                        }
                    }
                }
            }
            State::I => {
                for parent in current_element.parents {
                    match parent {
                        State::M => {
                            // gap open
                            // diag = diag + 1, offset = offset - 1
                        }
                        State::D => {
                            // impossible
                        }
                        State::I => {
                            // gap extend
                            // diag = diag + 1, offset = offset - 1
                        }
                    }
                }
            }
        }
        alignments
    }
}

#[derive(Debug)]
struct Alignment {
    seq1: Vec<u8>,
    seq2: Vec<u8>,
}

impl Display for Alignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{}",
            self.seq1
                .iter()
                .rev()
                .map(|item| char::from_u32(*item as u32).unwrap())
                .collect::<String>()
        )?;
        for (s1, s2) in self.seq1.iter().zip(self.seq2.iter()).rev() {
            if s1 != s2 {
                write!(f, " ")?;
            } else {
                write!(f, "|")?;
            }
        }

        writeln!(
            f,
            "{}",
            self.seq2
                .iter()
                .rev()
                .map(|item| char::from_u32(*item as u32).unwrap())
                .collect::<String>()
        )?;
        Ok(())
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
                        _marker: PhantomData
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
                    _marker: PhantomData,
                })],
            }),
            m: Some(WaveFront {
                lo: -2,
                hi: 3,
                elements: vec![
                    Some(WaveFrontElement {
                        offset: 1,
                        parents: Vec::default(),
                        state: State::I,
                        _marker: PhantomData
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
                        state: State::I,
                        _marker: PhantomData
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
                        _marker: PhantomData
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
                    _marker: PhantomData,
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
                    _marker: PhantomData,
                })],
            }),
        };

        let true_res_o = WaveFrontTensor {
            i: Some(WaveFront {
                hi: 1,
                lo: 1,
                elements: vec![Some(WaveFrontElement {
                    offset: 1,
                    parents: vec![State::M],
                    state: State::I,
                    _marker: PhantomData,
                })],
            }),
            d: Some(WaveFront {
                hi: -1,
                lo: -1,
                elements: vec![Some(WaveFrontElement {
                    offset: 0,
                    parents: vec![State::M],
                    state: State::D,
                    _marker: PhantomData,
                })],
            }),
            m: Some(WaveFront {
                hi: 1,
                lo: -1,
                elements: vec![
                    Some(WaveFrontElement {
                        offset: 0,
                        parents: vec![State::D],
                        state: State::M,
                        _marker: PhantomData,
                    }),
                    None,
                    Some(WaveFrontElement {
                        offset: 1,
                        parents: vec![State::I],
                        state: State::M,
                        _marker: PhantomData,
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
                    parents: vec![State::M],
                    state: State::M,
                    _marker: PhantomData,
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
        assert!(initial.is_converged(query, db).is_none())
    }
}
