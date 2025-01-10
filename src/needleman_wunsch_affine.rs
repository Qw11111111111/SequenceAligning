use std::{
    cmp::{Eq, Ordering, PartialEq},
    collections::VecDeque,
    fmt::{Debug, Display, Write},
    rc::Rc,
    time::Instant,
};

use crate::errors::{AlignerError, Result};
use crate::parse::{Args, Mode, Record};
use crate::utils::vec_u8_to_str;

type Array<'a> = Vec<Vec<Rc<ArrayElement<'a>>>>;

const SCHEME: ScoringScheme = ScoringScheme {
    gap_opening: -8,
    gap_extension: -6,
    mismatch: -4,
    match_: 5,
};

#[derive(Default, Clone, Eq, PartialEq, Debug)]
struct ArrayElement<'a> {
    score: i32,
    parents: Vec<Rc<ArrayElement<'a>>>,
    state: State,
    _marker: std::marker::PhantomData<&'a u8>,
}

impl<'a> ArrayElement<'a> {
    fn new(score: i32, parents: Vec<Rc<ArrayElement<'a>>>, state: State) -> Self {
        Self {
            score,
            parents,
            state,
            _marker: std::marker::PhantomData,
        }
    }
}

impl<'a> Ord for ArrayElement<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score.cmp(&other.score)
    }
}

impl<'a> PartialOrd for ArrayElement<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

enum ScoreTensor<'a> {
    Global(ScoreMatrix<'a>),
    SemiGlobal(ScoreMatrix<'a>),
    Local(ScoreMatrix<'a>),
}

#[derive(Default)]
struct ScoreMatrix<'a> {
    m_scores: Array<'a>,
    i_scores: Array<'a>,
    d_scores: Array<'a>,
}

impl<'a> ScoreMatrix<'a> {
    fn new(x: usize, y: usize) -> Self {
        // seq1 == y seq2 == x
        Self {
            m_scores: default_array(x, y),
            i_scores: default_array(x, y), // insertion in seq1 aka gap in seq2
            d_scores: default_array(x, y), // deletion in seq1 aka gap in seq1
        }
    }

    fn m_score(&self, x: usize, y: usize, is_match: bool) -> i32 {
        self.m_scores[x - 1][y - 1]
            .score
            .max(self.i_scores[x - 1][y - 1].score)
            .max(self.d_scores[x - 1][y - 1].score)
            + if is_match {
                SCHEME.match_
            } else {
                SCHEME.mismatch
            }
    }
    fn d_score(&self, x: usize, y: usize) -> i32 {
        (self.m_scores[x - 1][y].score + SCHEME.gap_opening).max(self.d_scores[x - 1][y].score)
            + SCHEME.gap_extension
    }
    fn i_score(&self, x: usize, y: usize) -> i32 {
        (self.m_scores[x][y - 1].score + SCHEME.gap_opening).max(self.i_scores[x][y - 1].score)
            + SCHEME.gap_extension
    }

    fn d_pointer(&self, x: usize, y: usize) -> Vec<Rc<ArrayElement<'a>>> {
        let mut p = vec![];
        if self.d_score(x, y) == self.d_scores[x - 1][y].score + SCHEME.gap_extension {
            p.push(Rc::clone(&self.d_scores[x - 1][y]));
        }
        if self.d_score(x, y)
            == self.m_scores[x - 1][y].score + SCHEME.gap_opening + SCHEME.gap_extension
        {
            p.push(Rc::clone(&self.m_scores[x - 1][y]));
        }
        p
    }
    fn i_pointer(&self, x: usize, y: usize) -> Vec<Rc<ArrayElement<'a>>> {
        let mut p = Vec::default();
        if self.i_score(x, y) == self.i_scores[x][y - 1].score + SCHEME.gap_extension {
            p.push(Rc::clone(&self.i_scores[x][y - 1]));
        }
        if self.i_score(x, y)
            == self.m_scores[x][y - 1].score + SCHEME.gap_opening + SCHEME.gap_extension
        {
            p.push(Rc::clone(&self.m_scores[x][y - 1]));
        }
        p
    }
    fn m_pointer(&self, x: usize, y: usize, is_match: bool) -> Vec<Rc<ArrayElement<'a>>> {
        let mut p = Vec::default();
        if self.m_score(x, y, is_match)
            == self.m_scores[x - 1][y - 1].score
                + if is_match {
                    SCHEME.match_
                } else {
                    SCHEME.mismatch
                }
        {
            p.push(Rc::clone(&self.m_scores[x - 1][y - 1]));
        }
        if self.m_score(x, y, is_match)
            == self.i_scores[x - 1][y - 1].score
                + if is_match {
                    SCHEME.match_
                } else {
                    SCHEME.mismatch
                }
        {
            p.push(Rc::clone(&self.i_scores[x - 1][y - 1]));
        }
        if self.m_score(x, y, is_match)
            == self.d_scores[x - 1][y - 1].score
                + if is_match {
                    SCHEME.match_
                } else {
                    SCHEME.mismatch
                }
        {
            p.push(Rc::clone(&self.d_scores[x - 1][y - 1]));
        }
        p
    }
}

impl<'a> ScoreTensor<'a> {
    fn global(x: usize, y: usize) -> Self {
        ScoreTensor::Global(ScoreMatrix::new(x, y))
    }

    fn local(x: usize, y: usize) -> Self {
        ScoreTensor::Local(ScoreMatrix::new(x, y))
    }

    fn semi_global(x: usize, y: usize) -> Self {
        ScoreTensor::SemiGlobal(ScoreMatrix::new(x, y))
    }

    fn fill(&mut self, seq1: &[u8], seq2: &[u8]) {
        match self {
            ScoreTensor::Global(matrix) => {
                matrix.m_scores[0][0] = Rc::new(ArrayElement::new(0, Vec::default(), State::InM));
                matrix.d_scores[0][0] = Rc::new(ArrayElement::new(
                    i16::MIN as i32,
                    Vec::default(),
                    State::InD,
                ));
                matrix.i_scores[0][0] = Rc::new(ArrayElement::new(
                    i16::MIN as i32,
                    Vec::default(),
                    State::InI,
                ));
                (1..=seq1.len()).for_each(|i| {
                    matrix.m_scores[0][i] = Rc::new(ArrayElement::new(
                        i16::MIN as i32,
                        Vec::default(),
                        State::InM,
                    ));
                    matrix.i_scores[0][i] = Rc::new(ArrayElement::new(
                        i16::MIN as i32,
                        Vec::default(),
                        State::InI,
                    ));
                    matrix.d_scores[0][i] = Rc::new(ArrayElement::new(
                        (i as i32 + 1) * SCHEME.gap_extension + SCHEME.gap_opening,
                        vec![Rc::clone(&matrix.d_scores[0][i - 1])],
                        State::InD,
                    ));
                });
                (1..=seq2.len()).for_each(|i| {
                    matrix.m_scores[i][0] = Rc::new(ArrayElement::new(
                        i16::MIN as i32,
                        Vec::default(),
                        State::InM,
                    ));
                    matrix.i_scores[i][0] = Rc::new(ArrayElement::new(
                        SCHEME.gap_opening + (i as i32 + 1) * SCHEME.gap_extension,
                        vec![Rc::clone(&matrix.i_scores[i - 1][0])],
                        State::InI,
                    ));
                    matrix.d_scores[i][0] = Rc::new(ArrayElement::new(
                        i16::MIN as i32,
                        Vec::default(),
                        State::InD,
                    ));
                });
                for i in 1..=seq2.len() {
                    for j in 1..=seq1.len() {
                        matrix.m_scores[i][j] = Rc::new(ArrayElement::new(
                            matrix.m_score(i, j, seq1[j - 1] == seq2[i - 1]),
                            matrix.m_pointer(i, j, seq1[j - 1] == seq2[i - 1]),
                            State::InM,
                        ));

                        matrix.i_scores[i][j] = Rc::new(ArrayElement::new(
                            matrix.i_score(i, j),
                            matrix.i_pointer(i, j),
                            State::InI,
                        ));
                        matrix.d_scores[i][j] = Rc::new(ArrayElement::new(
                            matrix.d_score(i, j),
                            matrix.d_pointer(i, j),
                            State::InD,
                        ));
                    }
                }
            }
            ScoreTensor::Local(_matrix) => { /* local fill logic */ }
            ScoreTensor::SemiGlobal(_matrix) => { /* semi-global fill logic */ }
        }
    }
    fn traceback(&self, seq1: &[u8], seq2: &[u8]) {
        match self {
            ScoreTensor::Global(matrix) => {
                //println!("{:#?}", matrix);
                let mut queue: Vec<TraceBackInfo> = Vec::default();
                let max_val = matrix.i_scores[seq2.len()][seq1.len()]
                    .score
                    .max(matrix.d_scores[seq2.len()][seq1.len()].score)
                    .max(matrix.m_scores[seq2.len()][seq1.len()].score);
                if max_val == matrix.i_scores[seq2.len()][seq1.len()].score {
                    queue.push(TraceBackInfo {
                        seq1: Vec::default(),
                        seq2: Vec::default(),
                        current_cell: Rc::clone(&matrix.i_scores[seq2.len()][seq1.len()]),
                        current_state: State::InI,
                        x: seq2.len(),
                        y: seq1.len(),
                    });
                }
                if max_val == matrix.m_scores[seq2.len()][seq1.len()].score {
                    queue.push(TraceBackInfo {
                        seq1: Vec::default(),
                        seq2: Vec::default(),
                        current_cell: Rc::clone(&matrix.m_scores[seq2.len()][seq1.len()]),
                        current_state: State::InM,
                        x: seq2.len(),
                        y: seq1.len(),
                    });
                }
                if max_val == matrix.d_scores[seq2.len()][seq1.len()].score {
                    queue.push(TraceBackInfo {
                        seq1: Vec::default(),
                        seq2: Vec::default(),
                        current_cell: Rc::clone(&matrix.d_scores[seq2.len()][seq1.len()]),
                        current_state: State::InD,
                        x: seq2.len(),
                        y: seq1.len(),
                    });
                }
                while let Some(element) = queue.pop() {
                    //println!("e: {:#?}", element);
                    if element.x == 0 && element.y == 0 {
                        println!("alignment found");
                        println!("{}", element);
                    }
                    for parent in &element.current_cell.parents {
                        let mut s1: Vec<u8>;
                        let mut s2: Vec<u8>;
                        let mut x = element.x;
                        let mut y = element.y;
                        if element.current_state == State::InM {
                            s1 = vec![seq1[element.y - 1]];
                            s2 = vec![seq2[element.x - 1]];
                            s1.extend_from_slice(&element.seq1);
                            s2.extend_from_slice(&element.seq2);
                        } else if element.current_state == State::InD {
                            s1 = vec![b'-'];
                            s2 = vec![seq2[element.x - 1]];
                            s1.extend_from_slice(&element.seq1);
                            s2.extend_from_slice(&element.seq2);
                        } else {
                            s1 = vec![seq1[element.y - 1]];
                            s2 = vec![b'-'];
                            s1.extend_from_slice(&element.seq1);
                            s2.extend_from_slice(&element.seq2);
                        }
                        match element.current_state {
                            State::InM => {
                                x -= 1;
                                y -= 1;
                            }
                            State::InD => {
                                x -= 1;
                            }
                            State::InI => {
                                y -= 1;
                            }
                        }
                        queue.push(TraceBackInfo {
                            seq1: s1,
                            seq2: s2,
                            current_cell: Rc::clone(parent),
                            current_state: parent.state.clone(),
                            x,
                            y,
                        });
                    }
                }
            }
            ScoreTensor::SemiGlobal(_) => {}
            ScoreTensor::Local(_) => {}
        }
    }
}

impl<'a> Debug for ScoreMatrix<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "m: ")?;
        for i in 0..self.m_scores.len() {
            for j in 0..self.m_scores[0].len() {
                write!(f, " {} ", self.m_scores[i][j].score)?;
            }
            writeln!(f)?;
        }
        writeln!(f, "d: ")?;
        for i in 0..self.d_scores.len() {
            for j in 0..self.d_scores[0].len() {
                write!(f, " {} ", self.d_scores[i][j].score)?;
            }
            writeln!(f)?;
        }
        writeln!(f, "i: ")?;
        for i in 0..self.i_scores.len() {
            for j in 0..self.i_scores[0].len() {
                write!(f, " {} ", self.i_scores[i][j].score)?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
enum State {
    #[default]
    InM,
    InD,
    InI,
}

struct TraceBackInfo<'a> {
    seq1: Vec<u8>,
    seq2: Vec<u8>,
    current_cell: Rc<ArrayElement<'a>>,
    current_state: State,
    x: usize,
    y: usize,
}

#[derive(Default, Debug)]
struct ScoringScheme {
    gap_opening: i32,
    gap_extension: i32,
    mismatch: i32,
    match_: i32,
}

impl<'a> Display for TraceBackInfo<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\nseq1: {}", vec_u8_to_str(&self.seq1))?;
        write!(
            f,
            "\n      {}",
            self.seq1
                .iter()
                .zip(self.seq2.iter())
                .map(|(c1, c2)| {
                    if c1 == c2 {
                        '|'
                    } else {
                        ' '
                    }
                })
                .collect::<String>()
        )?;
        write!(f, "\nseq2: {}", vec_u8_to_str(&self.seq2))?;
        Ok(())
    }
}

impl<'a> Debug for TraceBackInfo<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{:#?}", self.current_cell)?;
        Ok(())
    }
}

fn default_array<'a>(x: usize, y: usize) -> Array<'a> {
    vec![vec![Rc::default(); y + 1]; x + 1]
}

pub fn n_w_align<'a>(seq1: &Record, seq2: &Record, _verbose: bool, mode: Mode) -> Result<'a, ()> {
    let now = Instant::now();
    match mode {
        Mode::Global => {
            let mut mat = ScoreTensor::global(seq2.seq.len(), seq1.seq.len());
            mat.fill(&seq1.seq, &seq2.seq);
            mat.traceback(&seq1.seq, &seq2.seq);
            println!("{:#?}", now.elapsed());
        }
        Mode::SemiGlobal => return Err(AlignerError::AlignmentError("not implemented")),
        Mode::Local => return Err(AlignerError::AlignmentError("not implemented")),
    }
    Ok(())
}
/*
fn argmax(mat: &Array) -> Vec<(usize, usize)> {
    let mut idc = Vec::default();
    let mut max = i32::MIN;
    for i in 0..mat.len() {
        for j in 0..mat[0].len() {
            match mat[i][j].cmp(&max) {
                Ordering::Greater => {
                    max = mat[i][j];
                    idc = vec![(i, j)];
                }
                Ordering::Equal => idc.push((i, j)),
                _ => (),
            }
        }
    }
    idc
}
*/

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_mat_fill() {
        //todo!();
    }

    #[test]
    fn test_traceback() {
        //todo!()
    }
}
