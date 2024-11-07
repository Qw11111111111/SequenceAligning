use std::{
    cmp::{Eq, Ordering, PartialEq},
    collections::VecDeque,
    fmt::{Display, Write},
    time::Instant,
};

use crate::errors::{AStarError, Result};
use crate::parse::{Args, Mode, Record};
use crate::utils::vec_u8_to_str;

type Array = Vec<Vec<i32>>;

const SCHEME: ScoringScheme = ScoringScheme {
    gap_opening: -8,
    gap_extension: -6,
    mismatch: -4,
    match_: 5,
};

enum ScoreTensor {
    Global(ScoreMatrix),
    SemiGlobal(ScoreMatrix),
    Local(ScoreMatrix),
}

#[derive(Default, Debug)]
struct ScoreMatrix {
    m_scores: Array,
    i_scores: Array,
    d_scores: Array,
}

impl ScoreMatrix {
    fn new(x: usize, y: usize) -> Self {
        // seq1 == y seq2 == x
        Self {
            m_scores: default_array(x, y),
            i_scores: default_array(x, y), // insertion in seq1 aka gap in seq2
            d_scores: default_array(x, y), // deletion in seq1 aka gap in seq1
        }
    }

    fn m_score(&self, x: usize, y: usize, char1: &u8, char2: &u8) -> i32 {
        self.m_scores[x - 1][y - 1]
            .max(self.i_scores[x - 1][y - 1])
            .max(self.d_scores[x - 1][y - 1])
            + if *char1 == *char2 {
                SCHEME.match_
            } else {
                SCHEME.mismatch
            }
    }
    fn d_score(&self, x: usize, y: usize) -> i32 {
        (self.m_scores[x - 1][y] + SCHEME.gap_opening).max(self.d_scores[x - 1][y])
            + SCHEME.gap_extension
    }
    fn i_score(&self, x: usize, y: usize) -> i32 {
        (self.m_scores[x][y - 1] + SCHEME.gap_opening).max(self.i_scores[x][y - 1])
            + SCHEME.gap_extension
    }
}

impl ScoreTensor {
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
                (0..seq1.len()).for_each(|i| {
                    matrix.m_scores[0][i] = i32::MIN;
                    matrix.i_scores[0][i] = i32::MIN;
                    matrix.d_scores[0][i] =
                        (i as i32 + 1) * SCHEME.gap_extension + SCHEME.gap_opening;
                });
                (0..seq2.len()).for_each(|i| {
                    matrix.m_scores[i][0] = i32::MIN;
                    matrix.i_scores[i][0] =
                        SCHEME.gap_opening + (i as i32 + 1) * SCHEME.gap_extension;
                    matrix.d_scores[i][0] = i32::MIN;
                });
                (1..=seq2.len()).for_each(|i| {
                    (1..=seq1.len()).for_each(|j| {
                        matrix.m_scores[i][j] = matrix.m_score(i, j, &seq1[j - 1], &seq2[i - 1]);
                        matrix.i_scores[i][j] = matrix.i_score(i, j);
                        matrix.d_scores[i][j] = matrix.d_score(i, j);
                    })
                });
            }
            ScoreTensor::Local(matrix) => { /* local fill logic */ }
            ScoreTensor::SemiGlobal(matrix) => { /* semi-global fill logic */ }
        }
    }

    fn traceback(&self, seq1: &[u8], seq2: &[u8]) {
        match self {
            ScoreTensor::Global(matrix) => {
                let mut queue = Vec::default();
                let max = matrix.m_scores[seq2.len()][seq1.len()]
                    .max(matrix.i_scores[seq2.len()][seq1.len()])
                    .max(matrix.d_scores[seq2.len()][seq1.len()]);
                if matrix.m_scores[seq2.len()][seq1.len()] == max {
                    queue.push(TraceBackInfo {
                        x: seq2.len(),
                        y: seq1.len(),
                        seq1: VecDeque::default(),
                        seq2: VecDeque::default(),
                        state: State::InM,
                    });
                }
                if matrix.d_scores[seq2.len()][seq1.len()] == max {
                    queue.push(TraceBackInfo {
                        x: seq2.len(),
                        y: seq1.len(),
                        seq1: VecDeque::default(),
                        seq2: VecDeque::default(),
                        state: State::InD,
                    });
                }
                if matrix.i_scores[seq2.len()][seq1.len()] == max {
                    queue.push(TraceBackInfo {
                        x: seq2.len(),
                        y: seq1.len(),
                        seq1: VecDeque::default(),
                        seq2: VecDeque::default(),
                        state: State::InI,
                    });
                }
                while traceback_step(
                    &mut queue,
                    &matrix.m_scores,
                    &matrix.i_scores,
                    &matrix.d_scores,
                    seq1,
                    seq2,
                )
                .is_ok()
                {}
            }
            ScoreTensor::Local(_) => {}
            _ => {}
        }
    }
}
#[derive(Eq, PartialEq, Debug)]
enum State {
    InM,
    InD,
    InI,
}

#[derive(Debug)]
struct TraceBackInfo {
    x: usize,
    y: usize,
    seq1: VecDeque<u8>,
    seq2: VecDeque<u8>,
    state: State,
}

fn traceback_step<'a>(
    queue: &mut Vec<TraceBackInfo>,
    m_scores: &Array,
    i_scores: &Array,
    d_scores: &Array,
    seq1: &[u8],
    seq2: &[u8],
) -> Result<'a, ()> {
    if let Some(info) = queue.pop() {
        if info.x == 0 && info.y == 0 {
            println!("{}", info);
            return Ok(());
        }
        println!("{:?}", info);
        match info.state {
            State::InM => {
                if info.x == 0 || info.y == 0 {
                    return Ok(());
                }
                push_next(
                    queue,
                    m_scores,
                    m_scores,
                    &info,
                    seq1,
                    seq2,
                    State::InM,
                    1,
                    1,
                    if seq1[info.y - 1] == seq2[info.x - 1] {
                        SCHEME.match_
                    } else {
                        SCHEME.mismatch
                    },
                );
                push_next(
                    queue,
                    m_scores,
                    d_scores,
                    &info,
                    seq1,
                    seq2,
                    State::InD,
                    1,
                    1,
                    if seq1[info.y - 1] == seq2[info.x - 1] {
                        SCHEME.match_
                    } else {
                        SCHEME.mismatch
                    },
                );
                push_next(
                    queue,
                    m_scores,
                    i_scores,
                    &info,
                    seq1,
                    seq2,
                    State::InI,
                    1,
                    1,
                    if seq1[info.y - 1] == seq2[info.x - 1] {
                        SCHEME.match_
                    } else {
                        SCHEME.mismatch
                    },
                );
            }
            State::InD => {
                push_next(
                    queue,
                    d_scores,
                    m_scores,
                    &info,
                    seq1,
                    seq2,
                    State::InM,
                    1,
                    0,
                    SCHEME.gap_opening + SCHEME.gap_extension,
                );
                push_next(
                    queue,
                    d_scores,
                    d_scores,
                    &info,
                    seq1,
                    seq2,
                    State::InD,
                    1,
                    0,
                    SCHEME.gap_extension,
                );
            }
            State::InI => {
                push_next(
                    queue,
                    i_scores,
                    m_scores,
                    &info,
                    seq1,
                    seq2,
                    State::InM,
                    0,
                    1,
                    SCHEME.gap_opening + SCHEME.gap_extension,
                );
                push_next(
                    queue,
                    i_scores,
                    i_scores,
                    &info,
                    seq1,
                    seq2,
                    State::InI,
                    0,
                    1,
                    SCHEME.gap_extension,
                );
            }
        }
        println!("ok");
        return Ok(());
    }
    Err(AStarError::AlignmentError("tried to step on empty queue"))
}

fn push_next(
    queue: &mut Vec<TraceBackInfo>,
    arr: &Array,
    arr2: &Array,
    info: &TraceBackInfo,
    seq1: &[u8],
    seq2: &[u8],
    state: State,
    delta_x: usize,
    delta_y: usize,
    delta_score: i32,
) {
    if delta_x > info.x || delta_y > info.y {
        return;
    }
    if arr2[info.x - delta_x][info.y - delta_y] + delta_score == arr[info.x][info.y] {
        let mut s1 = info.seq1.clone();
        let mut s2 = info.seq2.clone();
        if state != State::InD {
            s1.push_front(seq1[info.y - 1])
        } else {
            s1.push_front(b'-')
        };
        if state != State::InI {
            s2.push_front(seq2[info.x - 1])
        } else {
            s2.push_front(b'-')
        };
        queue.push(TraceBackInfo {
            x: info.x - delta_x,
            y: info.y - delta_y,
            seq1: s1,
            seq2: s2,
            state,
        })
    } else {
        println!("not pushed");
    }
}

impl Display for TraceBackInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\nseq1: {}", vec_u8_to_str(&self.seq1.as_slices().0))?;
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
        write!(f, "\nseq2: {}", vec_u8_to_str(&self.seq2.as_slices().0))?;
        Ok(())
    }
}

#[derive(Default, Debug)]
struct ScoringScheme {
    gap_opening: i32,
    gap_extension: i32,
    mismatch: i32,
    match_: i32,
}

#[derive(Default, Debug)]
struct Hit {
    query: String,
    db: String,
    start_in_query: usize,
    start_in_db: usize,
}

impl Display for Hit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "\nseq1: {}",
            self.query.chars().rev().collect::<String>()
        )?;
        write!(f, "\n      ")?;
        for (s1, s2) in self.query.chars().rev().zip(self.db.chars().rev()) {
            if s1 == s2 {
                write!(f, "|")?;
            } else {
                write!(f, " ")?;
            }
        }
        write!(f, "\nseq2: {}\n", self.db.chars().rev().collect::<String>())?;
        writeln!(
            f,
            "start in seq1: {}\nstart in seq2: {}\n",
            self.start_in_query, self.start_in_db
        )?;
        Ok(())
    }
}

fn default_array(x: usize, y: usize) -> Array {
    vec![vec![0; y + 1]; x + 1]
}

pub fn n_w_align<'a>(seq1: &Record, seq2: &Record, verbose: bool, mode: Mode) -> Result<'a, ()> {
    let now = Instant::now();
    match mode {
        Mode::Global => {
            let mut mat = ScoreTensor::global(seq2.seq.len(), seq1.seq.len());
            mat.fill(&seq1.seq, &seq2.seq);
            mat.traceback(&seq1.seq, &seq2.seq);
            println!("{:#?}", now.elapsed());
        }
        Mode::SemiGlobal => {}
        Mode::Local => {}
    }
    Ok(())
}

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
