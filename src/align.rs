use crate::errors::{AStarError, Result};
use crate::parse::Record;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::rc::Rc;

const SCHEME: ScoringScheme = ScoringScheme {
    match_score: 0,
    mismatch_score: 4,
    gap_opening: 2,
    gap_extension: 8,
    epsilon: 1.5,
};

pub fn align(seq1: &Record, seq2: &Record) -> Result<()> {
    let target_length = seq1.seq.len().max(seq2.seq.len());
    let mut queue = BinaryHeap::new();
    queue.push(State {
        reach_cost: 0,
        cost: get_h(&seq1.seq, &seq2.seq, 0, 0, target_length),
        position: Position { x: 0, y: 0 },
        parent: None,
    });
    while let Some(s) = queue.pop() {
        //TODO calculate cost of node as heuristic + cost_to reach
        // also need to somehow save the cumulative cost up to that point. and the path to the point
        if is_converged(&s, &seq1.seq, &seq2.seq) {
            println!(
                "Alignment for db {} and query {} found",
                seq2.name.iter().map(|&i| i as char).collect::<String>(),
                seq1.name.iter().map(|&i| i as char).collect::<String>()
            );
            pprint(s, &seq1.seq, &seq2.seq);
            break;
        }
        let p = Rc::from(s);
        let pos = &p.position;
        if pos.x < seq2.seq.len() {
            queue.push(State {
                cost: get_h(&seq1.seq, &seq2.seq, pos.x, pos.y, target_length),
                reach_cost: p.reach_cost + SCHEME.gap_extension,
                position: Position {
                    x: pos.x + 1,
                    y: pos.y,
                },
                parent: Some(p.clone()),
            });
        }
        if pos.y < seq1.seq.len() {
            queue.push(State {
                cost: get_h(&seq1.seq, &seq2.seq, pos.x, pos.y, target_length),
                reach_cost: p.reach_cost + SCHEME.gap_extension,
                position: Position {
                    x: pos.x,
                    y: pos.y + 1,
                },
                parent: Some(p.clone()),
            })
        }
        if pos.y < seq1.seq.len() && pos.x < seq2.seq.len() {
            queue.push(State {
                cost: get_h(&seq1.seq, &seq2.seq, pos.x, pos.y, target_length),
                reach_cost: p.reach_cost + get_cost(&seq1.seq[pos.y], &seq2.seq[pos.x]),
                position: Position {
                    x: pos.x + 1,
                    y: pos.y + 1,
                },
                parent: Some(p.clone()),
            })
        }
    }
    Ok(())
}

struct ScoringScheme {
    match_score: usize,
    mismatch_score: usize,
    gap_opening: usize,
    gap_extension: usize,
    epsilon: f64,
}

fn get_h(seq1: &[u8], seq2: &[u8], x: usize, y: usize, target_length: usize) -> usize {
    ((1. + SCHEME.epsilon * dynamic_weight(x, y, target_length)) * heuristic_d(seq1, seq2, x, y))
        as usize
}

fn dynamic_weight(x: usize, y: usize, target_length: usize) -> f64 {
    if x.max(y) <= target_length {
        1. - x.max(y) as f64 / target_length as f64
    } else {
        0.
    }
}

fn heuristic_d(seq1: &[u8], seq2: &[u8], x: usize, y: usize) -> f64 {
    // ((seq1.len() - y).abs_diff(seq2.len() - x)) * SCHEME.gap_extension
    ((seq1.len() - y) + (seq2.len() - x)) as f64
}

fn _heuristic_d(seq1: &[u8], seq2: &[u8], x: usize, y: usize) -> f64 {
    avg_step_cost(seq1, seq2) * (seq1.len() - y + seq2.len() - x) as f64
}

fn avg_step_cost(seq1: &[u8], seq2: &[u8]) -> f64 {
    (((seq1.len().abs_diff(seq2.len()) * SCHEME.gap_extension) as f64
        + seq1.len().min(seq2.len()) as f64
            * (0.75 * SCHEME.mismatch_score as f64 + 0.25 * SCHEME.match_score as f64))
        / seq1.len().max(seq2.len()) as f64)
}

fn is_converged(current: &State, seq1: &[u8], seq2: &[u8]) -> bool {
    current.position.x == seq2.len() && current.position.y == seq1.len()
}

fn pprint(state: State, seq1: &[u8], seq2: &[u8]) {
    let mut db = String::with_capacity(seq1.len().max(seq2.len()));
    let mut q = String::with_capacity(seq1.len().max(seq2.len()));
    let mut x = state.position.x;
    let mut y = state.position.y;
    let mut current = state.parent;
    while let Some(parent) = current {
        if parent.position.x == x {
            y -= 1;
            db.push('-');
            q.push(seq1[y] as char);
        } else if parent.position.y == y {
            x -= 1;
            db.push(seq2[x] as char);
            q.push('-');
        } else {
            x -= 1;
            y -= 1;
            db.push(seq2[x] as char);
            q.push(seq1[y] as char);
        }
        current = parent.parent.clone();
    }
    println!("{}", db);
    println!(
        "{}",
        q.chars()
            .zip(db.chars())
            .map(|(q, d)| if q == d { '|' } else { ' ' })
            .collect::<String>()
    );
    println!("{}", q);
}

#[derive(Debug, Default, PartialEq, Eq, Clone)]
struct State {
    cost: usize,
    reach_cost: usize,
    position: Position,
    parent: Option<Rc<State>>,
}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        (other.cost + other.reach_cost)
            .cmp(&(self.cost + self.reach_cost))
            .then_with(|| self.position.cmp(&other.position))
            .then_with(|| self.parent.cmp(&other.parent))
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Default, Clone)]
struct Position {
    x: usize,
    y: usize,
}

fn get_cost(c1: &u8, c2: &u8) -> usize {
    if *c1 == *c2 || *c1 == b'N' || *c2 == b'N' {
        SCHEME.match_score
    } else {
        SCHEME.mismatch_score
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_heuristic() {
        let seq1 = b"AATG";
        let seq2 = b"AATGAA";
        assert!(
            heuristic_d(seq1, seq2, 0, 0) as usize
                <= 2 * SCHEME.gap_extension + 4 * SCHEME.match_score,
            "expected: {}, actual: {}",
            2 * SCHEME.gap_extension + 4 * SCHEME.match_score,
            heuristic_d(seq1, seq2, 0, 0)
        );
    }
    #[test]
    fn test_queue() {
        let mut q = BinaryHeap::new();
        q.push(State {
            cost: 10,
            reach_cost: 0,
            position: Position { x: 0, y: 0 },
            parent: None,
        });
        q.push(State {
            cost: 5,
            reach_cost: 4,
            position: Position { x: 2, y: 3 },
            parent: None,
        });
        assert_eq!(
            q.peek(),
            Some(&State {
                cost: 5,
                reach_cost: 4,
                position: Position { x: 2, y: 3 },
                parent: None,
            })
        )
    }
}
