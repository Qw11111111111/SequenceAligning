use crate::errors::{AStarError, Result};
use crate::parse::{Record, Records};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

type EdgeArray = Vec<Vec<Vec<Edge>>>;

const SCHEME: ScoringScheme = ScoringScheme {
    match_score: 0,
    mismatch_score: 2,
    gap_opening: 2,
    gap_extension: 4,
};

pub fn align(seq1: &Record, seq2: &Record) -> Result<()> {
    let mut graph = match Graph::from_seq(&seq1.seq, &seq2.seq) {
        Err(AStarError::AlignmentError(_)) => {
            eprintln!(
                "One of the provided sequences in records {:#?} or {:#?} was empty. Skipping",
                seq1.name, seq2.name
            );
            return Ok(());
        }
        Ok(res) => res,
        _ => {
            return Err(AStarError::AlignmentError(
                "Something unexpected happend during grao constructio".into(),
            ))
        }
    };
    let mut queue = BinaryHeap::new();
    queue.push(State {
        reach_cost: 0,
        cost: heuristic_d(&seq1.seq, &seq2.seq, 0, 0),
        position: Position { x: 0, y: 0 },
        moves: Vec::with_capacity(seq1.seq.len().max(seq2.seq.len())),
    });
    while let Some(s) = queue.pop() {
        //graph.extend(&s.position);
        //TODO calculate cost of node as heuristic + cost_to reach
        // also need to somehow save the cumulative cost up to that point. and the path to the point
        if is_converged(&s, &seq1.seq, &seq2.seq) {
            println!(
                "Alignment for db {} and query {} found",
                seq2.name.iter().map(|&i| i as char).collect::<String>(),
                seq1.name.iter().map(|&i| i as char).collect::<String>()
            );
            pprint(&s, &seq1.seq, &seq2.seq);
            break;
        }
        if s.position.x < seq1.seq.len() {
            let mut mvs = s.moves.clone();
            mvs.push(Move::Horizontal);
            queue.push(State {
                cost: heuristic_d(&seq1.seq, &seq2.seq, s.position.x, s.position.y),
                reach_cost: s.reach_cost,
                position: Position {
                    x: s.position.x + 1,
                    y: s.position.y,
                },
                moves: mvs,
            });
        }
        if s.position.y < seq2.seq.len() {
            let mut mvs = s.moves.clone();
            mvs.push(Move::Vertical);
            queue.push(State {
                cost: heuristic_d(&seq1.seq, &seq2.seq, s.position.x, s.position.y),
                reach_cost: s.reach_cost,
                position: Position {
                    x: s.position.x,
                    y: s.position.y + 1,
                },
                moves: mvs,
            })
        }
        if s.position.y < seq2.seq.len() && s.position.x < seq1.seq.len() {
            let mut mvs = s.moves.clone();
            mvs.push(Move::Diagonal);
            queue.push(State {
                cost: heuristic_d(&seq1.seq, &seq2.seq, s.position.x, s.position.y),
                reach_cost: s.reach_cost
                    + if seq1.seq[s.position.y] == seq2.seq[s.position.x] {
                        SCHEME.match_score
                    } else {
                        SCHEME.mismatch_score
                    },
                position: Position {
                    x: s.position.x + 1,
                    y: s.position.y + 1,
                },
                moves: mvs,
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
}

fn heuristic_d(seq1: &[u8], seq2: &[u8], x: usize, y: usize) -> usize {
    0
}

fn is_converged(current: &State, seq1: &[u8], seq2: &[u8]) -> bool {
    current.position.x == seq2.len() && current.position.y == seq1.len()
}

fn pprint(state: &State, seq1: &[u8], seq2: &[u8]) {
    let mut db = String::with_capacity(state.moves.len());
    let mut q = String::with_capacity(state.moves.len());
    let mut x = state.position.x;
    let mut y = state.position.y;
    for mv in state.moves.iter().rev() {
        match mv {
            Move::Diagonal => {
                db.push(seq2[y] as char);
                q.push(seq1[x] as char);
                y -= 1;
                x -= 1;
            }
            Move::Horizontal => {
                db.push(seq2[x] as char);
                q.push('-');
                x -= 1;
            }
            Move::Vertical => {
                db.push('-');
                q.push(seq1[y] as char);
                y -= 1;
            }
        }
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

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
enum Move {
    Diagonal,
    Horizontal,
    Vertical,
}

#[derive(Debug, Default, PartialEq, Eq)]
struct State {
    cost: usize,
    reach_cost: usize,
    position: Position,
    moves: Vec<Move>,
}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        (other.cost + other.reach_cost)
            .cmp(&(self.cost + self.reach_cost))
            .then_with(|| self.position.cmp(&other.position))
            .then_with(|| self.moves.cmp(&other.moves))
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

#[derive(Default, Debug, Clone)]
struct Edge {
    node: Position,
    cost: usize,
}

#[derive(Debug)]
struct Graph<'a> {
    query: &'a [u8],
    db: &'a [u8],
    edges: EdgeArray,
}

impl<'a> Graph<'a> {
    fn from_seq(seq1: &'a [u8], seq2: &'a [u8]) -> Result<Self> {
        // initializes a graph with the first three edges starting at 0 0
        if seq1.is_empty() || seq2.is_empty() {
            return Err(AStarError::AlignmentError(
                "Tried to align with empty sequence".into(),
            ));
        }

        let edges = vec![vec![Vec::default(); seq1.len()]; seq2.len()];
        /*
        edges[0][0] = vec![
            Edge {
                node: Position { x: 1, y: 0 },
                cost: SCHEME.gap_extension + SCHEME.gap_opening,
            },
            Edge {
                node: Position { x: 1, y: 1 },
                cost: get_cost(&seq1[0], &seq2[0]),
            },
            Edge {
                node: Position { x: 0, y: 1 },
                cost: SCHEME.gap_extension + SCHEME.gap_opening,
            },
        ];
        */
        Ok(Self {
            query: seq1,
            db: seq2,
            edges,
        })
    }

    fn extend(&mut self, node: &Position) {
        if !self.edges[node.x][node.y].is_empty()
            || node.y == self.query.len() && node.x == self.db.len()
        {
            return;
        }
        if node.y == self.query.len() {
            self.edges[node.x][node.y] = vec![Edge {
                node: Position {
                    x: node.x + 1,
                    y: node.y,
                },
                cost: SCHEME.gap_extension,
            }];
        } else if node.x == self.db.len() {
            self.edges[node.x][node.y] = vec![Edge {
                node: Position {
                    x: node.x,
                    y: node.y + 1,
                },
                cost: SCHEME.gap_extension,
            }];
        } else {
            self.edges[node.x][node.y] = vec![
                Edge {
                    node: Position {
                        x: node.x + 1,
                        y: node.y,
                    },
                    cost: SCHEME.gap_extension,
                },
                Edge {
                    node: Position {
                        x: node.x + 1,
                        y: node.y + 1,
                    },
                    cost: get_cost(&self.query[node.y], &self.query[node.x]),
                },
                Edge {
                    node: Position {
                        x: node.x,
                        y: node.y + 1,
                    },
                    cost: SCHEME.gap_extension,
                },
            ];
        }
    }
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
    fn test_heuristuc() {}
    #[test]
    fn good_graph() {
        let s1 = b"ATGCGATG";
        let s2 = b"CAT";
        assert!(Graph::from_seq(s1, s2).is_ok());
    }
    #[test]
    fn bad_graph() {
        let s1 = b"";
        let s2 = b"ATGTG";
        assert!(Graph::from_seq(s1, s2).is_err());
    }
    #[test]
    fn test_graph_extend() {}
}
