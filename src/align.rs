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
        cost: 0,
        position: Position { x: 0, y: 0 },
    });
    while let Some(s) = queue.pop() {
        graph.extend(&s.position);
        //TODO calculate cost of node as heuristic + cost_to reach
        // also need to somehow save the cumulative cost up to that point.
    }
    Ok(())
}

struct ScoringScheme {
    match_score: i16,
    mismatch_score: i16,
    gap_opening: i16,
    gap_extension: i16,
}

fn heuristic_d(seq1: &Vec<u8>, seq2: &Vec<u8>, x: usize, y: usize) -> u16 {
    0
}

#[derive(Debug, Default, PartialEq, Eq)]
struct State {
    cost: usize,
    position: Position,
}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .cmp(&self.cost)
            .then_with(|| self.position.cmp(&other.position))
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
    cost: i16,
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

fn get_cost(c1: &u8, c2: &u8) -> i16 {
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
