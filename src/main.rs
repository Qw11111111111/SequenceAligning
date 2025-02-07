mod align;
pub mod errors;
#[allow(dead_code, unused_imports)]
//mod needleman_wunsch;
mod needleman_wunsch_affine;
mod parse;
pub mod utils;
mod wfa;

use align::align;
use clap::Parser;
use errors::AlignerError;
use needleman_wunsch_affine::n_w_align;
//use needleman_wunsch::n_w_align;
use parse::{parse_fasta, Algo, Args};
use wfa::wfa_align;

//TODO handle errors appropriatly
fn main() {
    let args = Args::parse();

    let db = match parse_fasta(args.db_file) {
        Ok(rec) => rec,
        Err(AlignerError::FastaError(e)) => {
            eprintln!("DB fasta could not be opened: {}", e);
            eprintln!("aborting");
            return;
        }
        Err(AlignerError::CharError { chars, res }) => {
            eprintln!(
                "Invalid character '{:#?}' detected in db fasta; continuing by ignoring it",
                chars
            );
            res
        }
        Err(e) => {
            eprintln!("Unexpected error in db fasta: {}", e);
            return;
        }
    };

    let query = match parse_fasta(args.query_file) {
        Ok(q) => q,
        Err(AlignerError::FastaError(e)) => {
            eprintln!("Query fasta could not be opened: {}", e);
            eprintln!("aborting");
            return;
        }
        Err(AlignerError::CharError { chars, res }) => {
            eprintln!(
                "Invalid character '{:#?}' detected in query fasta; continuing by ignoring it",
                chars
            );
            res
        }
        Err(e) => {
            eprintln!("Unexpected error in query fasta: {}", e);
            return;
        }
    };
    for d in db.records.iter() {
        for q in query.records.iter() {
            match match args.algo {
                Algo::AStar => align(q, d, args.verbose, false),
                Algo::NeedlemanWunsch => n_w_align(q, d, args.verbose, args.mode.clone()),
                Algo::Wfa => wfa_align(q, d, args.mode.clone()),
            } {
                Err(AlignerError::AlignmentError(e)) => {
                    eprintln!(
                        "An error occured during alignment of {} and {}\n{e}",
                        q.name.iter().map(|i| *i as char).collect::<String>(),
                        d.name.iter().map(|i| *i as char).collect::<String>()
                    )
                }
                Ok(_) => {}
                _ => {}
            }
        }
    }
}
