mod align;
pub mod errors;
mod parse;

use clap::Parser;
use errors::{AStarError, Result};
use parse::{parse_fasta, Args, Records};

fn main() {
    let args = Args::parse();

    let db = match parse_fasta(args.db_file) {
        Ok(rec) => rec,
        Err(AStarError::FastaError(e)) => {
            eprintln!("DB fasta could not be opened: {}", e);
            eprintln!("aborting");
            return;
        }
        Err(AStarError::CharError { chars, res }) => {
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
        Err(AStarError::FastaError(e)) => {
            eprintln!("Query fasta could not be opened: {}", e);
            eprintln!("aborting");
            return;
        }
        Err(AStarError::CharError { chars, res }) => {
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
}
