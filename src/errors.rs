use thiserror::Error;

pub type Result<'a, T> = std::result::Result<T, AStarError<'a, T>>;

//TODO more error types (unrecoverable, recoverable,...)

#[derive(Error, Debug)]
pub enum AStarError<'a, T> {
    #[error("Fasta could not be opened with err: {0}")]
    FastaError(#[from] std::io::Error),
    #[error("Error in alignment: {0}")]
    AlignmentError(&'a str),
    #[error("Invalid character: {chars:?}")]
    CharError { res: T, chars: Vec<char> },
}
