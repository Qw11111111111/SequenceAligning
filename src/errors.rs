use thiserror::Error;

pub type Result<T> = std::result::Result<T, AStarError<T>>;

#[derive(Error, Debug)]
pub enum AStarError<T> {
    #[error("Fasta could not be opened with err: {0}")]
    FastaError(#[from] std::io::Error),
    #[error("Error in alignment: {0}")]
    AlignmentError(String),
    #[error("Invalid character: {chars:?}")]
    CharError { res: T, chars: Vec<char> },
}
