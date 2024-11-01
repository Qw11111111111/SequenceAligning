use crate::errors::{AStarError, Result};
use clap::Parser;
use std::fmt::Display;
use std::fs::{read, write};
use std::io;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// Path to query sequence
    #[arg(short, long)]
    pub query_file: PathBuf,

    /// path to db sequence
    #[arg(short, long)]
    pub db_file: PathBuf,

    /// out path
    #[arg(short, long, default_value = "./results")]
    pub out_path: PathBuf,

    /// verbose
    #[arg(short, long, default_value = "false")]
    pub verbose: bool,
    //TODO: enable more modi
    /// modus
    #[arg(short, long, default_value = "false")]
    pub local: bool,
}

const ALLOWED_CHARS: [u8; 5] = [b'A', b'G', b'C', b'T', b'N'];

pub fn parse_fasta(path: PathBuf) -> Result<Records> {
    if !path.ends_with("fa") && !path.ends_with("fasta") {
        return Err(AStarError::FastaError(io::Error::from(
            io::ErrorKind::InvalidInput,
        )));
    }
    let mut recs = Records::default();
    let contents: Vec<u8> = read(path)?;
    let mut current_rec = Record::default();
    let mut in_name = false;
    let mut err_chars = vec![];
    for &c in &contents {
        if c == b'>' {
            recs.records.push(current_rec);
            current_rec = Record {
                seq: Vec::default(),
                name: vec![c],
            };
            in_name = true;
        }
        if in_name {
            if c == b'\n' {
                in_name = false;
                continue;
            }
            current_rec.name.push(c);
        } else if c == b'\n' {
            continue;
        } else if ALLOWED_CHARS.contains(&c) {
            err_chars.push(char::from_u32(c as u32).unwrap_or('?'));
        } else {
            current_rec.seq.push(c);
        }
    }
    recs.records.remove(0);
    if !err_chars.is_empty() {
        return Err(AStarError::CharError {
            res: recs,
            chars: err_chars,
        });
    }
    Ok(recs)
}

#[derive(Default)]
pub struct Records {
    pub records: Vec<Record>,
}

impl Display for Records {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for rec in &self.records {
            write!(f, "{}", rec)?;
        }
        Ok(())
    }
}

impl Iterator for Records {
    type Item = Record;
    fn next(&mut self) -> Option<Self::Item> {
        self.records.pop()
    }
}

impl Records {
    fn write_to_fa(&self, path: PathBuf) -> Result<()> {
        write(path, format!("{}", self))?;
        Ok(())
    }
}

#[derive(Default, PartialEq, PartialOrd, Ord, Eq)]
pub struct Record {
    seq: Vec<u8>,
    name: Vec<u8>,
}

impl Display for Record {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, ">")?;
        for c in &self.name {
            write!(f, "{}", char::from_u32(*c as u32).unwrap())?;
        }
        writeln!(f)?;
        for c in &self.seq {
            write!(f, "{}", char::from_u32(*c as u32).unwrap())?;
        }
        writeln!(f)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::panic;
    use std::fs;
    use std::io;
    use std::io::Write;
    use tempfile::NamedTempFile;

    //TODO
    #[test]
    fn parse_good_fasta() -> io::Result<()> {
        let mut t_file = NamedTempFile::new()?;
        write!(
            t_file,
            ">Record1\nATGCATGCATGCATGCATGCATGCATGC\n>Record2\nATGCATGCGTGCAGTGACCACA"
        )?;
        t_file.flush()?;
        let p = t_file.path().to_path_buf();
        assert!(parse_fasta(p).is_ok());
        Ok(())
    }

    #[test]
    fn parse_bad_header() -> io::Result<()> {
        let mut t_file = NamedTempFile::new()?;
        write!(
            t_file,
            ">Record1\nATGCATGCATGCATGCATGCATGCATGC\nRecord2\nATGCATGCGTGCAGTGACCACA"
        )?;
        t_file.flush()?;
        let p = t_file.path().to_path_buf();
        match parse_fasta(p) {
            Err(AStarError::CharError { res, chars }) => {
                assert_eq!(chars, vec!['R', 'e', 'c', 'o', 'r', 'd', '2']);
                if let Some(rec) = res.records.first() {
                    assert_eq!(rec.name, b">Record1");
                    assert_eq!(
                        rec.seq,
                        b"ATGCATGCATGCATGCATGCATGCATGCATGCATGCGTGCAGTGACCACA"
                    )
                } else {
                    panic!("records was empty, but should conatin 1 element");
                }
            }
            _ => panic!("expected char error"),
        }
        Ok(())
    }

    #[test]
    fn parse_bad_nt() -> io::Result<()> {
        let mut t_file = NamedTempFile::new()?;
        write!(t_file, ">Record1\nATGCATGCAKGCATGCATGCANNNGCATGC")?;
        t_file.flush()?;
        let p = t_file.path().to_path_buf();
        match parse_fasta(p) {
            Err(AStarError::CharError { res, chars }) => {
                assert_eq!(chars, vec!['K']);
                if let Some(rec) = res.records.first() {
                    assert_eq!(rec.name, b">Record1");
                    assert_eq!(rec.seq, b"ATGCATGCAGCATGCATGCANNNGCATGC")
                } else {
                    panic!("records was empty, but should conatin 1 element");
                }
            }
            _ => panic!("expected char error"),
        }
        Ok(())
    }

    #[test]
    fn parse_false_file() -> io::Result<()> {
        let t_file = NamedTempFile::new()?;
        let original_p = t_file.path().to_path_buf();
        let new_p = original_p.with_extension("txt");
        fs::rename(original_p, &new_p)?;
        match parse_fasta(new_p) {
            Err(AStarError::FastaError(_)) => {}
            _ => panic!("should panic on non fasta file"),
        }
        Ok(())
    }
}