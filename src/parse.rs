use crate::errors::{AlignerError, Result};
use clap::{Parser, ValueEnum};
use std::fmt::Display;
use std::fs::{read, write};
use std::io;
use std::path::{Path, PathBuf};

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
    #[arg(short, long, value_enum, default_value_t = Mode::Global)]
    pub mode: Mode,

    /// algo
    #[arg(short, long, value_enum, default_value_t = Algo::AStar)]
    pub algo: Algo,
}

#[derive(ValueEnum, Debug, Default, Clone)]
pub enum Algo {
    #[default]
    AStar,
    NeedlemanWunsch,
    Wfa,
}

#[derive(ValueEnum, Debug, Default, Clone)]
pub enum Mode {
    #[default]
    Global,
    Local,
    SemiGlobal,
}

const ALLOWED_CHARS: [u8; 5] = [b'A', b'G', b'C', b'T', b'N'];

pub fn parse_fasta<'a>(path: PathBuf) -> Result<'a, Records> {
    if !(has_extension(&path, "fa") || has_extension(&path, "fasta") || has_extension(&path, "fna"))
    {
        return Err(AlignerError::FastaError(io::Error::from(
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
            continue;
        }
        if in_name {
            if c == b'\n' {
                in_name = false;
                continue;
            }
            current_rec.name.push(c);
        } else if c == b'\n' {
            continue;
        } else if !ALLOWED_CHARS.contains(&c) {
            err_chars.push(char::from_u32(c as u32).unwrap_or('?'));
        } else {
            current_rec.seq.push(c);
        }
    }
    recs.records.push(current_rec);
    recs.records.remove(0);
    if !err_chars.is_empty() {
        return Err(AlignerError::CharError {
            res: recs,
            chars: err_chars,
        });
    }
    Ok(recs)
}

fn has_extension(path: &Path, ext: &str) -> bool {
    match path.extension() {
        Some(extension) => extension == ext,
        None => false,
    }
}
#[derive(Default, Debug)]
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
    fn _write_to_fa(&self, path: PathBuf) -> Result<()> {
        write(path, format!("{}", self))?;
        Ok(())
    }
}

#[derive(Default, PartialEq, PartialOrd, Ord, Eq, Debug)]
pub struct Record {
    pub seq: Vec<u8>,
    pub name: Vec<u8>,
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
        let original_p = t_file.path().to_path_buf();
        let new_p = original_p.with_extension("fa");
        fs::rename(original_p, &new_p)?;
        write!(
            t_file,
            ">Record1\nATGCATGCATGCATGCATGCATGCATGC\n>Record2\nATGCATGCGTGCAGTGACCACA"
        )?;
        t_file.flush()?;
        match parse_fasta(new_p) {
            Ok(res) => {
                assert_eq!(res.records.len(), 2);
                assert_eq!(res.records[0].name.len(), 8, "{:#?}", res.records[0].name);
                assert_eq!(res.records[0].seq.len(), 28);
            }
            Err(e) => panic!("should return Ok not {:#?}", e),
        }
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
        let original_p = t_file.path().to_path_buf();
        let new_p = original_p.with_extension("fa");
        fs::rename(original_p, &new_p)?;
        match parse_fasta(new_p) {
            Err(AlignerError::CharError { res, chars }) => {
                assert_eq!(chars, vec!['R', 'e', 'c', 'o', 'r', 'd', '2']);
                if let Some(rec) = res.records.first() {
                    assert_eq!(rec.name, b">Record1");
                    assert_eq!(
                        rec.seq,
                        b"ATGCATGCATGCATGCATGCATGCATGCATGCATGCGTGCAGTGACCACA"
                    )
                } else {
                    panic!("records was empty, but should contain 1 element");
                }
            }
            e => panic!("expected char error, got {:#?}", e),
        }
        Ok(())
    }

    #[test]
    fn parse_bad_nt() -> io::Result<()> {
        let mut t_file = NamedTempFile::new()?;
        write!(t_file, ">Record1\nATGCATGCAKGCATGCATGCANNNGCATGC")?;
        t_file.flush()?;
        let original_p = t_file.path().to_path_buf();
        let new_p = original_p.with_extension("fa");
        fs::rename(original_p, &new_p)?;
        match parse_fasta(new_p) {
            Err(AlignerError::CharError { res, chars }) => {
                assert_eq!(chars, vec!['K']);
                if let Some(rec) = res.records.first() {
                    assert_eq!(rec.name, b">Record1");
                    assert_eq!(rec.seq, b"ATGCATGCAGCATGCATGCANNNGCATGC")
                } else {
                    panic!("records was empty, but should conatin 1 element");
                }
            }
            e => panic!("expected char error, got: {:#?}", e),
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
            Err(AlignerError::FastaError(_)) => {}
            _ => panic!("should panic on non fasta file"),
        }
        Ok(())
    }
}
