use std::fmt;
use std::io;

/// Error categories and details for MAGeCKrs.
#[derive(Debug)]
pub enum MageckError {
    /// Configuration or CLI argument validation error.
    Config(String),
    /// Library parsing or index construction error.
    Library(String),
    /// FASTQ stream or format error.
    Fastq(String),
    /// Standard I/O error.
    Io(io::Error),
    /// Design matrix format or rank error.
    DesignMatrix(String),
    /// Pathway GMT or parsing error.
    Pathway(String),
    /// Math or modeling error.
    Math(String),
}

impl fmt::Display for MageckError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MageckError::Config(msg) => write!(f, "Configuration error: {msg}"),
            MageckError::Library(msg) => write!(f, "Library error: {msg}"),
            MageckError::Fastq(msg) => write!(f, "FASTQ error: {msg}"),
            MageckError::Io(err) => write!(f, "I/O error: {err}"),
            MageckError::DesignMatrix(msg) => write!(f, "Design matrix error: {msg}"),
            MageckError::Pathway(msg) => write!(f, "Pathway error: {msg}"),
            MageckError::Math(msg) => write!(f, "Math/Modeling error: {msg}"),
        }
    }
}

impl std::error::Error for MageckError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            MageckError::Io(err) => Some(err),
            _ => None,
        }
    }
}

impl From<io::Error> for MageckError {
    fn from(err: io::Error) -> Self {
        MageckError::Io(err)
    }
}

impl From<crate::stats::cnv::CnvError> for MageckError {
    fn from(error: crate::stats::cnv::CnvError) -> Self {
        Self::Math(error.to_string())
    }
}
