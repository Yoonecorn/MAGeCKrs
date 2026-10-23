pub mod count_matrix;
pub mod packed_dna;
pub mod parsed_table;
pub mod umi;

pub use count_matrix::{CountMatrix, CountMatrixError, StridedIter};
pub use packed_dna::{ByteGuideIndex, EncodingError, MatchResult, PackedDna, PackedGuideIndex};
pub use parsed_table::ParsedCountTable;
pub use umi::{UmiDeduper, UmiError};
