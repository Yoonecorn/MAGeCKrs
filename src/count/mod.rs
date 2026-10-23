pub mod matcher;
pub mod reader;

pub use matcher::{
    build_index, parse_library, parse_library_file, GuideIndex, Library, SgRnaRecord,
};
pub use reader::{
    count_all_fastq, count_grouped_fastq, GroupedCountResult, PhysicalFastqSummary, ReadSummary,
};
pub mod normalize;
pub use normalize::normalize;
pub mod writer;
pub use writer::{write_count_outputs, CountOutputPaths};
pub mod control;
pub use control::{parse_control_sgrna, parse_control_sgrna_file, ResolvedControlGuides};
pub mod run;
pub use run::run;
