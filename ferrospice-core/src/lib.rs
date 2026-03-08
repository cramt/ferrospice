pub mod mna;
pub mod sparse;

pub use mna::{MnaError, MnaSystem, assemble_mna};
pub use sparse::{LinearSystem, SparseMatrix, SparseMatrixError};
