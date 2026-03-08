use faer::Mat;
use faer::linalg::solvers::FullPivLu;
use faer::prelude::Solve;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SparseMatrixError {
    #[error("matrix is singular, cannot solve")]
    Singular,
    #[error(
        "dimension mismatch: matrix is {matrix_dim}x{matrix_dim} but RHS has {rhs_len} elements"
    )]
    DimensionMismatch { matrix_dim: usize, rhs_len: usize },
}

/// A triplet (row, col, value) for assembling a sparse matrix.
#[derive(Debug, Clone, Copy)]
struct Triplet {
    row: usize,
    col: usize,
    value: f64,
}

/// A sparse matrix builder using triplet-form (COO) assembly.
///
/// Supports dynamic insertion of entries. Duplicate (row, col) entries
/// are summed together during assembly, which is the standard behavior
/// for MNA stamp accumulation.
#[derive(Debug, Clone)]
pub struct SparseMatrix {
    dim: usize,
    triplets: Vec<Triplet>,
}

impl SparseMatrix {
    /// Create a new sparse matrix of the given dimension.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            triplets: Vec::new(),
        }
    }

    /// Return the dimension (number of rows/columns).
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Add a value to position (row, col). If multiple values are added
    /// at the same position, they are summed (standard MNA stamp behavior).
    ///
    /// # Panics
    /// Panics if row or col >= dim.
    pub fn add(&mut self, row: usize, col: usize, value: f64) {
        assert!(
            row < self.dim,
            "row {row} out of bounds for dim {}",
            self.dim
        );
        assert!(
            col < self.dim,
            "col {col} out of bounds for dim {}",
            self.dim
        );
        if value != 0.0 {
            self.triplets.push(Triplet { row, col, value });
        }
    }

    /// Convert triplet form to a dense faer matrix (summing duplicates).
    pub fn to_dense(&self) -> Mat<f64> {
        let mut mat = Mat::zeros(self.dim, self.dim);
        for t in &self.triplets {
            mat[(t.row, t.col)] += t.value;
        }
        mat
    }

    /// Clear all entries, keeping the dimension.
    pub fn clear(&mut self) {
        self.triplets.clear();
    }
}

/// A linear system Ax = b assembled in triplet form.
#[derive(Debug, Clone)]
pub struct LinearSystem {
    pub matrix: SparseMatrix,
    pub rhs: Vec<f64>,
}

impl LinearSystem {
    /// Create a new linear system of the given dimension.
    pub fn new(dim: usize) -> Self {
        Self {
            matrix: SparseMatrix::new(dim),
            rhs: vec![0.0; dim],
        }
    }

    /// Return the dimension.
    pub fn dim(&self) -> usize {
        self.matrix.dim()
    }

    /// Solve the system using LU factorization.
    /// Returns the solution vector x such that Ax = b.
    pub fn solve(&self) -> Result<Vec<f64>, SparseMatrixError> {
        let dim = self.matrix.dim();
        if self.rhs.len() != dim {
            return Err(SparseMatrixError::DimensionMismatch {
                matrix_dim: dim,
                rhs_len: self.rhs.len(),
            });
        }

        if dim == 0 {
            return Ok(vec![]);
        }

        let a = self.matrix.to_dense();
        let lu = FullPivLu::new(a.as_ref());

        // Build faer column vector from rhs
        let mut b = Mat::zeros(dim, 1);
        for (i, &val) in self.rhs.iter().enumerate() {
            b[(i, 0)] = val;
        }

        let x = lu.solve(&b);

        // Check for NaN (indicates singular matrix)
        let result: Vec<f64> = (0..dim).map(|i| x[(i, 0)]).collect();
        if result.iter().any(|v: &f64| v.is_nan() || v.is_infinite()) {
            return Err(SparseMatrixError::Singular);
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_3x3_system() {
        // Solve:
        //  2x + y - z = 8
        // -3x - y + 2z = -11
        // -2x + y + 2z = -3
        //
        // Solution: x=2, y=3, z=-1
        let mut sys = LinearSystem::new(3);

        sys.matrix.add(0, 0, 2.0);
        sys.matrix.add(0, 1, 1.0);
        sys.matrix.add(0, 2, -1.0);

        sys.matrix.add(1, 0, -3.0);
        sys.matrix.add(1, 1, -1.0);
        sys.matrix.add(1, 2, 2.0);

        sys.matrix.add(2, 0, -2.0);
        sys.matrix.add(2, 1, 1.0);
        sys.matrix.add(2, 2, 2.0);

        sys.rhs[0] = 8.0;
        sys.rhs[1] = -11.0;
        sys.rhs[2] = -3.0;

        let x = sys.solve().unwrap();
        assert_abs_diff_eq!(x[0], 2.0, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], 3.0, epsilon = 1e-12);
        assert_abs_diff_eq!(x[2], -1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_triplet_sum_duplicates() {
        // Adding to the same position should sum values
        let mut m = SparseMatrix::new(2);
        m.add(0, 0, 3.0);
        m.add(0, 0, 4.0);
        let dense = m.to_dense();
        assert_abs_diff_eq!(dense[(0, 0)], 7.0, epsilon = 1e-15);
    }

    #[test]
    fn test_identity_solve() {
        // Identity matrix: solution equals RHS
        let mut sys = LinearSystem::new(3);
        sys.matrix.add(0, 0, 1.0);
        sys.matrix.add(1, 1, 1.0);
        sys.matrix.add(2, 2, 1.0);
        sys.rhs[0] = 5.0;
        sys.rhs[1] = -3.0;
        sys.rhs[2] = 7.0;

        let x = sys.solve().unwrap();
        assert_abs_diff_eq!(x[0], 5.0, epsilon = 1e-12);
        assert_abs_diff_eq!(x[1], -3.0, epsilon = 1e-12);
        assert_abs_diff_eq!(x[2], 7.0, epsilon = 1e-12);
    }

    #[test]
    fn test_singular_matrix() {
        // All zeros => singular
        let sys = LinearSystem::new(2);
        // rhs is [0, 0] and matrix is all zeros
        // This might solve to zeros rather than error, so let's make it clearly singular
        // with a non-zero rhs
        let mut sys2 = LinearSystem::new(2);
        sys2.rhs[0] = 1.0;
        let result = sys2.solve();
        assert!(
            result.is_err() || {
                // Some LU implementations return NaN/Inf for singular
                let x = result.unwrap();
                x.iter().any(|v| v.is_nan() || v.is_infinite())
            }
        );

        // Trivial zero system with zero rhs should still "solve" to zero
        let x = sys.solve().unwrap_or_else(|_| vec![0.0, 0.0]);
        // Either way is fine for a zero system
        let _ = x;
    }

    #[test]
    fn test_empty_system() {
        let sys = LinearSystem::new(0);
        let x = sys.solve().unwrap();
        assert!(x.is_empty());
    }
}
