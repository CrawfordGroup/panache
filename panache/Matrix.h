/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef PANACHE_MATRIX_H
#define PANACHE_MATRIX_H

#include <memory>
#include "Exception.h"
#include "Dimension.h"

namespace panache {

class PSIO;
class Matrix;
class Vector;
class Vector3;
typedef std::shared_ptr<Vector> SharedVector;

class SimpleVector;
class MatrixFactory;
class SimpleMatrix;
class Molecule;

enum diagonalize_order {
    evals_only_ascending = 0,
    ascending = 1,
    evals_only_descending = 2,
    descending = 3
};

class Matrix;
typedef std::shared_ptr<Matrix> SharedMatrix;


/*! \ingroup MINTS
 *  \class Matrix
 *  \brief Makes using matrices just a little easier.
 *
 * Using a matrix factory makes creating these a breeze.
 */
class Matrix {
protected:
    /// Matrix data
    double ***matrix_;
    /// Number of irreps
    int nirrep_;
    /// Rows per irrep array
    Dimension rowspi_;
    /// Columns per irrep array
    Dimension colspi_;
    /// Name of the matrix
    std::string name_;
    /// Symmetry of this matrix (in most cases this will be 0 [totally symmetric])
    int symmetry_;

    /// Allocates matrix_
    void alloc();
    /// Release matrix_
    void release();

    /// Copies data from the passed matrix to this matrix_
    void copy_from(double ***);

    /// allocate a block matrix -- analogous to libciomr's block_matrix
    static double** matrix(int nrow, int ncol);
    /// free a (block) matrix -- analogous to libciomr's free_block
    static void free(double** Block);


public:

    /// Default constructor, zeros everything out
    Matrix();
    /**
     * Constructor, zeros everything out, sets name_
     *
     * @param name Name of the matrix, used in saving and printing.
     */
    Matrix(const std::string& name, int symmetry = 0);
    /// copy reference constructor
    Matrix(const Matrix& copy);
    /// Explicit shared point copy constructor
    explicit Matrix(const SharedMatrix& copy);
    /// copy pointer constructor
    explicit Matrix(const Matrix* copy);
    /**
     * Constructor, sets up the matrix
     *
     * @param nirreps Number of blocks.
     * @param rowspi Array of length nirreps giving row dimensionality.
     * @param colspi Array of length nirreps giving column dimensionality.
     */
    Matrix(int nirrep, const int *rowspi, const int *colspi, int symmetry = 0);
    /**
     * Constructor, sets name_, and sets up the matrix
     *
     * @param name Name of the matrix.
     * @param nirreps Number of blocks.
     * @param rowspi Array of length nirreps giving row dimensionality.
     * @param colspi Array of length nirreps giving column dimensionality.
     */
    Matrix(const std::string& name, int nirrep, const int *rowspi, const int *colspi, int symmetry = 0);
    /**
     * Constructor, forms non-standard matrix.
     * @param nirreps Number of blocks.
     * @param rows Singular value. All blocks have same number of rows.
     * @param colspi Array of length nirreps. Defines blocking scheme for columns.
     */
    Matrix(int nirrep, int rows, const int *colspi);

    /**
     * Constructor, forms non-standard matrix.
     * @param nirreps Number of blocks.
     * @param rowspi Array of length nirreps. Defines blocking scheme for rows.
     * @param cols Singular value. All blocks have same number of columns.
     */
    Matrix(int nirrep, const int* rowspi, int cols);

    /**
     * Constructor, sets up the matrix
     * Convenience case for 1 irrep
     * Note: You should be using SimpleMatrix
     *
     * @param rows Row dimensionality.
     * @param cols Column dimensionality.
     */
    Matrix(int rows, int cols);
    /**
     * Constructor, sets up the matrix
     * Convenience case for 1 irrep
     * Note: You should be using SimpleMatrix
     *
     * @param name Name of the matrix.
     * @param rows Row dimensionality.
     * @param cols Column dimensionality.
     */
    Matrix(const std::string&, int rows, int cols);

    /**
     * Constructor using Dimension objects to define order and dimensionality.
     *
     * @param name Name of the matrix.
     * @param rows Dimension object providing row information.
     * @param cols Dimension object providing column information.
     */
    Matrix(const std::string& name, const Dimension& rows, const Dimension& cols, int symmetry = 0);

    /**
     * Constructor using Dimension objects to define order and dimensionality.
     *
     * @param name Name of the matrix.
     * @param rows Dimension object providing row information.
     * @param cols Dimension object providing column information.
     */
    Matrix(const Dimension& rows, const Dimension& cols, int symmetry = 0);

    /// Destructor, frees memory
    ~Matrix();

    /**
     * Initializes a matrix
     *
     * @param nirreps Number of blocks in this matrix.
     * @param rowspi Array of length nirreps giving row dimensionality.
     * @param colspi Array of length nirreps giving column dimensionality.
     */
    void init(int nirrep, const int *rowspi, const int *colspi, const std::string& name = "", int symmetry = 0);

    void init(const Dimension& rowspi, const Dimension& colspi, const std::string& name = "", int symmetry = 0);

    /// Creates an exact copy of the matrix and returns it.
    SharedMatrix clone() const;

    /**
     * Convenient creation function return SharedMatrix
     */
    static SharedMatrix create(const std::string& name,
                               const Dimension& rows,
                               const Dimension& cols);

    /**
     * @{
     * Copies data onto this
     * @param cp Object to copy from.
     */
    void copy(const SharedMatrix& cp);
    void copy(const Matrix& cp);
    void copy(const Matrix* cp);
    /** @} */

    /**
    * Horizontally concatenate matrices
    * @param mats std::vector of Matrix objects to concatenate
    */
    static SharedMatrix horzcat(const std::vector<SharedMatrix >& mats);

    /**
    * Vertically concatenate matrices
    * @param mats std::vector of Matrix objects to concatenate
    */
    static SharedMatrix vertcat(const std::vector<SharedMatrix >& mats);

    /**
    ** For a matrix of 3D vectors (ncol==3), rotate a set of points around an
    ** arbitrary axis.  Vectors are the rows of the matrix.
    **
    ** @param  axis  Vector3   : axis around which to rotate (need not be normalized)
    ** @param  phi   double    : magnitude of rotation in rad
    ** @param  Sn    bool      : if true, then also reflect in plane through origin and
    **                           perpendicular to rotation
    ** @returns SharedMatrix with rotated points (rows)
    */
    SharedMatrix matrix_3d_rotation(Vector3 axis, double phi, bool Sn);

    /// Copies data to the row specified. Assumes data is of correct length.
    void copy_to_row(int h, int row, double const * const data);

    enum SaveType {
        Full,
        SubBlocks,
        LowerTriangle
    };

    /**
     * Set every element of matrix_ to val
     *
     * @param val Value to apply to entire matrix.
     */
    void set(double val);

    /**
     * Copies lower triangle tri to matrix_, calls tri_to_sq
     *
     * @param tri Lower triangle matrix to set to.
     */
    void set(const double * const tri);

    /**
     * @{
     * Copies sq to matrix_
     *
     * @param sq Double matrix to copy over.
     */
    void set(const double * const * const sq);
    /** @} */

    /**
     * @{
     * Copies sq to matrix_
     *
     * @param sq SimpleMatrix object to set this matrix to.
     */
    void set(const SimpleMatrix * const sq);
    void set(const std::shared_ptr<SimpleMatrix>& sq);
    /** @} */

    /**
     * Set a single element of matrix_
     *
     * @param h Subblock to address
     * @param m Row
     * @param n Column
     * @param val Value
     */
    void set(int h, int m, int n, double val) { matrix_[h][m][n] = val; }

    /**
     * Set a single element of matrix_
     *
     * @param m Row
     * @param n Column
     * @param val Value
     */
    void set(int m, int n, double val) { matrix_[0][m][n] = val; }

    /**
     * @{
     * Set the diagonal of matrix_ to vec
     *
     * @param vec Vector to apply to the diagonal.
     */
    void set_diagonal(const Vector * const vec);
    void set_diagonal(const Vector& vec);
    void set_diagonal(const std::shared_ptr<Vector>& vec);
    /** @} */

    /**
     * Returns a single element of matrix_
     *
     * @param h Subblock
     * @param m Row
     * @param n Column
     * @returns value at position (h, m, n)
     */
    double get(const int& h, const int& m, const int& n) const { return matrix_[h][m][n]; }

    /**
     * Returns a single element of matrix_
     *
     * @param h Subblock
     * @param m Row
     * @param n Column
     * @returns value at position (h, m, n)
     */
    double get(const int& m, const int& n) const { return matrix_[0][m][n]; }

    /**
     * Returns a single row of matrix_
     *
     * @param h Subblock
     * @param m Row
     * @returns SharedVector object
     */
    SharedVector get_row(int h, int m);

    /**
     * Returns a single column of matrix_
     *
     * @param h Subblock
     * @param m Column
     * @returns SharedVector object
     */
    SharedVector get_column(int h, int m);

    /**
     * Set a single row of matrix_
     *
     * @param h Subblock
     * @param m Row
     * @returns SharedVector object
     */
    void set_row(int h, int m, SharedVector vec);

    /**
     * Set a single column of matrix_
     *
     * @param h Subblock
     * @param m Column
     * @returns SharedVector object
     */
    void set_column(int h, int m, SharedVector vec);

    /**
     * Returns the double** pointer to the h-th irrep block matrix
     * NOTE: This method is provided for convenience in advanced
     * BLAS/LAPACK calls, and should be used with caution. In particular,
     * operations performed with these pointers should be scoped to avoid
     * erroneous alteration of the objects primitive data. Moreover,
     * the memory location/size of the double** obtained with this method
     * should NEVER be resized, moved, or freed.
     *
     * @param h Subblock
     * @returns pointer to h-th subblock in block-matrix form
     */
    double** pointer(const int& h = 0) const { return matrix_[h]; }
    const double** const_pointer(const int& h=0) const { return const_cast<const double**>(matrix_[h]); }

    /**
     * Returns the double* pointer to the h-th irrep block matrix
     * NOTE: This method is provided for convenience in advanced
     * BLAS/LAPACK calls, and should be used with caution. In particular,
     * operations performed with these pointers should be scoped to avoid
     * erroneous alteration of the objects primitive data. Moreover,
     * the memory location/size of the double* obtained with this method
     * should NEVER be resized, moved, or freed.
     *
     * @param h Subblock
     * @returns pointer to h-th subblock in block-matrix form
     */
    double* get_pointer(const int& h = 0) const {
        if(rowspi_[h]*colspi_[h] > 0)
           return &(matrix_[h][0][0]);
        else
           return 0;}
    const double* get_const_pointer(const int& h=0) const {
        if(rowspi_[h]*colspi_[h] > 0)
           return const_cast<const double*>(&(matrix_[h][0][0]));
        else
           return 0;}

    size_t size(const int &h=0) const { return colspi_[h] * rowspi_[h]; }

    /// apply_denominators a matrix to this
    void apply_denominator(const Matrix* const);
    /// apply_denominators a matrix to this
    void apply_denominator(const Matrix&);
    /// apply_denominators a matrix to this
    void apply_denominator(const SharedMatrix&);

    /**
     * Returns a copy of the current matrix.
     *
     * @returns the matrix
     */
    double **to_block_matrix() const;
    /**
     * Returns a copy of the current matrix in lower triangle form.
     *
     * @returns the matrix
     */
    double *to_lower_triangle() const;

    /**
     * Converts this to a full non-symmetry-block matrix
     *
     * @returns The SimpleMatrix copy of the current matrix.
     */
    SimpleMatrix *to_simple_matrix() const;

    /**
     * Sets the name of the matrix, used in print(...) and save(...)
     *
     * @param name New name to use.
     */
    void set_name(const std::string& name) { name_ = name; }

    /**
     * Gets the name of the matrix.
     */
    const std::string& name() const { return name_; }

    /// Returns the rows in irrep h
    int rowdim(const int& h = 0) const { return rowspi_[h]; }
    /// Returns the cols in irrep h
    int coldim(const int& h = 0) const { return colspi_[h]; }

    /// Returns the rows per irrep array
    const Dimension& rowspi() const {
        return rowspi_;
    }
    /// Returns the rows per irrep array
    int rowspi(const int& h) const {
        return rowdim(h);
    }
    /// Returns the columns per irrep array
    const Dimension& colspi() const {
        return colspi_;
    }
    /// Returns the columns per irrep array
    int colspi(const int& h) const {
        return coldim(h);
    }
    /// Returns the number of irreps
    const int& nirrep() const {
        return nirrep_;
    }

    /// Returns the total number of rows.
    int nrow() const {
        int rows = 0;
        for (int h=0; h<nirrep(); ++h)
            rows += rowdim(h);
        return rows;
    }

    /// Returns the total number of columns.
    int ncol() const {
        int cols = 0;
        for (int h=0; h<nirrep(); ++h)
            cols += coldim(h);
        return cols;
    }

    /// Returns the row size of the largest block.
    int max_nrow() const {
        int row = 0;
        for (int h=0; h<nirrep(); ++h)
            if (row < rowdim(h))
                row = rowdim(h);
        return row;
    }

    /// Returns the column size of the largest block.
    int max_ncol() const {
        int col = 0;
        for (int h=0; h<nirrep(); ++h)
            if (col < coldim(h))
                col = coldim(h);
        return col;
    }

    /**
     * Returns the overall symmetry of the matrix.
     * For a totally-symmetric matrix this will be 0.
     * The value returned is compatible with bitwise XOR (^) math.
     */
    const int& symmetry() const {
        return symmetry_;
    }

    /**
     * Symmetrizes the matrix using information from the given Molecule.
     */
    void symmetrize(std::shared_ptr<Molecule> mol);

    /// Set this to identity
    void identity();
    /// Zeros this out
    void zero();
    /// Zeros the diagonal
    void zero_diagonal();

    // Math routines
    /// Returns the trace of this
    double trace();
    /// Creates a new matrix which is the transpose of this
    SharedMatrix transpose();

    /// In place transposition
    void transpose_this();

    /// Adds a matrix to this
    void add(const Matrix* const);
    /// Adds a matrix to this
    void add(const Matrix&);
    /// Adds a matrix to this
    void add(const SharedMatrix&);

    /// Subtracts a matrix from this
    void subtract(const Matrix* const);
    /// Subtracts a matrix from this
    void subtract(const SharedMatrix&);
    /// Multiplies the two arguments and adds their result to this
    void accumulate_product(const Matrix* const, const Matrix* const);
    void accumulate_product(const SharedMatrix&, const SharedMatrix&);
    /// Scales this matrix
    void scale(double);
    /// Returns the sum of the squares of this
    double sum_of_squares();
    /// Returns the rms of this
    double rms();
    /// Add val to an element of this
    void add(int h, int m, int n, double val) {
        matrix_[h][m][n] += val;
    }
    /// Add val to an element of this
    void add(int m, int n, double val) {
        matrix_[0][m][n] += val;
    }

    void element_add_mirror() {
        for (int h=0; h<nirrep_; ++h) {
            for (int i=0; i<rowspi_[h]; ++i) {
                for (int j=0; j<i; ++j) {
                    matrix_[h][i][j] = matrix_[h][j][i] = (matrix_[h][i][j] + matrix_[h][j][i]);
                }
            }
        }
    }

    /// Scale row m of irrep h by a
    void scale_row(int h, int m, double a);
    /// Scale column n of irrep h by a
    void scale_column(int h, int n, double a);

    /** Special function to transform a SimpleMatrix (no symmetry) into
     *  a symmetry matrix.
     *
     *  \param a SimpleMatrix to transform
     *  \param transformer The matrix returned by PetiteList::aotoso() that acts as the transformer
     */
    void apply_symmetry(const SharedMatrix& a, const SharedMatrix& transformer);

    /** Special function to transform a SimpleMatrix (no symmetry) into
     *  a symmetry matrix.
     *
     *  \param a SimpleMatrix to transform
     *  \param transformer The matrix returned by PetiteList::sotoao() that acts as the transformer
     */
    void remove_symmetry(const SharedMatrix& a, const SharedMatrix& SO2AO);
    /** Performs a the transformation L^ F R. Result goes to this.
     *
     * \param L left transformation matrix (will be transposed)
     * \param F matrix to apply transformation to
     * \param R right transformation matrix (will not be transposed)
     */
    void transform(const SharedMatrix& L,
                   const SharedMatrix& F,
                   const SharedMatrix& R);

    /// @{
    /// Transform a by transformer save result to this
    void transform(const Matrix* const a, const Matrix* const transformer);
    void transform(const SharedMatrix& a, const SharedMatrix& transformer);
    /// @}

    /// @{
    /// Transform this by transformer
    void transform(const Matrix* const transformer);
    void transform(const SharedMatrix& transformer);
    /// @}

    /// @{
    /// Back transform a by transformer save result to this
    void back_transform(const Matrix* const a, const Matrix* const transformer);
    void back_transform(const SharedMatrix& a, const SharedMatrix& transformer);
    /// @}

    /// @{
    /// Back transform this by transformer
    void back_transform(const Matrix* const transformer);
    void back_transform(const SharedMatrix& transformer);
    /// @}

    /// Returns the vector dot product of this by rhs
    double vector_dot(const Matrix* const rhs);
    double vector_dot(const SharedMatrix& rhs);
    double vector_dot(const Matrix& rhs);

    /// @{
    /** General matrix multiply, saves result to this
     * \param transa Transpose the left matrix
     * \param transb Transpose the right matrix
     * \param alpha Prefactor for the matrix multiplication
     * \param a Left matrix
     * \param b Right matrix
     * \param beta Prefactor for the resulting matrix
     */
    void gemm(bool transa, bool transb, double alpha, const Matrix* const a, const Matrix* const b, double beta);
    void gemm(bool transa, bool transb, double alpha, const SharedMatrix& a, const SharedMatrix& b, double beta);
    void gemm(bool transa, bool transb, double alpha, const SharedMatrix& a, const Matrix& b, double beta);
    void gemm(bool transa, bool transb, double alpha, const Matrix& a, const SharedMatrix& b, double beta);
    /// @}

    /// @{
    /** Raw access to the underlying dgemm call. Saves result to this.
     * \param transa Transpose the left matrix
     * \param transb Transpose the right matrix
     * \param
     */
    void gemm(const char& transa, const char& transb,
              const std::vector<int>& m,
              const std::vector<int>& n,
              const std::vector<int>& k,
              const double& alpha,
              const SharedMatrix& a, const std::vector<int>& lda,
              const SharedMatrix& b, const std::vector<int>& ldb,
              const double& beta,
              const std::vector<int>& ldc,
              const std::vector<unsigned long>& offset_a = std::vector<unsigned long>(),
              const std::vector<unsigned long>& offset_b = std::vector<unsigned long>(),
              const std::vector<unsigned long>& offset_c = std::vector<unsigned long>());
    void gemm(const char& transa, const char& transb,
              const int& m,
              const int& n,
              const int& k,
              const double& alpha,
              const SharedMatrix& a, const int& lda,
              const SharedMatrix& b, const int& ldb,
              const double& beta,
              const int& ldc,
              const unsigned long& offset_a = 0,
              const unsigned long& offset_b = 0,
              const unsigned long& offset_c = 0);
    /// @}

    /** Simple doublet GEMM with on-the-fly allocation
    * \param A The first matrix
    * \param B The second matrix
    * \param transA Transpose the first matrix
    * \param transB Transpose the second matrix
    */
    static SharedMatrix doublet(const SharedMatrix& A, const SharedMatrix& B, bool transA = false, bool transB = false);

    /** Simple triplet GEMM with on-the-fly allocation
    * \param A The first matrix
    * \param B The second matrix
    * \param C The third matrix
    * \param transA Transpose the first matrix
    * \param transB Transpose the second matrix
    * \param transC Transpose the third matrix
    */
    static SharedMatrix triplet(const SharedMatrix& A, const SharedMatrix& B, const SharedMatrix& C, bool transA = false, bool transB = false, bool transC = false);

    /** Summation collapse along either rows (0) or columns (1), always producing a column matrix
    * \param dim 0 (row sum) or 1 (col sum)
    * \return \sum_{i} M_{ij} => T_j if dim = 0 or 
    *         \sum_{j} M_{ij} => T_i if dim = 1
    */
    SharedMatrix collapse(int dim = 0);

    /// @{
    /// Diagonalizes this, eigvectors and eigvalues must be created by caller.  Only for symmetric matrices.
    void diagonalize(Matrix* eigvectors, Vector* eigvalues, diagonalize_order nMatz = ascending);
    void diagonalize(SharedMatrix& eigvectors, std::shared_ptr<Vector>& eigvalues, diagonalize_order nMatz = ascending);
    void diagonalize(SharedMatrix& eigvectors, Vector& eigvalues, diagonalize_order nMatz = ascending);
    /// @}

    /// @{
    /// Diagonalizes this, applying supplied metric, eigvectors and eigvalues must be created by caller.  Only for symmetric matrices.
    void diagonalize(SharedMatrix& metric, SharedMatrix& eigvectors, std::shared_ptr<Vector>& eigvalues, diagonalize_order nMatz = ascending);
    /// @}


    /*! Computes the inverse of a real symmetric positive definite
     *  matrix A using the Cholesky factorization A = L*L**T
     *  computed by cholesky_factorize().
     */
    void invert();


    // Reference versions of the above functions
    /// Transform a by transformer save result to this
    void transform(const Matrix& a, const Matrix& transformer);
    /// Transform this by transformer
    void transform(const Matrix& transformer);
    /// Back transform a by transformer save result to this
    void back_transform(const Matrix& a, const Matrix& transformer);
    /// Back transform this by transformer
    void back_transform(const Matrix& transformer);

    /**
     * Expands the row dimension by one, and then orthogonalizes vector v against
     * the current rows, before setting the new row to the orthogonalized copy of v
     */
    bool add_and_orthogonalize_row(const SharedVector v);

    /*! @{
     * Assume this is a orthogonal matrix.  This function Gram-Schmidt
     * orthogonalizes a new vector v and adds it to matrix A. This must contain
     * a free row pointer for a new row.  Don't add orthogonalized v' if
     * norm(v') < NORM_TOL.
     *
     * Adapted from libqt's version by David Sherrill, Feb 1994
     *
     * \param rows current number of valid rows in this
     *             (this must have space for 'rows+1' row.)
     * \param v vector to add to A after it has been made orthogonal
     *             to rest of A
     *
     * \returns true if a vector is added, false otherwise
    */
    bool schmidt_add_row(int h, int rows, Vector& v) throw();
    bool schmidt_add_row(int h, int rows, double* v) throw();
    /// @}

    /*! Calls libqt schmidt function */
    void schmidt();

    /*! Schmidt orthogonalize this. S is the overlap matrix.
     *  n is the number of columns to orthogonalize. */
    void schmidt_orthog(SharedMatrix S, int n);

    /*! Schmidt orthogonalize this. You'll likely want to View this matrix afterwards
     *  using the result to obtain a properly sized Matrix.
     *  \param S overlap matrix.
     *  \param tol is the tolerance.
     *  \returns A Dimension object tell you how many were removed in each irrep.
     */
    Dimension schmidt_orthog_columns(SharedMatrix S, double tol, double*res=0);

    /*!
     * Project out the row vectors in the matrix provided out of this matrix.
     * Assumes all matrices are C1 in nature. Future version will handle irreps.
     * Note: this is destroyed.
     *
     * \param v Matrix to project out
     */
    void project_out(Matrix& v);

    /// General matrix multiply, saves result to this
    void gemm(bool transa, bool transb, double alpha, const Matrix& a, const Matrix& b, double beta);
    /// Diagonalize a symmetric matrix. Eigvectors and eigvalues must be created by caller.
    void diagonalize(Matrix& eigvectors, Vector& eigvalues, int nMatz = 1);

    /// @{
    /// Retrieves the i'th irrep
    double** operator[](int i) { return matrix_[i]; }
    double& operator()(int i, int j) { return matrix_[0][i][j]; }
    const double& operator()(int i, int j) const { return matrix_[0][i][j]; }
    double& operator()(int h, int i, int j) { return matrix_[h][i][j]; }
    const double& operator()(int h, int i, int j) const { return matrix_[h][i][j]; }
    /// @}

    // Serializable pure virtual functions:
    void send();
    void recv();

    /// @{
    /// Checks matrix equality.
    /// @param rhs Matrix to compare to.
    /// @returns true if equal, otherwise false.
    bool equal(const Matrix& rhs);
    bool equal(const SharedMatrix& rhs);
    bool equal(const Matrix* rhs);
    /// @}

    /// @{
    /// Checks matrix equality, but allows rows to be in a different order.
    /// @param rhs Matrix to compare to.
    /// @returns true if equal, otherwise false.
    bool equal_but_for_row_order(const Matrix& rhs, double TOL=1.0e-10);
    bool equal_but_for_row_order(const SharedMatrix& rhs, double TOL=1.0e-10);
    bool equal_but_for_row_order(const Matrix* rhs, double TOL=1.0e-10);
    /// @}

    /**
     * Rotates columns i and j in irrep h, by an angle theta
     * @param h - the irrep in which the rotation will be applied
     * @param i - the zero-based (within irrep) column number for i
     * @param j - the zero-based (within irrep) column number for j
     * @param theta - the angle (in radians) about which to rotate
     */
    void rotate_columns(int h, int i, int j, double theta);
    friend class Vector;

    double * give_up();
};

}

#endif //PANACHE_MATRIX_H
