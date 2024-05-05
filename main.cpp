#include "matrix.hpp"
#include "vector"
#include <random>
#include "chrono.hpp"

int main() {
    // Start the timer
    Timings::Chrono clock;
    
    // Create two matrices
    algebra::Matrix<double,algebra::StorageOrder::RowMajor> matrix1;
    algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> matrix2;
    
    // Read matrix1 from file
    matrix1.readMatrix("lnsp_131.mtx");
    std::cout << "Reading matrix1..." << std::endl;
    std::cout << "Number of columns: " << matrix1.get_cols() << std::endl;
    std::cout << "Number of rows: " << matrix1.get_rows() << std::endl;
    std::cout << "Number of non-zero elements: " << matrix1.get_nnz() << std::endl;

    // Read matrix2 from file
    matrix2.readMatrix("lnsp_131.mtx");
    std::cout << "Reading matrix2..." << std::endl;
    std::cout << "Number of columns: " << matrix2.get_cols() << std::endl;
    std::cout << "Number of rows: " << matrix2.get_rows() << std::endl;
    std::cout << "Number of non-zero elements: " << matrix2.get_nnz() << std::endl;

    // Generate a vector of the right dimension
    unsigned int seed = 1234;
    std::default_random_engine generator{seed};
    std::uniform_real_distribution<double> uniform(-10000., 10000.);
    std::vector<double> v(matrix1.get_cols());
    for (std::size_t i=0; i<matrix1.get_cols(); ++i)
        v[i] = uniform(generator);

    std::cout << "\n\nTESTING MATRIX-VECTOR PRODUCT" << std::endl;

    // Result of the product
    std::vector<double> b(matrix1.get_rows());

    // Test in the uncompressed, row-wise storage
    clock.start();
    b = matrix1 * v;
    clock.stop();
    std::cout << "Time for Matrix-Vector multiplication with uncompressed, row-wise storage: "
              << clock << std::endl;

    // Test in the compressed, row-wise storage
    matrix1.compress();
    clock.start();
    b = matrix1 * v;
    clock.stop();
    std::cout << "Time for Matrix-Vector multiplication with compressed, row-wise storage: "
              << clock << std::endl;

    // Test in the uncompressed, column-wise storage
    clock.start();
    b = matrix2 * v;
    clock.stop();
    std::cout << "Time for Matrix-Vector multiplication with uncompressed, column-wise storage: "
              << clock << std::endl;

    // Test in the compressed, column-wise storage
    matrix2.compress();
    clock.start();
    b = matrix2 * v;
    clock.stop();
    std::cout << "Time for Matrix-Vector multiplication with compressed, column-wise storage: "
              << clock << std::endl;

    return 0;
}


