#include "matrix.hpp"
#include "vector"
#include <random>
#include "chrono.hpp"

int main() {
    Timings::Chrono clock;
    
    // Create two matrices
    algebra::Matrix<double,algebra::StorageOrder::RowMajor> matrix1;
    algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> matrix2;

    // Read matrix from file
    matrix1.readMatrix("lnsp_131.mtx");

    matrix2.readMatrix("lnsp_131.mtx");

    // Generate a vector of the right dimension
    unsigned int seed = 1234;
    std::default_random_engine generator{seed};
    std::uniform_real_distribution<double> uniform(-10000., 10000.);
    algebra::Matrix<double,algebra::StorageOrder::RowMajor> rowvec(matrix1.get_cols(),1);
    algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> colvec(matrix1.get_cols(),1);
    for (std::size_t i=0; i<matrix1.get_cols(); ++i){
        rowvec(i,0) = uniform(generator);
        colvec(i,0) = uniform(generator);
    }

    std::cout << "\n\nTESTING MATRIX - VECTOR PRODUCT" << std::endl;

    // Result of the product
    std::vector<double> b(matrix1.get_rows());

    // Test in the uncompressed, row-wise storage
    clock.start();
    b = matrix1 * rowvec;
    clock.stop();
    std::cout << "Time for Matrix-Vector multiplication with uncompressed, row-wise storage: "
              << clock << std::endl;

    // Test in the compressed, row-wise storage
    matrix1.compress();
    clock.start();
    b = matrix1 * rowvec;
    clock.stop();
    std::cout << "Time for Matrix-Vector multiplication with compressed, row-wise storage: "
              << clock << std::endl;

    // Test in the uncompressed, column-wise storage
    clock.start();
    b = matrix2 * colvec;
    clock.stop();
    std::cout << "Time for Matrix-Vector multiplication with uncompressed, column-wise storage: "
              << clock << std::endl;

    // Test in the compressed, column-wise storage
    matrix2.compress();
    clock.start();
    b = matrix2 * colvec;
    clock.stop();
    std::cout << "Time for Matrix-Vector multiplication with compressed, column-wise storage: "
              << clock << std::endl;

    std::cout << "\n\nTESTING MATRIX NORMS" << std::endl;

    matrix1.uncompress();
    matrix2.uncompress();

    // Uncompressed, row-wise storage
    std::cout << "Uncompressed, row-wise storage:" << std::endl;
    std::cout << "Infinity norm: " << matrix1.norm<algebra::Norm::InfinityNorm>() << std::endl;
    std::cout << "1-norm: " << matrix1.norm<algebra::Norm::OneNorm>() << std::endl;
    std::cout << "Frobenius-norm: " << matrix1.norm<algebra::Norm::FrobeniusNorm>() << std::endl;

    // Compressed, row-wise storage
    std::cout << "\nCompressed, row-wise storage:" << std::endl;
    matrix1.compress();
    std::cout << "Infinity norm: " << matrix1.norm<algebra::Norm::InfinityNorm>() << std::endl;
    std::cout << "1-norm: " << matrix1.norm<algebra::Norm::OneNorm>() << std::endl;
    std::cout << "Frobenius-norm: " << matrix1.norm<algebra::Norm::FrobeniusNorm>() << std::endl;

    // Uncompressed, column-wise storage
    std::cout << "\nUncompressed, column-wise storage:" << std::endl;
    std::cout << "Infinity norm: " << matrix2.norm<algebra::Norm::InfinityNorm>() << std::endl;
    std::cout << "1-norm: " << matrix2.norm<algebra::Norm::OneNorm>() << std::endl;
    std::cout << "Frobenius-norm: " << matrix2.norm<algebra::Norm::FrobeniusNorm>() << std::endl;

    // Compressed, column-wise storage
    std::cout << "\nCompressed, column-wise storage:" << std::endl;
    matrix2.compress();
    std::cout << "Infinity norm: " << matrix2.norm<algebra::Norm::InfinityNorm>() << std::endl;
    std::cout << "1-norm: " << matrix2.norm<algebra::Norm::OneNorm>() << std::endl;
    std::cout << "Frobenius-norm: " << matrix2.norm<algebra::Norm::FrobeniusNorm>() << std::endl;

    return 0;
}




