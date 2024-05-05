#include "matrix.hpp"
#include "vector"
#include <random>
#include "chrono.hpp"

int main() {
Timings::Chrono clock;
    
algebra::Matrix<double,algebra::StorageOrder::RowMajor> matrix3;
algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> matrix4;

matrix3.readMatrix("lnsp_131.mtx");
std::cout << matrix3.get_cols() << std::endl;
std::cout << matrix3.get_rows() << std::endl;
std::cout << matrix3.get_nnz() << std::endl;

matrix4.readMatrix("lnsp_131.mtx");
std::cout << matrix4.get_cols() << std::endl;
std::cout << matrix4.get_rows() << std::endl;
std::cout << matrix4.get_nnz() << std::endl;



//Generate a vector of the right dimension
    unsigned int seed = 1234;
    std::default_random_engine generator{seed};
    std::uniform_real_distribution<double> uniform(-10000., 10000.);
    algebra::Matrix<double,algebra::StorageOrder::RowMajor> rowvec(matrix3.get_cols(),1);
    algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> colvec(matrix3.get_cols(),1);
    for (std::size_t i=0; i<matrix3.get_cols(); ++i){
        rowvec(i,0) = uniform(generator);
        colvec(i,0) = uniform(generator);
    }

    std::cout << "\n\nTESTING MATRIX - VECTOR PRODUCT" << std::endl;

    // result of the product
    std::vector<double> b(matrix3.get_rows());

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    clock.start();
    b = matrix3 * rowvec;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with uncompressed, row-wise storage: "
              << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    matrix3.compress();
    clock.start();
    b = matrix3 * rowvec;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with compressed, row-wise storage: "
              << clock << std::endl;


    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    clock.start();
    b = matrix4 * colvec;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with uncompressed, column-wise storage: "
              << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    matrix4.compress();
    clock.start();
    b = matrix4 * colvec;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with compressed, column-wise storage: "
              << clock << std::endl;


            std::cout << "\n\nTESTING MATRIX NORMS" << std::endl;

            matrix3.uncompress();
            matrix4.uncompress();

            // Uncompressed, row-wise storage
            std::cout << "Uncompressed, row-wise storage:" << std::endl;
            std::cout << "Infinity norm: " << matrix3.norm<algebra::Norm::InfinityNorm>() << std::endl;
            std::cout << "1-norm: " << matrix3.norm<algebra::Norm::OneNorm>() << std::endl;
            std::cout << "Frobenius-norm: " << matrix3.norm<algebra::Norm::FrobeniusNorm>() << std::endl;

            // Compressed, row-wise storage
            std::cout << "Compressed, row-wise storage:" << std::endl;
            matrix3.compress();
            std::cout << "Infinity norm: " << matrix3.norm<algebra::Norm::InfinityNorm>() << std::endl;
            std::cout << "1-norm: " << matrix3.norm<algebra::Norm::OneNorm>() << std::endl;
            std::cout << "Frobenius-norm: " << matrix3.norm<algebra::Norm::FrobeniusNorm>() << std::endl;

            // Uncompressed, column-wise storage
            std::cout << "Uncompressed, column-wise storage:" << std::endl;
            std::cout << "Infinity norm: " << matrix4.norm<algebra::Norm::InfinityNorm>() << std::endl;
            std::cout << "1-norm: " << matrix4.norm<algebra::Norm::OneNorm>() << std::endl;
            std::cout << "Frobenius-norm: " << matrix4.norm<algebra::Norm::FrobeniusNorm>() << std::endl;

            // Compressed, column-wise storage
            std::cout << "Compressed, column-wise storage:" << std::endl;
            matrix4.compress();
            std::cout << "Infinity norm: " << matrix4.norm<algebra::Norm::InfinityNorm>() << std::endl;
            std::cout << "1-norm: " << matrix4.norm<algebra::Norm::OneNorm>() << std::endl;
            std::cout << "Frobenius-norm: " << matrix4.norm<algebra::Norm::FrobeniusNorm>() << std::endl;

return 0;
}




