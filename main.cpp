#include "matrix.hpp"
#include "vector"
#include <random>
#include "chrono.hpp"

int main() {
Timings::Chrono clock;
    // Create an instance of your matrix class
//algebra::Matrix<double,algebra::StorageOrder::RowMajor> matrix1(4,6);
//algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> matrix2(3,3);
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




// Read the matrix from the file
// matrix1(0,0)=10;
// matrix1(0,1)=20;
// matrix1(1,1)=30;
// matrix1(1,3)=40;
// matrix1(2,2)=50;
// matrix1(2,3)=60;
// matrix1(2,4)=70;
// matrix1(3,5)=80;

// matrix2(0,0)=1;
// matrix2(2,0)=2;
// matrix2(0,1)=4;
// matrix2(1,2)=3;

// matrix1.printMatrix();
// matrix2.printMatrix();
// matrix1.compress();
// matrix2.compress();
// matrix1.printMatrix();
// matrix2.printMatrix();
// matrix1.uncompress();
// matrix2.uncompress();
// matrix1.printMatrix();
// matrix2.printMatrix();
// Generate a vector of the right dimension
    unsigned int seed = 1234;
    std::default_random_engine generator{seed};
    std::uniform_real_distribution<double> uniform(-10000., 10000.);
    std::vector<double> v(matrix3.get_cols());
    for (std::size_t i=0; i<matrix3.get_cols(); ++i)
        v[i] = uniform(generator);

    std::cout << "\n\nTESTING MATRIX - VECTOR PRODUCT" << std::endl;

    // result of the product
    std::vector<double> b(matrix3.get_rows());

    // TEST IN THE UNCOMPRESSED, ROW-WISE STORAGE
    clock.start();
    b = matrix3 * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with uncompressed, row-wise storage: "
              << clock << std::endl;

    // TEST IN THE COMPRESSED, ROW-WISE STORAGE
    matrix3.compress();
    clock.start();
    b = matrix3 * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with compressed, row-wise storage: "
              << clock << std::endl;


    // TEST IN THE UNCOMPRESSED, COLUMN-WISE STORAGE
    clock.start();
    b = matrix4 * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with uncompressed, column-wise storage: "
              << clock << std::endl;

    // TEST IN THE COMPRESSED, COLUMN-WISE STORAGE
    matrix4.compress();
    clock.start();
    b = matrix4 * v;
    clock.stop();
    std::cout << "Time for Matrix-Vecor multiplication with compressed, column-wise storage: "
              << clock << std::endl;

//matrix2.compress();
//std::vector<double> result4 = matrix2*vector;
//matrix2.printVec(result4);
//auto end = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

//std::cout << "Matrix-vector multiplication took " << duration << " microseconds." << std::endl;

return 0;
}


