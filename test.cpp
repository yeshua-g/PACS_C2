#include "matrix.hpp"
#include "vector"
#include <chrono>

int main() {
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
std::vector<double> vector(matrix3.get_cols());
for (std::size_t i = 0; i < matrix3.get_cols(); i++) {
    vector[i] = i;
}

// Perform matrix-vector multiplication and time the executionon
auto start = std::chrono::high_resolution_clock::now();
std::vector<double> result1 = matrix3 * vector;
auto end = std::chrono::high_resolution_clock::now();
auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

start = std::chrono::high_resolution_clock::now();
std::vector<double> result2 = matrix4 * vector;
end = std::chrono::high_resolution_clock::now();
auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

matrix3.compress();
start = std::chrono::high_resolution_clock::now();
std::vector<double> result3 = matrix3 * vector;
end = std::chrono::high_resolution_clock::now();
auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

matrix4.compress();
start = std::chrono::high_resolution_clock::now();
std::vector<double> result4 = matrix4 * vector;
end = std::chrono::high_resolution_clock::now();
auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

std::cout << "Matrix-vector multiplication 1 took " << duration1 << " microseconds." << std::endl;
std::cout << "Matrix-vector multiplication 2 took " << duration2 << " microseconds." << std::endl;
std::cout << "Matrix-vector multiplication 3 took " << duration3 << " microseconds." << std::endl;
std::cout << "Matrix-vector multiplication 4 took " << duration4 << " microseconds." << std::endl;
//matrix2.compress();
//std::vector<double> result4 = matrix2*vector;
//matrix2.printVec(result4);
//auto end = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

//std::cout << "Matrix-vector multiplication took " << duration << " microseconds." << std::endl;

if (result1 == result2 && result2 == result3 && result3 == result4) {
    std::cout << "All four results are equal." << std::endl;
} else {
    std::cout << "The results are not equal." << std::endl;
}
return 0;
}


