#include "matrix.hpp"
#include "vector"

int main() {
    // Create an instance of your matrix class
algebra::Matrix<double,algebra::StorageOrder::RowMajor> matrix1(4,6);
algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> matrix2(3,3);

// Read the matrix from the file
matrix1(0,0)=10;
matrix1(0,1)=20;
matrix1(1,1)=30;
matrix1(1,3)=40;
matrix1(2,2)=50;
matrix1(2,3)=60;
matrix1(2,4)=70;
matrix1(3,5)=80;

matrix2(0,0)=1;
matrix2(2,0)=2;
matrix2(0,1)=4;
matrix2(1,2)=3;

matrix1.printMatrix();
matrix2.printMatrix();
matrix1.compress();
matrix2.compress();
matrix1.printMatrix();
matrix2.printMatrix();
matrix1.uncompress();
matrix2.uncompress();
matrix1.printMatrix();
matrix2.printMatrix();
// Generate a vector of the right dimension
std::vector<double> vector(matrix1.get_cols());
for (std::size_t i = 0; i < matrix1.get_cols(); i++) {
    vector[i] = i;
}


// Perform matrix-vector multiplication and time the execution
//auto start = std::chrono::high_resolution_clock::now();
std::vector<double> result = matrix1*vector;
matrix1.printVec(result);
matrix1.compress();
std::vector<double> result2 = matrix1*vector;
matrix1.printVec(result2);

std::vector<double> result3 = matrix2*vector;
matrix2.printVec(result3);
matrix2.compress();
std::vector<double> result4 = matrix2*vector;
matrix2.printVec(result4);
//auto end = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

//std::cout << "Matrix-vector multiplication took " << duration << " microseconds." << std::endl;

return 0;
}


