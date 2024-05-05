#include <iostream>
#include <vector>
#include <complex>
#include "matrix.hpp"

int main() {
    // Create a 2x2 matrix with complex numbers
    algebra::Matrix<std::complex<double>, algebra::StorageOrder::RowMajor> mat(2, 2);
    mat(0, 0) = {1, 1};
    mat(0, 1) = {2, 2};
    mat(1, 0) = {3, 3};
    mat(1, 1) = {4, 4};

    // Print the matrix
    mat.printMatrix();
    std::cout << std::endl;

    // Create a vector with complex numbers
    std::vector<std::complex<double>> vec = {{1, 1}, {2, 2}};

    // Multiply the matrix and the vector
    std::vector<std::complex<double>> result = mat * vec;

    // Print the result
    std::cout << "Matrix times vector Result:" << std::endl;
    for (const auto& val : result) {
        std::cout << val << ' ';
    }
    std::cout << std::endl;

    return 0;
}