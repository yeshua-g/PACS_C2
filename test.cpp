#include "matrix.hpp"

int main() {
    
    // Create two matrices: matrix1 with row-major storage order and matrix2 with column-major storage order
    algebra::Matrix<double,algebra::StorageOrder::RowMajor> matrix1(4,6);
    algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> matrix2(3,3);

    // Set values in matrix1
    matrix1(0,0)=10;
    matrix1(0,1)=20;
    matrix1(1,1)=30;
    matrix1(1,3)=40;
    matrix1(2,2)=50;
    matrix1(2,3)=60;
    matrix1(2,4)=70;
    matrix1(3,5)=80;

    // Set values in matrix2
    matrix2(0,0)=1;
    matrix2(2,0)=2;
    matrix2(0,1)=4;
    matrix2(1,2)=3;

    // Testing the printMatrix function
    std::cout << "Printing matrix1:" << std::endl;
    matrix1.printMatrix();
    std::cout << std::endl; // Add newline
    std::cout << "Printing matrix2:" << std::endl;
    matrix2.printMatrix();
    std::cout << std::endl; // Add newline

    // Compress the matrices
    matrix1.compress();
    matrix2.compress();

    // Print the compressed matrices
    std::cout << "Printing compressed matrix1:" << std::endl;
    matrix1.printMatrix();
    std::cout << std::endl; // Add newline
    std::cout << "Printing compressed matrix2:" << std::endl;
    matrix2.printMatrix();
    std::cout << std::endl; // Add newline

    // Uncompress the matrices
    matrix1.uncompress();
    matrix2.uncompress();

    // Print the uncompressed matrices
    std::cout << "Printing uncompressed matrix1:" << std::endl;
    matrix1.printMatrix();
    std::cout << std::endl; // Add newline
    std::cout << "Printing uncompressed matrix2:" << std::endl;
    matrix2.printMatrix();
    std::cout << std::endl; // Add newline

    // Resize matrix2 to 3x4
    matrix2.resize(3,4);
    std::cout << "Printing resized matrix2:" << std::endl;
    matrix2.printMatrix();
    std::cout << std::endl; // Add newline

    // Resize matrix1 to 4x3
    matrix1.resize(4,3);
    std::cout << "Printing resized matrix1:" << std::endl;
    matrix1.printMatrix();
    std::cout << std::endl; // Add newline

    // Create vectors vec1, vec2, vec3, vec4
    std::vector<double> vec1 = {1,2,3,4};
    std::vector<double> vec2(vec1.size());
    std::vector<double> vec3 = {1,2,3};
    std::vector<double> vec4(vec3.size());

    // Multiply matrix2 with vec1 and store the result in vec2
    vec2 = matrix2 * vec1;
    std::cout << "Printing matrix2*vec1:" << std::endl;
    for (const auto& value : vec2) {
        std::cout << value << " ";
    }
    std::cout << std::endl; // Add newline

    // Multiply matrix1 with vec3 and store the result in vec4
    vec4 = matrix1 * vec3;
    std::cout << "Printing matrix1*vec3:" << std::endl;
    for (const auto& value : vec4) {
        std::cout << value << " ";
    }
    std::cout << std::endl; // Add newline

    return 0;
}