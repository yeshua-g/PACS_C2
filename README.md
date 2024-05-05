# Pacs Challenge 2 (Academic Year 2023-2024)

This project aims to create a matrix class capable of storing sparse matrices. The matrix class provides efficient storage and manipulation of sparse matrices, which are matrices that are mostly filled with zeros.

The matrix class allows for flexibility in the storage of the matrix. Users can choose between row-wise storage and column-wise storage depending on their needs. Additionally, the matrix class provides the option to compress the matrix in Compressed Sparse Column (CSC) or Compressed Sparse Row (CSR) format. This compression can significantly speed up matrix-vector multiplication, making the matrix class highly efficient for operations involving sparse matrices. This is the notation to instantiate an element of the class:


```
algebra::Matrix<double,algebra::StorageOrder::RowMajor> matrix1; //for a rowmajor ordering
algebra::Matrix<double,algebra::StorageOrder::ColMajor> matrix2; //for a colmajor ordering

matrix1.compress(); //CSR compressing technique
matrix2.compress(); //CSC compressing technique
```

To let the program work the User can just type `make all` in the terminal and this will generate four executables called **main**, **test**, **test2** and **extras**

Before starting it is !!!**EXTREMELY IMPORTANT**!!! to say that the code is written using the *chrono* utility, that is not included in this folder. So the User is required to change the first lines of the **Makefile** in order to set the proper path on his local machine.
```
# Lines the User should change for setting the right path on his machine 
# Compiler flags
CXXFLAGS = -O3 -Wall -I. -Wno-conversion-null -Wno-deprecated-declarations -I(//yourpathtopacs-examples/include) -std=c++20

# Linker flags
LDFLAGS = -L(//yourpathtopacs-examples/lib)  -Wl,-rpath=(//yourpathtopacs-examples/lib) 
```
## MAIN
the main is the requested code that tests the implementation and measure the time employed for the matrix-vector multiplication using the chrono utility
## TEST
test is an example with small matrixes to just visualize the implementation 
## TEST2
test 2 is to show that the matrix vector product works also with std::complex<T>
## EXTRAS
extras is for the matrix time one column matrix and norms. This is the norms notation:
```
matrix.norm<algebra::Norm::InfinityNorm>(); //infinity norm
matrix1.norm<algebra::Norm::OneNorm>(); // one norm
matrix1.norm<algebra::Norm::FrobeniusNorm>(); // Frobenius norm
```
