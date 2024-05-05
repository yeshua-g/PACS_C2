# Pacs Challenge 2 (Academic Year 2023-2024)

This README.md is written in Markdown so I suggest to read it [online](https://github.com/yeshua-g/PACS_C2/blob/main/README.md).

This project aims to create a matrix class capable of storing sparse matrices. The matrix class provides efficient storage and manipulation of sparse matrices, which are matrices that are mostly filled with zeros.

The matrix class allows for flexibility in the storage of the matrix. Users can choose between row-wise storage and column-wise storage depending on their needs. Additionally, the matrix class provides the option to compress the matrix in Compressed Sparse Column (CSC) or Compressed Sparse Row (CSR) format. This compression can significantly speed up matrix-vector multiplication, making the matrix class highly efficient for operations involving sparse matrices.

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
the main is the requested code that tests the implementation
## TEST
test is an example with small matrixes to visualize the implementation 
## TEST2
test 2 is to show that the matrix vector product works also with std::complex<T>
## EXTRAS
extras is for the matrix time one column matrix and norms
