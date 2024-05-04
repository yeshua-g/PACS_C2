#include "matrix.hpp"

int main() {
    
algebra::Matrix<double,algebra::StorageOrder::RowMajor> matrix1(4,6);
algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> matrix2(3,3);

 
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

 //Testing the printMatrix function

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

 matrix2.resize(3,4);
 matrix2.printMatrix();

 matrix1.resize(4,3);
 matrix1.printMatrix();

std::vector<double> vec1 = {1,2,3,4};
std::vector<double> vec2(vec1.size());
std::vector<double> vec3 = {1,2,3};
std::vector<double> vec4(vec3.size());

vec2=matrix2*vec1;
matrix1.printVec(vec2);

vec4=matrix1*vec3;
matrix2.printVec(vec4);





    return 0;
}