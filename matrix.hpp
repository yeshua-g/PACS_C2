#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <map>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>



namespace algebra {

    enum class StorageOrder {
        RowMajor,
        ColumnMajor
    };

    template<typename T, StorageOrder Order>
    class Matrix {
    private:
        std::map<std::array<std::size_t, 2>, T> data; // Map for dynamic construction
        std::vector<std::size_t> innerIndex; // Inner indexes for compressed storage
        std::vector<std::size_t> outerIndex; // Outer indexes for compressed storage
        std::vector<T> values; // Values for compressed storage
        bool compressed; // Flag to indicate compression state
        std::size_t rows;
        std::size_t cols;

    public:

        //method do check range
        bool checkRange(std::size_t i, std::size_t j) const {
            return  i < rows && j < cols;
        }

        // Constructor
        Matrix(std::size_t rows = 0, std::size_t cols = 0) : compressed(false), rows(rows), cols(cols)  {
            if(Order==StorageOrder::ColumnMajor){
                std::swap(rows, cols);
            }
        }

        // Method to resize the matrix
        void resize(std::size_t n_rows, std::size_t n_cols) {
            if(compressed)
            uncompress();
            if(n_rows==0 || n_cols==0)
            throw std::invalid_argument("Invalid dimensions");
            this->rows = n_rows;
            this->cols = n_cols;
            for (auto it = data.begin(); it != data.end();) {
                std::array<std::size_t, 2> indices = it->first;
                if (indices[0] >= rows || indices[1] >= cols) {
                    it = data.erase(it);
                } else {
                    ++it;
                }
            }
        }

        // Method to insert or update elements in uncompressed state
        T& operator()(std::size_t i, std::size_t j) {
            // Implementation to insert or update elements in data map
            if(!checkRange(i, j)) {
            throw std::out_of_range("Index out of range");
            }
            if (compressed) {
            // In compressed state, only change the value of existing non-zero elements
            if(Order == StorageOrder::RowMajor) {
                std::size_t rowStart = outerIndex[i];
                std::size_t rowEnd = outerIndex[i + 1];
                // Search within this range in the innerIndex vector to find the column index j
                for (std::size_t k = rowStart; k < rowEnd; ++k) {
                if (innerIndex[k] == j) {
                    return values[k]; // Element found, return its value
                }
                }
                // Element not found, throw an error
                throw std::runtime_error("Element not found in compressed matrix.");
            } else {
                // Find the range of elements for column j
                std::size_t colStart = outerIndex[j];
                std::size_t colEnd = outerIndex[j + 1];
                // Search within this range in the innerIndex vector to find the row index i
                for (std::size_t k = colStart; k < colEnd; ++k) {
                if (innerIndex[k] == i) {
                    return values[k]; // Element found, return its value
                }
                }
                // Element not found, throw an error
                throw std::runtime_error("Element not found in compressed matrix.");
            }
            } else {
            // In uncompressed state, insert or update element in the data map
            if(Order==StorageOrder::ColumnMajor)
                std::swap(i, j);
            return data[std::array<std::size_t, 2>{i, j}];
            }
        }
        
        // Method to convert to compressed format
        void compress() {
            // Implementation to convert data to compressed format (CSR or CSC)
            if(!compressed){
            if(Order==StorageOrder::ColumnMajor)
                CSC_compr();
            else if(Order==StorageOrder::RowMajor)
                CSR_compr();
            compressed = true; // Set compressed flag to true
            }
            else {
            std::cout << "Matrix is already compressed." << std::endl;
            }
        }

        //Method to compress matrix in CSR format
        void CSR_compr(){
            // Clear compressed storage vectors
            innerIndex.clear();
            outerIndex.clear();
            values.clear();
            
            // Iterate over the data map and populate the compressed storage vectors
            std::size_t nnz = 0; // Number of non-zero elements
            std::size_t currentRow = 0; // Current row index
            outerIndex.push_back(0); // Start index of the first row
            
            for (const auto& entry : data) {
                std::array<std::size_t, 2> indices = entry.first;
                T value = entry.second;
                
                std::size_t row = indices[0];
                std::size_t col = indices[1];
                
                if (row != currentRow) {
                    // Start a new row
                    currentRow = row;
                    outerIndex.push_back(nnz);
                }
                
                innerIndex.push_back(col);
                values.push_back(value);
                nnz++;
            }
            
            // Add the end index of the last row
            outerIndex.push_back(nnz);
            
            // Clear the data map
            data.clear();
        }

        //Transpose the matrix in uncompressed format
        std::map<std::array<std::size_t, 2>, T> transpose() const{

            if(!compressed){
            

            std::map<std::array<std::size_t, 2>, T> transposedData;
    
            for (const auto& entry : data) {
            std::array<std::size_t, 2> indices = entry.first;
            T value = entry.second;
            
            std::size_t row = indices[0];
            std::size_t col = indices[1];
            
            std::array<std::size_t, 2> transposedIndices = {col, row};
            transposedData[transposedIndices] = value;
            }

            return transposedData;
            }
            else{
                std::cerr << "Matrix is compressed. Please uncompress the matrix first." << std::endl;
                return data;
            }
        }

        //Method to compress matrix in CSC format
        void CSC_compr(){

            data=transpose();

            CSR_compr();
        }

        // Method to bring back to uncompressed state
        void uncompress() {
            // Implementation to bring back to uncompressed state
             if(compressed){
            if(Order==StorageOrder::ColumnMajor)
                CSC_uncompr();
            else if(Order==StorageOrder::RowMajor)
                CSR_uncompr();
            compressed = false; // Set compressed flag to true
            }
            else {
            std::cout << "Matrix is already uncompressed." << std::endl;
            }
        }

        //Method to uncompress matrix in CSR format
        void CSR_uncompr(){
            // Clear data map
            data.clear();
            
            // Iterate over the compressed storage vectors and populate the data map
            std::size_t nnz = 0; // Number of non-zero elements
            
            for (std::size_t i = 0; i < outerIndex.size() - 1; ++i) {
            std::size_t rowStart = outerIndex[i];
            std::size_t rowEnd = outerIndex[i + 1];
            
            for (std::size_t k = rowStart; k < rowEnd; ++k) {
                std::array<std::size_t, 2> indices = {i, innerIndex[k]};
                data[indices] = values[k];
                nnz++;
            }
            }
        }

        //Method to uncompress matrix in CSC format
        void CSC_uncompr(){
            // Clear data map
            data.clear();
            
            // Iterate over the compressed storage vectors and populate the data map
            std::size_t nnz = 0; // Number of non-zero elements
            
            for (std::size_t i = 0; i < outerIndex.size() - 1; ++i) {
            std::size_t colStart = outerIndex[i];
            std::size_t colEnd = outerIndex[i + 1];
            
            for (std::size_t k = colStart; k < colEnd; ++k) {
                std::array<std::size_t, 2> indices = {innerIndex[k], i};
                data[indices] = values[k];
                nnz++;
            }
            }
        }

        // Method to check if the matrix is in compressed state
        bool isCompressed() const {
            return compressed;
        }

        // Call operator to return elements of the matrix
        const T& operator()(std::size_t i, std::size_t j) const {
            // Implementation to return elements of the matrix
            if(!checkRange(i, j)) {
            throw std::out_of_range("Index out of range");
            }
            if(compressed) {
            // Implementation to return elements in compressed state
            if(Order == StorageOrder::RowMajor) {
                std::size_t colStart = outerIndex[j];
                std::size_t colEnd = outerIndex[j + 1];
                for(std::size_t k = colStart; k < colEnd; ++k) {
                if(innerIndex[k] == i) {
                    return values[k];
                }
                }
                static T zero{};
                return zero;
            } else {
                std::size_t rowStart = outerIndex[i];
                std::size_t rowEnd = outerIndex[i + 1];
                for(std::size_t k = rowStart; k < rowEnd; ++k) {
                if(innerIndex[k] == j) {
                    return values[k];
                }
                }
                static T zero{};
                return zero;
            }
            } else {
            // Implementation to return elements in uncompressed state
            if(Order==StorageOrder::ColumnMajor)
                std::swap(i, j);
            auto it = data.find(std::array<std::size_t, 2>{i, j});
            if(it != data.end()) {
                return it->second;
            } else {
                static T zero{};
                return zero;
            }
            }
        }

        // Friend operator for matrix-vector multiplication
        friend std::vector<T> operator*(const Matrix<T, Order>& mat, const std::vector<T>& vec) {
            std::vector<T> result(mat.rows, T{});
            
            if (mat.compressed) {
            if (Order == StorageOrder::RowMajor) {
                for (std::size_t i = 0; i < mat.rows; ++i) {
                std::size_t rowStart = mat.outerIndex[i];
                std::size_t rowEnd = mat.outerIndex[i + 1];
                for (std::size_t k = rowStart; k < rowEnd; ++k) {
                    std::size_t col = mat.innerIndex[k];
                    result[i] += mat.values[k] * vec[col];
                }
                }
            } else {
                for (std::size_t j = 0; j < mat.cols; ++j) {
                std::size_t colStart = mat.outerIndex[j];
                std::size_t colEnd = mat.outerIndex[j + 1];
                for (std::size_t k = colStart; k < colEnd; ++k) {
                    std::size_t row = mat.innerIndex[k];
                    result[row] += mat.values[k] * vec[j];
                }
                }
            }
            } else {
            for (const auto& entry : mat.data) {
                std::array<std::size_t, 2> indices = entry.first;
                T value = entry.second;
                std::size_t row = indices[0];
                std::size_t col = indices[1];
                result[row] += value * vec[col];
            }
            }
            
            return result;
        }


        // Method to read a matrix from a file in Matrix Market format
        

        // Method to return the number of columns
        std::size_t get_cols() const {
            return cols;
        }

        // Method to print the matrix
        void printMatrix() const {
            if(compressed){
                std::cout << "Outer indexes: ";
                for (std::size_t i = 0; i < outerIndex.size(); ++i) {
                    std::cout << outerIndex[i] << " ";
                }
                std::cout << std::endl;
                std::cout << "Inner indexes: ";
                for (std::size_t i = 0; i < innerIndex.size(); ++i) {
                    std::cout << innerIndex[i] << " ";
                }
                std::cout << std::endl;
                std::cout << "Values: ";
                for (const auto& value : values) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;
            }
            else if (Order == StorageOrder::RowMajor) {
                for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t j = 0; j < cols; ++j) {
                    std::cout << (*this)(i, j) << " ";
                }
                std::cout << std::endl;
                }
            } else if (Order == StorageOrder::ColumnMajor) {
                for (std::size_t j = 0; j < cols; ++j) {
                for (std::size_t i = 0; i < rows; ++i) {
                    std::cout << (*this)(i, j) << " ";
                }
                std::cout << std::endl;
                }
            }
             
        }

        //Method to print a vector
        void printVec(std::vector<T> vec) const{
            for (std::size_t i = 0; i < vec.size(); ++i) {
                std::cout << vec[i] << " ";
            }
            std::cout << std::endl;
        }
    };
} // namespace algebra

    
#endif // MATRIX_HPP