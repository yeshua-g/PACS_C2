#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <map>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <numeric>
#include <cmath>



namespace algebra {

    // Enum class to specify the storage order of the matrix
    enum class StorageOrder {
        RowMajor,
        ColumnMajor
    };

    // Enum class to specify the norm type
    enum class Norm 
    {
        OneNorm,
        FrobeniusNorm,
        InfinityNorm
    };

    // Struct to define the comparison operator for the map
    template <StorageOrder Order>
    struct myOP{
        bool operator ()
        (const std::array<std::size_t,2> & lhs, const std::array<std::size_t,2> & rhs) const
        {
            if constexpr (Order == StorageOrder::ColumnMajor)
            return (lhs[1]<rhs[1] || (lhs[1]==rhs[1] && lhs[0]<rhs[0]));
            else 
            return (lhs < rhs);
        };
    };

    template<typename T, StorageOrder Order=StorageOrder::RowMajor>
    class Matrix {
    private:
        std::map<std::array<std::size_t, 2>, T, myOP<Order>> data; // Map for dynamic construction
        std::vector<std::size_t> innerIndex; // Inner indexes for compressed storage
        std::vector<std::size_t> outerIndex; // Outer indexes for compressed storage
        std::vector<T> values; // Values for compressed storage
        bool compressed; // Flag to indicate compression state
        // Variables to store the number of rows and columns
        std::size_t rows;
        std::size_t cols;
        
    public:
    
        //method do check range
        bool checkRange(std::size_t i, std::size_t j) const {
            return  i < rows && j < cols;
        }

        // Constructor
        Matrix(std::size_t rows = 1, std::size_t cols = 1) : compressed(false), rows(rows), cols(cols) {
            if(rows==0 || cols==0)
            throw std::invalid_argument("Invalid dimensions");
        }
        

        // Method to resize the matrix
        void resize(std::size_t n_rows, std::size_t n_cols) {
            if(compressed)
            uncompress();
            if(n_rows==0 || n_cols==0)
            throw std::invalid_argument("Invalid dimensions");
            this->rows = n_rows;
            this->cols = n_cols;
            for(auto it=data.begin();it!=data.end();){
                if(it->first[0]>=n_rows || it->first[1]>=n_cols)
                it=data.erase(it);
                else
                ++it;
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
            if constexpr(Order == StorageOrder::RowMajor) {
                std::size_t rowStart = innerIndex[i];
                std::size_t rowEnd = innerIndex[i + 1];
                // Search within this range in the innerIndex vector to find the column index j
                for (std::size_t k = rowStart; k < rowEnd; ++k) {
                if (outerIndex[k] == j) {
                    return values[k]; // Element found, return its value
                }
                }
                // Element not found, throw an error
                throw std::runtime_error("Element not found in compressed matrix.");
            } else {
                // Find the range of elements for column j
                std::size_t colStart = innerIndex[j];
                std::size_t colEnd = innerIndex[j + 1];
                // Search within this range in the innerIndex vector to find the row index i
                for (std::size_t k = colStart; k < colEnd; ++k) {
                if (outerIndex[k] == i) {
                    return values[k]; // Element found, return its value
                }
                }
                // Element not found, throw an error
                throw std::runtime_error("Element not found in compressed matrix.");
            }
            } else {
            // In uncompressed state, insert or update element in the data map
            return data[std::array<std::size_t, 2>{i, j}];
            }
        }
        
        // Method to convert to compressed format
        void compress() {
            // Implementation to convert data to compressed format (CSR or CSC)
            if(!compressed){
            if constexpr(Order==StorageOrder::ColumnMajor)
                CSC_compr();
            else if constexpr(Order==StorageOrder::RowMajor)
                CSR_compr();
            compressed = true; // Set compressed flag to true
            }
            else {
            std::cout << "Matrix is already compressed." << std::endl;
            return;
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
            innerIndex.push_back(0); // Start index of the first row
            
            for (const auto& entry : data) {
                std::array<std::size_t, 2> indices = entry.first;
                T value = entry.second;
                
                std::size_t row = indices[0];
                std::size_t col = indices[1];
                
                if (row != currentRow) {
                    // Start a new row
                    currentRow = row;
                    innerIndex.push_back(nnz);
                }
                
                outerIndex.push_back(col);
                values.push_back(value);
                nnz++;
            }
            
            // Add the end index of the last row
            innerIndex.push_back(nnz);
            
            // Clear the data map
            data.clear();
        }

        //Method to compress matrix in CSC format
        void CSC_compr(){
            // Clear compressed storage vectors
            innerIndex.clear();
            outerIndex.clear();
            values.clear();
            
            // Iterate over the data map and populate the compressed storage vectors
            std::size_t nnz = 0; // Number of non-zero elements
            std::size_t currentCol = 0; // Current col index
            innerIndex.push_back(0); // Start index of the first col
            
            for (const auto& entry : data) {
                std::array<std::size_t, 2> indices = entry.first;
                T value = entry.second;
                
                std::size_t row = indices[0];
                std::size_t col = indices[1];
                
                if (col != currentCol) {
                    // Start a new col
                    currentCol = col;
                    innerIndex.push_back(nnz);
                }
                
                outerIndex.push_back(row);
                values.push_back(value);
                nnz++;
            }
            
            // Add the end index of the last col
            innerIndex.push_back(nnz);
            
            // Clear the data map
            data.clear();
        }

        // Method to bring back to uncompressed state
        void uncompress() {
            // Implementation to bring back to uncompressed state
            if(compressed){
            if constexpr(Order==StorageOrder::ColumnMajor)
                CSC_uncompr();
            else if constexpr(Order==StorageOrder::RowMajor)
                CSR_uncompr();
            compressed = false; // Set compressed flag to true
            innerIndex.clear();
            outerIndex.clear();
            values.clear();
            }
            else {
            std::cout << "Matrix is already uncompressed." << std::endl;
            }
            return;
        }

        //Method to uncompress matrix in CSR format
        void CSR_uncompr(){
            // Clear data map
            data.clear();
            
            // Iterate over the compressed storage vectors and populate the data map
            std::size_t nnz = 0; // Number of non-zero elements
            
            for (std::size_t i = 0; i < innerIndex.size() - 1; ++i) {
            std::size_t rowStart = innerIndex[i];
            std::size_t rowEnd = innerIndex[i + 1];
            
            for (std::size_t k = rowStart; k < rowEnd; ++k) {
                std::array<std::size_t, 2> indices = {i, outerIndex[k]};
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
            
            for (std::size_t j = 0; j < innerIndex.size() - 1; ++j) {
            std::size_t colStart = innerIndex[j];
            std::size_t colEnd = innerIndex[j + 1];
            
            for (std::size_t k = colStart; k < colEnd; ++k) {
                std::array<std::size_t, 2> indices = {outerIndex[k], j};
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
            if constexpr(Order == StorageOrder::RowMajor) {
                std::size_t rowStart = innerIndex[i];
                std::size_t rowEnd = innerIndex[i + 1];
                for(std::size_t k = rowStart; k < rowEnd; ++k) {
                if(outerIndex[k] == j) {
                    return values[k];
                }
                }
                static T zero{};
                return zero;
            } else {
                std::size_t colStart = innerIndex[j];
                std::size_t colEnd = innerIndex[j + 1];
                for(std::size_t k = colStart; k < colEnd; ++k) {
                if(outerIndex[k] == i) {
                    return values[k];
                }
                }
                static T zero{};
                return zero;
            }
            } else {
            // Implementation to return elements in uncompressed state
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

            if(vec.size() != mat.cols) 
            throw std::invalid_argument("Invalid dimensions for matrix-vector multiplication");

            std::vector<T> result(mat.rows, T{});
            
            if (mat.compressed) {
            if constexpr(Order == StorageOrder::RowMajor) {
                for (std::size_t i = 0; i < mat.rows; ++i) {
                std::size_t rowStart = mat.innerIndex[i];
                std::size_t rowEnd = mat.innerIndex[i + 1];
                for (std::size_t k = rowStart; k < rowEnd; ++k) {
                    std::size_t col = mat.outerIndex[k];
                    result[i] += mat.values[k] * vec[col];
                }
                }
            } else {
                for (std::size_t j = 0; j < mat.cols; ++j) {
                std::size_t colStart = mat.innerIndex[j];
                std::size_t colEnd = mat.innerIndex[j + 1];
                for (std::size_t k = colStart; k < colEnd; ++k) {
                    std::size_t row = mat.outerIndex[k];
                    result[row] += mat.values[k] * vec[j];
                }
                }
            }
            } else {
            if constexpr(Order == StorageOrder::RowMajor) {
                for (std::size_t i=0;i<mat.rows;i++) {
                std::array<std::size_t, 2> indices = {i, 0};
                auto lower=mat.data.lower_bound(indices);
                indices={i,mat.cols-1};
                auto upper=mat.data.upper_bound(indices);
                for (auto it=lower;it!=upper;++it) {
                    std::size_t col = it->first[1];
                    result[i] += it->second * vec[col];
                    }
                }
            }
             else {
                for (std::size_t j=0;j<mat.cols;j++) {
                std::array<std::size_t, 2> indices = {0, j};
                auto lower=mat.data.lower_bound(indices);
                indices={mat.rows-1,j};
                auto upper=mat.data.upper_bound(indices);
                for (auto it=lower;it!=upper;++it) {
                    std::size_t row = it->first[0];
                    result[row] += it->second * vec[j];
                    }
                }
            }
        }
            return result;
        }


        // Method to read a matrix from a file in Matrix Market format
        void readMatrix(const std::string& filename) {
            std::ifstream file(filename);
            if (!file.is_open()) {
                throw std::runtime_error("Failed to open file: " + filename);
            }

            std::string line;
            std::getline(file, line);
            if (line != "%%MatrixMarket matrix coordinate real general") {
                throw std::runtime_error("Invalid Matrix Market format");
            }

            std::getline(file, line);
            std::istringstream iss(line);
            std::size_t rows, cols, nnz;
            iss >> rows >> cols >> nnz;

            // Resize the matrix
            resize(rows, cols);

            // Read the matrix entries
            for (std::size_t i = 0; i < nnz; ++i) {
                std::getline(file, line);
                std::istringstream iss(line);
                std::size_t row, col;
                T value;
                iss >> row >> col >> value;

                // Insert the value into the matrix
               
                (*this)(row - 1, col - 1) = value;
        
            }

            file.close();
        }

        
        std::size_t get_cols() const {
            return cols;
        }

        std::size_t get_rows() const {
            return rows;
        }

        std::size_t get_nnz() const {
            if(compressed)
            return values.size();
            else
            return data.size();
        }

        // Method to print the matrix
        void printMatrix() const {
            if(compressed){
                std::cout << "Compressed matrix:" << std::endl;
                std::cout << "Rows: " << rows << std::endl;
                std::cout << "Cols: " << cols << std::endl;
                std::cout << "NNZ: " << get_nnz() << std::endl;
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
            else{
                std::cout << "Uncompressed matrix:" << std::endl;
                std::cout << "Rows: " << rows << std::endl;
                std::cout << "Cols: " << cols << std::endl;
                std::cout << "NNZ: " << get_nnz() << std::endl;
                for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t j = 0; j < cols; ++j) {
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

    // Friend operator for matrix-matrix(1 col) multiplication
        friend std::vector<T> operator*(const Matrix<T, Order>& mat, const Matrix<T, Order>& mat2) {

            std::size_t vec_size = mat2.get_rows();;
            if(vec_size != mat.cols || mat2.get_cols() != 1)
            throw std::invalid_argument("Invalid dimensions for matrix-vector multiplication");

            std::vector<T> result(mat.rows, T{});
            
            if (mat.compressed) {
            if constexpr(Order == StorageOrder::RowMajor) {
                for (std::size_t i = 0; i < mat.rows; ++i) {
                std::size_t rowStart = mat.innerIndex[i];
                std::size_t rowEnd = mat.innerIndex[i + 1];
                for (std::size_t k = rowStart; k < rowEnd; ++k) {
                    std::size_t col = mat.outerIndex[k];
                    result[i] += mat.values[k] * mat2(col, 0);
                }
                }
            } else {
                for (std::size_t j = 0; j < mat.cols; ++j) {
                std::size_t colStart = mat.innerIndex[j];
                std::size_t colEnd = mat.innerIndex[j + 1];
                for (std::size_t k = colStart; k < colEnd; ++k) {
                    std::size_t row = mat.outerIndex[k];
                    result[row] += mat.values[k] * mat2(j, 0);
                }
                }
            }
            } else {
            if constexpr(Order == StorageOrder::RowMajor) {
                for (std::size_t i=0;i<mat.rows;i++) {
                std::array<std::size_t, 2> indices = {i, 0};
                auto lower=mat.data.lower_bound(indices);
                indices={i,mat.cols-1};
                auto upper=mat.data.upper_bound(indices);
                for (auto it=lower;it!=upper;++it) {
                    std::size_t col = it->first[1];
                    result[i] += it->second * mat2(col, 0);
                    }
                }
            }
             else {
                for (std::size_t j=0;j<mat.cols;j++) {
                std::array<std::size_t, 2> indices = {0, j};
                auto lower=mat.data.lower_bound(indices);
                indices={mat.rows-1,j};
                auto upper=mat.data.upper_bound(indices);
                for (auto it=lower;it!=upper;++it) {
                    std::size_t row = it->first[0];
                    result[row] += it->second * mat2(j, 0);
                    }
                }
            }
        }
            return result;
        }

    // Method to calculate the norm of the matrix
    template <Norm n>
    T norm() const {
        if constexpr (n == Norm::OneNorm) {
            auto sum_abs = [](double value, const auto& p) { return std::abs(value) + (std::abs(p.second)); };
            std::vector<T> col_sums(cols, T{});

            if constexpr (Order == StorageOrder::RowMajor) {
                if (!compressed) {
                    for (const auto& val : data) {
                        col_sums[val.first[1]] += std::abs(val.second);
                    }
                } else {
                    for (std::size_t i = 0; i < outerIndex.size(); ++i) {
                        col_sums[outerIndex[i]] += std::abs(values[i]);
                    }
                }
            } else {
                if (!compressed) {
                    for (std::size_t c = 0; c < cols; ++c) {
                        col_sums[c] = std::accumulate(data.lower_bound({ 0, c }),
                            data.upper_bound({ rows - 1, c }),
                            0.,
                            sum_abs);
                    }
                } else {
                    for (std::size_t c = 0; c < cols; ++c) {
                        for (std::size_t i = innerIndex[c]; i < innerIndex[c + 1]; ++i) {
                            col_sums[c] += std::abs(values[i]);
                        }
                    }
                }
            }

            return *std::max_element(col_sums.begin(), col_sums.end());
        }

        if constexpr (n == Norm::InfinityNorm) {
            auto sum_abs = [](double value, const auto& p) { return std::abs(value) + (std::abs(p.second)); };
            std::vector<T> row_sums(rows, T{});

            if constexpr (Order == StorageOrder::RowMajor) {
                if (!compressed) {
                    for (std::size_t r = 0; r < rows; ++r) {
                        row_sums[r] = std::accumulate(data.lower_bound({ r, 0 }),
                            data.upper_bound({ r, cols - 1 }),
                            0.,
                            sum_abs);
                    }
                } else {
                    for (std::size_t r = 0; r < rows; ++r) {
                        for (std::size_t i = innerIndex[r]; i < innerIndex[r + 1]; ++i) {
                            row_sums[r] += std::abs(values[i]);
                        }
                    }
                }
            } else {
                if (!compressed) {
                    for (const auto& val : data) {
                        row_sums[val.first[0]] += std::abs(val.second);
                    }
                } else {
                    for (std::size_t i = 0; i < outerIndex.size(); ++i) {
                        row_sums[outerIndex[i]] += std::abs(values[i]);
                    }
                }
            }

            return *std::max_element(row_sums.begin(), row_sums.end());
        }

        if constexpr (n == Norm::FrobeniusNorm) {
            if (!compressed) {
                auto sum_square = [](double value, const auto& p) { return std::abs(value) + ((std::abs(p.second)) * (std::abs(p.second))); };
                return std::sqrt(std::accumulate(data.cbegin(),
                    data.cend(),
                    0.,
                    sum_square));
            } else {
                auto sum_square = [](double value, const auto& p) { return std::abs(value) + ((std::abs(p)) * (std::abs(p))); };
                return std::sqrt(std::accumulate(values.cbegin(),
                    values.cend(),
                    0.,
                    sum_square));
            }
        }
    }
    };

    
    
} // namespace algebra



    
#endif // MATRIX_HPP