#include <iostream>
#include <Eigen/Sparse>

int main() {
    // Create a sparse vector (sparse matrix with a single column)
    Eigen::SparseMatrix<double> sparseVector(10, 1);

    // Insert a value at position 2
    std::cout << "Original dimensions:\n" << sparseVector.rows() << std::endl;
    sparseVector.insert(2, 0) = 42.0;

    // Display the result
    std::cout << "Sparse Vector:\n" << sparseVector << std::endl;
    std::cout << "Dimensions:\n" << sparseVector.rows() << std::endl;

    // Dynamically resize the sparse vector
    sparseVector.conservativeResize(11, 1);

    // Display the resized result
    std::cout << "Resized Sparse Vector:\n" << sparseVector << std::endl;
    std::cout << "Resized Dimensions:\n" << sparseVector.rows() << std::endl;

    return 0;
}
