#include "dfsurface.h"
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::MatrixXd, Eigen::MatrixXcd, Eigen::VectorXd, Eigen::VectorXcd,
Eigen::Vector3d, Eigen::SparseMatrix, Eigen::SparseVector, Eigen::Triplet;

//output function for UV coordinate class
inline std::ostream &operator<<(std::ostream &os, SILT::PointUV &t) {
	return os<<t.u<<'\t'<<t.v;
}

//output function for Triangle  class
inline std::ostream &operator<<(std::ostream &os, SILT::Triangle &t) {
	return os<<t.a<<'\t'<<t.b<<'\t'<<t.c;
}

//output function for 3D point class
inline std::ostream &operator<<(std::ostream &os, SILT::Point3D &p) {
	return os<<p.x<<'\t'<<p.y<<'\t'<<p.z;
}

typedef std::complex<double> Complex;

SparseVector<Complex> CGMethod(SparseMatrix<Complex>& A, SparseVector<Complex>& b, SparseVector<Complex>& x0, double tol) {
    b = A.transpose() * b;
    A = A.transpose() * A;

    int maxIterations = 1000;
    
    Eigen::ConjugateGradient<SparseMatrix<Complex>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver initialization failed!" << std::endl;
        return x0;
    }

    solver.setMaxIterations(maxIterations);
    solver.setTolerance(tol);

    x0 = solver.solve(b);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver failed to converge!" << std::endl;
    } else {
        std::cout << "Converged after " << solver.iterations() << " iterations" << std::endl;
        std::cout << "Estimated error: " << solver.error() << std::endl;
    }

    return x0;
}

void fixedPoints(SILT::DFSurface surface) {
    // Vector3d vertexP1, vertexP2;
    double diameter, distance = (0,0);
    int index1, index2;

    // Indices of the two fixed vertices
    for (int i = 0; i < surface.vertices.size(); i++) {
        for (int j = i+1; j < surface.vertices.size(); j++) {
            distance = (surface.vertices[i].x - surface.vertices[j].x)*(surface.vertices[i].x - surface.vertices[j].x) + (surface.vertices[i].y - surface.vertices[j].y)*(surface.vertices[i].y - surface.vertices[j].y) + (surface.vertices[i].z - surface.vertices[j].z)*(surface.vertices[i].z - surface.vertices[j].z);
            if (distance > diameter) {
                diameter = distance;
                index1 = i;
                index2 = j;
                std::cout << "i = " << i << ", j = " << j << ", diameter = " << diameter << std::endl;
            }
        }
    }

    std::cout << "Diametro: " << sqrt(diameter) << std::endl;
    std::cout << "Indices: " << index1 << ", " << index2 << std::endl;
}

void process(SILT::DFSurface surface) {

    const int M = surface.triangles.size();
    const int N = surface.vertices.size();

    std::cout << "M = " << M << ", N = " << N << std::endl;

    SparseMatrix<Complex> M_c(M, N);
    SparseMatrix<Complex> M_cf(M, N-2);
    SparseMatrix<Complex> M_cp(M, 2);

    std::vector<Triplet<Complex>> triplets_c;
    triplets_c.reserve(M * 3); // Assuming each triangle contributes 3 non-zero elements

    VectorXcd u_p(4,1);
    SparseVector<Complex> u_f(2*(N-2),1);
    u_f.setZero();

    // Indices of the two fixed vertices
    Vector3d vertexP1, vertexP2;
    std::vector<int> p_indices = {41846, 60298};
    vertexP1 << surface.vertices[p_indices[0]].x, surface.vertices[p_indices[0]].y, surface.vertices[p_indices[0]].z;
    vertexP2 << surface.vertices[p_indices[1]].x, surface.vertices[p_indices[1]].y, surface.vertices[p_indices[1]].z;
    double diameter = (vertexP1 - vertexP2).norm();
    u_p << Complex(-diameter/2,0), Complex(diameter/2,0), Complex(0,0), Complex(0,0);

    for (int i = 0; i < M; ++i) {
        const auto& triangle = surface.triangles[i];
        Vector3d vertexA, vertexB, vertexC;
        Vector3d vecAlfa, vecBeta, vecGamma;
        int vertices_index[] = {triangle.a, triangle.b, triangle.c};

        vertexA << surface.vertices[triangle.a].x, surface.vertices[triangle.a].y, surface.vertices[triangle.a].z;
        vertexB << surface.vertices[triangle.b].x, surface.vertices[triangle.b].y, surface.vertices[triangle.b].z;
        vertexC << surface.vertices[triangle.c].x, surface.vertices[triangle.c].y, surface.vertices[triangle.c].z;

        // Vectores entre vertices
        vecAlfa = vertexB - vertexA;
        vecBeta = vertexC - vertexA;
        vecGamma = vertexC - vertexB;

        // Normas de los vectores
        double alfa = vecAlfa.norm();
        double beta = vecBeta.norm();
        double gamma = vecGamma.norm();

        double theta = acos(vecAlfa.dot(vecBeta)/(alfa*beta));
        
        // Mapeo a Base Local Ortonormal
        double mappedA[] = {0, 0};
        double mappedB[] = {alfa, 0};
        double mappedC[] = {beta*cos(theta), beta*sin(theta)};
        // valor absoluto del area
        double area = fabs(beta*alfa*sin(theta)/2);
        Complex value;

        // Calcular G[j,i] segun como viene en el paper
        for (int j = 0; j < 3; ++j) {
            if (j==0) {
                value = Complex(mappedC[0]-mappedB[0]/ sqrt(2*area), mappedC[1]-mappedB[1]/ sqrt(2*area));
            } else if (j==1) {
                value = Complex(mappedC[0]-mappedA[0] / sqrt(2*area), mappedC[1]-mappedA[1] / sqrt(2*area));
            } else {
                value = Complex(mappedB[0]-mappedA[0] / sqrt(2*area), mappedB[1]-mappedA[1] / sqrt(2*area));
            }
            triplets_c.push_back(Triplet<Complex>(i, vertices_index[j], value));
            // if (i % 10000 == 0) std::cout << "i = " << i << ", j = " << j << ", value = " << value << std::endl;
        }   
    }

    M_c.setFromTriplets(triplets_c.begin(), triplets_c.end());
    M_cf = (M_c.leftCols(p_indices[0])+M_c.middleCols(p_indices[0]+1, p_indices[1]-p_indices[0]-1)+M_c.rightCols(N-p_indices[1]-1)).pruned();
    M_cp = (M_c-M_cf).pruned();
    // std::cout << "Si lees esto, no se rompio" << std::endl;
    // Separate real and imaginary parts for M_cf
    SparseMatrix<Complex> M_cf_real(M_cf.rows(), M_cf.cols());
    SparseMatrix<Complex> M_cf_imag(M_cf.rows(), M_cf.cols());

    for (int k = 0; k < M_cf.outerSize(); ++k) {
        for (SparseMatrix<Complex>::InnerIterator it(M_cf, k); it; ++it) {
            M_cf_real.coeffRef(it.row(), it.col()) = Complex(it.value().real(), 0);
            M_cf_imag.coeffRef(it.row(), it.col()) = Complex(0, it.value().imag());
        }
    }
    // Separate real and imaginary parts for M_cp
    SparseMatrix<Complex> M_cp_real(M_cp.rows(), M_cp.cols());
    SparseMatrix<Complex> M_cp_imag(M_cp.rows(), M_cp.cols());

    for (int k = 0; k < M_cp.outerSize(); ++k) {
        for (SparseMatrix<Complex>::InnerIterator it(M_cp, k); it; ++it) {
            M_cp_real.coeffRef(it.row(), it.col()) = Complex(it.value().real(), 0);
            M_cp_imag.coeffRef(it.row(), it.col()) = Complex(0, it.value().imag());
        }
    }

    // std::cout << "M_cf_real" << std::endl;
    // for (int k = 0; k < 5; ++k) {
    //     for (SparseMatrix<Complex>::InnerIterator it(M_cf_real, k); it; ++it) {
    //         if (it.value() != Complex(0,0)) std::cout << it.value() << std::endl;
    //     }
    // }

    // std::cout << "M_cf_imag" << std::endl;
    // for (int k = 0; k < 5; ++k) {
    //     for (SparseMatrix<Complex>::InnerIterator it(M_cf_imag, k); it; ++it) {
    //         if (it.value() != Complex(0,0)) std::cout << it.value() << std::endl;
    //     }
    // }

    SparseMatrix<Complex> A(2*M,2*(N-2)), B(2*M,4);
    SparseVector<Complex> b(2*M);
    std::vector<Triplet<Complex>> triplets_A, triplets_B;
    
    // Populate triplets_A for matrix A
    for (int k = 0; k < M_cf_real.outerSize(); ++k) {
        for (SparseMatrix<Complex>::InnerIterator it(M_cf_real, k); it; ++it) {
            triplets_A.push_back(Triplet<Complex>(it.row(), it.col(), it.value()));
        }

        // if (k == M_cf_real.outerSize() - 1) {
        //     std::cout << "Listo" << std::endl;
        // }
    }

    for (int k = 0; k < M_cf_imag.outerSize(); ++k) {
        for (SparseMatrix<Complex>::InnerIterator it(M_cf_imag, k); it; ++it) {
            triplets_A.push_back(Triplet<Complex>(it.row(), it.col() + M_cf_real.cols(), it.value()));
        }

        // if (k == M_cf_imag.outerSize() - 1) {
        //     std::cout << "Listo" << std::endl;
        // }
    }

    // Populate triplets_B for matrix 
    for (int k = 0; k < M_cp_real.outerSize(); ++k) {
        for (SparseMatrix<Complex>::InnerIterator it(M_cp_real, k); it; ++it) {
            triplets_B.push_back(Triplet<Complex>(it.row(), it.col(), it.value()));
        }

        // if (k == M_cp_real.outerSize() - 1) {
        //     std::cout << "Listo" << std::endl;
        // }
    }

    for (int k = 0; k < M_cp_imag.outerSize(); ++k) {
        for (SparseMatrix<Complex>::InnerIterator it(M_cp_imag, k); it; ++it) {
            triplets_B.push_back(Triplet<Complex>(it.row(), it.col() + M_cp_real.cols(), it.value()));
        }

        // if (k == M_cp_imag.outerSize() - 1) {
        //     std::cout << "Listo" << std::endl;
        // }
    }

    A.setFromTriplets(triplets_A.begin(), triplets_A.end());
    B.setFromTriplets(triplets_B.begin(), triplets_B.end());

    // ver elementos de A
    // for (int k = 0; k < 5; ++k) {
    //     for (SparseMatrix<Complex>::InnerIterator it(A, k); it; ++it) {
    //         if (it.value() != Complex(0,0)) std::cout << it.value() << std::endl;
    //     }
    // }

    // for (int k = 0; k < 5; ++k) {
    //     for (SparseMatrix<Complex>::InnerIterator it(B, k); it; ++it) {
    //         if (it.value() != Complex(0,0)) std::cout << it.value() << std::endl;
    //     }
    // }

    b = -1 * (B * u_p).sparseView();
    SparseVector<Complex> mappedCoords(N);
    SparseVector<Complex> freePoints = CGMethod(A, b, u_f, 1e-3);

    mappedCoords.insert(p_indices[0]) = u_p[0]+u_p[2];
    mappedCoords.insert(p_indices[1]) = u_p[1]+u_p[3];

    int correction = 0;
    for (int i = 0; i < N; i++) {
        if (i == p_indices[0] || i == p_indices[1]) {
            correction += 1;
            continue;
        }
        // if (i%10==0) std::cout << "i = " << i << std::endl;
        mappedCoords.insert(i) = freePoints.coeff(i - correction);
    }

    for (int i = 0; i < M; i++) {
        const auto& triangle = surface.triangles[i];
        Vector3d vertexA, vertexB, vertexC;
    }

    SILT::DFSurface mappedBrain = surface;

   for (int i = 0; i < N; i++) {
        double X, Y;
        X = mappedCoords.coeff(i).real();
        Y = mappedCoords.coeff(i).imag();
        mappedBrain.vertices[i].x = diameter*X / (1 + X*X + Y*Y);
        mappedBrain.vertices[i].y = diameter*Y / (1 + X*X + Y*Y);
        mappedBrain.vertices[i].z = diameter*(-1 + X*X + Y*Y) / (2*(1 + X*X + Y*Y));
    }

    SparseMatrix<Complex> M_spring(3*M,N);
    std::vector<Triplet<Complex>> triplets_spring;
    triplets_spring.reserve(3*M); // Assuming each triangle contributes 3 non-zero elements

    for (int i=0; i < M; i++) {
        int edgeIndex = 0;
        const auto& triangle = surface.triangles[i];
        Vector3d vertexA, vertexB, vertexC;
        Vector3d oldA, oldB, oldC;
        Vector3d vecAlfa, vecBeta, vecGamma;
        Vector3d oldVecAlfa, oldVecBeta, oldVecGamma;

        vertexA << mappedBrain.vertices[triangle.a].x, mappedBrain.vertices[triangle.a].y, mappedBrain.vertices[triangle.a].z;
        vertexB << mappedBrain.vertices[triangle.b].x, mappedBrain.vertices[triangle.b].y, mappedBrain.vertices[triangle.b].z;
        vertexC << mappedBrain.vertices[triangle.c].x, mappedBrain.vertices[triangle.c].y, mappedBrain.vertices[triangle.c].z;
        
        vecAlfa = vertexB - vertexA;
        vecBeta = vertexC - vertexA;
        vecGamma = vertexC - vertexB;

        oldA << surface.vertices[triangle.a].x, surface.vertices[triangle.a].y, surface.vertices[triangle.a].z;
        oldB << surface.vertices[triangle.b].x, surface.vertices[triangle.b].y, surface.vertices[triangle.b].z;
        oldC << surface.vertices[triangle.c].x, surface.vertices[triangle.c].y, surface.vertices[triangle.c].z;

        oldVecAlfa = oldB - oldA;
        oldVecBeta = oldC - oldA;
        oldVecGamma = oldC - oldB;

        Complex edge1, edge2, edge3;

        edge1 = Complex(sqrt(vecAlfa.norm()/oldVecAlfa.norm()),0);
        edge2 = Complex(sqrt(vecBeta.norm()/oldVecBeta.norm()),0);
        edge3 = Complex(sqrt(vecGamma.norm()/oldVecGamma.norm()),0);

        triplets_spring.push_back(Triplet<Complex>(edgeIndex, triangle.a, edge1));
        triplets_spring.push_back(Triplet<Complex>(edgeIndex+1, triangle.b, edge2));
        triplets_spring.push_back(Triplet<Complex>(edgeIndex+2, triangle.c, edge3));
        edgeIndex += 3;

    }

    M_spring.setFromTriplets(triplets_spring.begin(), triplets_spring.end());

    SparseMatrix<Complex> matrix(4*M, N);
    std::vector<Triplet<Complex>> triplets;
    triplets.reserve(4*M); // Assuming each triangle contributes 3 non-zero elements

    for (int k = 0; k < M_c.outerSize(); ++k) {
        for (SparseMatrix<Complex>::InnerIterator it(M_c, k); it; ++it) {
            triplets.push_back(Triplet<Complex>(it.row(), it.col(), it.value()));
        }
    }

    // Triplets for M_spring
    for (int k = 0; k < M_spring.outerSize(); ++k) {
        for (SparseMatrix<Complex>::InnerIterator it(M_spring, k); it; ++it) {
            triplets.push_back(Triplet<Complex>(M + it.row(), it.col(), it.value()));
        }
    }

    matrix.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "matrix" << std::endl;

    SparseMatrix<Complex> matrix_real(matrix.rows(), matrix.cols());
    SparseMatrix<Complex> matrix_imag(matrix.rows(), matrix.cols());
    std::vector<Triplet<Complex>> triplets_real, triplets_imag;
    triplets_real.reserve(4*M); // Assuming each triangle contributes 3 non-zero elements
    triplets_imag.reserve(4*M); // Assuming each triangle contributes 3 non-zero elements

    for (int k = 0; k < matrix.outerSize(); ++k) {
        for (SparseMatrix<Complex>::InnerIterator it(matrix, k); it; ++it) {
            triplets_real.push_back(Triplet<Complex>(it.row(), it.col(), Complex(it.value().real(), 0)));
            triplets_imag.push_back(Triplet<Complex>(it.row(), it.col(), Complex(0, it.value().imag())));
        }
    }

    matrix_real.setFromTriplets(triplets_real.begin(), triplets_real.end());
    matrix_imag.setFromTriplets(triplets_imag.begin(), triplets_imag.end());


    //mappedBrain.writeDFS("leftMappedBrain.dfs");

    std::cout << "Codigo Completo" << std::endl;
}

int main() {
    SILT::DFSurface leftBrain, rightBrain;
    leftBrain.readDFS("atlas.left.mid.cortex.svreg.dfs");
    //rightBrain.readDFS("atlas.right.mid.cortex.svreg.dfs");

    //fixedPoints(leftBrain);
    process(leftBrain);
    //process(rightBrain);

}