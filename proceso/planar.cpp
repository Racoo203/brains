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
    std::vector<VectorXcd> x;
    std::vector<VectorXcd> r;
    std::vector<VectorXcd> p;
    std::vector<double> alfa, beta;

    std::vector<Complex> finalCoords;

    // dimensiones de todo

    // std::cout << "A dimensiones: " << A.rows() << "x" << A.cols() << std::endl;
    // std::cout << "b dimensiones: " << b.rows() << "x" << b.cols() << std::endl;
    // std::cout << "x0 dimensiones: " << x0.rows() << "x" << x0.cols() << std::endl;
    int maxIterations = 10;
    b = A.transpose() * b;
    // std::cout << "Nueva b dimensiones: " << b.rows() << "x" << b.cols() << std::endl;

    A = A.transpose() * A;
    // std::cout << "Nueva A dimensiones: " << A.rows() << "x" << A.cols() << std::endl;
    

    x.push_back(x0); // x0
    r.push_back(b - A * x0); // r0

    if (r[0].norm() < tol) {
        std::cout << "Converged after 0 iterations" << std::endl;
        return x0;
    }

    p.push_back(r[0]); // p0

    for (int k = 0; k < maxIterations; k++) {
        // alfa 0
        alfa.push_back((r[k].transpose()*r[k]).value().real() / (p[k].transpose()*A*p[k]).value().real());
        // x 1
        x.push_back(x[k] + alfa[k] * p[k]);
        r.push_back(r[k] - alfa[k] * A * p[k]);

        if (r[k+1].norm() < tol) {
            std::cout << "Converged after " << k+1 << " iterations" << std::endl;
            break;
        }

        beta.push_back((r[k+1].transpose()*r[k+1]).value().real() / (r[k].transpose()*r[k]).value().real());
        p.push_back(r[k+1] + beta[k] * p[k]);

        if (k % 10 == 0) {
            std::cout << "k = " << k << " res = " << r.back().norm() << std::endl;
        }
    }

    // for (int i = 0; i < x.back().size(); i++) {
    //     if (x.back()[i] != Complex(0,0)) std::cout << x.back()[i] << std::endl;
    // }

    int splitIndex = x.back().size() / 2;

    // Add each part separately
    VectorXcd result = x.back().head(splitIndex) + x.back().tail(splitIndex);
    int zeroCounter = 0;
    for (int i = 0; i < result.size(); i++) {
        if (result[i] != Complex(0,0)) {
            std::cout << result[i] << std::endl;
        } else zeroCounter++;
    }

    std::cout << "X dimensiones: " << x.back().rows() << "x" << x.back().cols() << std::endl;
    std::cout << "R norma: " << r.back().norm() << std::endl;
    std::cout << "Result dimensiones: " << result.rows() << "x" << result.cols() << std::endl;
    std::cout << "Number of zeros: " << zeroCounter << std::endl;
    // llenar vector finalCoords con los valores de result y con las coordenadas fijas en sus respectivos indices

    return result.sparseView();
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
    u_p << Complex(-diameter/2 * 10,0), Complex(diameter/2 * 10,0), Complex(0,0), Complex(0,0);

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
            if (i % 10000 == 0) std::cout << "i = " << i << ", j = " << j << ", value = " << value << std::endl;
        }   
    }

    M_c.setFromTriplets(triplets_c.begin(), triplets_c.end());

    // Extract M_cf by removing the columns corresponding to the fixed vertices
    M_cf = M_c.block(0, 0, M, p_indices[0]).pruned();
    M_cf.conservativeResize(M, N - 2);

    // Extract M_cp by selecting only the columns corresponding to the fixed vertices
    M_cp = M_c.block(0, p_indices[0], M, 2).pruned();

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
    
    mappedCoords = freePoints.head(N-2);

    std::cout << "mappedCoords dimensiones: " << mappedCoords.rows() << "x" << mappedCoords.cols() << std::endl;
    std::cout << "numero de vertices " << surface.vertices.size() << std::endl;
    std::cout << "Codigo Completo" << std::endl;
}

int main() {
    SILT::DFSurface leftBrain, rightBrain;
    leftBrain.readDFS("atlas.left.mid.cortex.svreg.dfs");
    rightBrain.readDFS("atlas.right.mid.cortex.svreg.dfs");

    //fixedPoints(leftBrain);
    process(leftBrain);
    //process(rightBrain);

}