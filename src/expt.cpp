#include <iostream>
#include <memory>
#include <csignal>

template<typename T>
struct SharedMatrix {
    std::shared_ptr<T> backing;
    const size_t cols, rows;

    SharedMatrix(size_t cols, size_t rows) :
        // backing(std::make_shared<T[]>(cols * rows))
        backing(new T[cols * rows]),
        cols(cols), rows(rows)
    {}

    T* operator[](const size_t row) const {
        return backing.get() + row*cols;
    }
};

template<typename T>
class Matrix {
    SharedMatrix<T> storage;
    const size_t off_col, off_row;
public:
    const size_t cols, rows;

    Matrix(const Matrix<T> &other, size_t off_col, size_t off_row, size_t cols, size_t rows) :
        storage(other.storage),
        off_col(other.off_col + off_col),
        off_row(other.off_row + off_row),
        cols(cols), rows(rows)
    {
        assert(off_col + cols <= storage.cols);
        assert(off_row + rows <= storage.rows);
    }

    Matrix(size_t cols, size_t rows) :
        storage(cols, rows),
        off_col(0), off_row(0),
        cols(cols), rows(rows)
    {}

    Matrix(size_t dimension) :
        Matrix(dimension, dimension)
    {}

    Matrix<T> submatrix(size_t off_col, size_t off_row, size_t cols, size_t rows) {
        return Matrix(this, off_col, off_row, cols, rows);
    }

    Matrix<T> submatrix(size_t inset) const {
        return SubMatrix(1, 1, cols-1, rows-1);
    }

    T* operator[](const size_t row) const {
        return storage[row + off_row] + off_col;
    }

    T sum() const {
        T ret = 0;
        size_t row, col;
        for (row = 0; row < rows; row++)
            for (col = 0; col < cols; col++)
                ret += (*this)[row][col];
        return ret;
    }

    void fill(const T value) {
        size_t row, col;
        for (row = 0; row < rows; row++)
            for (col = 0; col < cols; col++)
                (*this)[row][col] = value;
    }

    void zero() {
        fill(0);
    }

    void print(const T thresh) const {
        size_t row, col;
        for (row = 0; row < rows; row++) {
            for (col = 0; col < cols; col++) {
                T val = (*this)[row][col];
                if (val*val <= thresh*thresh)
                    printf("             ");
                else
                    printf("% 0.5e ", val);
            }
            std::cout << std::endl;
        }
    }

    void print2() const {
        size_t row, col;
        std::ios_base::fmtflags flags( std::cout.flags() );
        std::cout << std::hexfloat;

        for (row = 0; row < rows; row++) {
            for (col = 0; col < cols; col++)
                std::cout << (*this)[row][col] << " ";
            std::cout << std::endl;
        }

        std::cout.flags(flags);
    }

    void print() const {
        print(0);
    }
};

class Lattice {
    double bias;
    size_t width, distance, dimension;
    Matrix<double> mat1_, mat2_; // padded backing
    Matrix<double> mat1,  mat2;  // 'real' view
    bool phase;

public:

    Lattice(double bias, size_t width, size_t distance, size_t dimension) :
        bias(bias), width(width), distance(distance), dimension(dimension),
        mat1_(dimension+2), mat2_(dimension+2),
        mat1(mat1_, 1, 1, dimension, dimension),
        mat2(mat2_, 1, 1, dimension, dimension),
        phase(true)
    {
        mat1_.zero();
        mat2_.zero();
        assert(width + distance < dimension);
    }

    double p() const {
        return 0.5 * (1 + bias);
    }

    double q() const {
        return 0.5 * (1 - bias);
    }

    void evolve() {
        Matrix<double> &m1 = phase ? mat1 : mat2,
                       // &m1_= phase ? mat1_: mat2_,
                       &m2 = phase ? mat2 : mat1,
                       &m2_= phase ? mat2_: mat1_;

        m2_.zero();
        double p2 = p() * 0.5;
        double q2 = q() * 0.5;
        double cur = 0;
        ssize_t x, y;
        size_t row, z = dimension - 1;

        for (y = 0; y < dimension; y++) {
            for (x = 0; x < dimension; x++) {
                cur = m1[y][x];
                m2[y-1][x] += p2 * cur;
                m2[y][x-1] += p2 * cur;
                m2[y+1][x] += q2 * cur;
                m2[y][x+1] += q2 * cur;
            }
        }

        // reflecting boundary
        for (x = 0; x < dimension; x++) {
            m2[x][0] += m2[x][ -1];
            m2[x][z] += m2[x][z+1];
            m2[0][x] += m2[ -1][x];
            m2[z][x] += m2[z+1][x];

            m2[x][ -1] = 0;
            m2[x][z+1] = 0;
            m2[ -1][x] = 0;
            m2[z+1][x] = 0;
        }

        // inject
        // row = width + distance;
        // x = row - 1;
        // y = 0;
        // double val = (double)1.0 / (double)row;
        // while (x >= 0) m2[y++][x--] += val;

        // inject
        row = width + distance;
        x = row - 1;
        y = 0;
        double val = (double)1.0 / (double)(row - 1);
                       m2[y++][x--] += val * 0.5;
        while (x >= 1) m2[y++][x--] += val;
                       m2[y++][x--] += val * 0.5;

        // absorb
        row = width;
        x = row - 1;
        y = 0;
        while (x >= 0) m2[y++][x--] = 0;

        phase = !phase;
    }

    double error() const {
        double ret = 0, cur;
        for (size_t y = 0; y < dimension; y++) {
            for (size_t x = 0; x < dimension; x++) {
                cur = mat2[y][x] - mat1[y][x];
                ret += cur * cur;
            }
        }
        return ret;
    }

    double mfpt() const {
        return (phase ? mat1 : mat2).sum();
    }

    void print() const {
        (phase ? mat1 : mat2).print2();
    }

    void print_currents() const {
        const Matrix<double> &mat = (phase ? mat1 : mat2);
        Matrix<double> sx(dimension), sy(dimension), divs(dimension-1);
        sx.zero(); sy.zero(); divs.zero();

        double p2 = p() * 0.5;
        double q2 = q() * 0.5;
        for (size_t y = 0; y < dimension; y++) {
            for (size_t x = 0; x < dimension; x++) {
                sx[y][x] = p2*mat[y][x+1] - q2*mat[y][x];
                sy[y][x] = p2*mat[y+1][x] - q2*mat[y][x];
            }
        }
        for (size_t y = 1; y < dimension; y++)
            for (size_t x = 1; x < dimension; x++)
                divs[y-1][x-1] = sx[y][x]-sx[y][x-1]
                               + sy[y][x]-sy[y-1][x];

        std::cout << "Current (left)" << std::endl;
        sx.print(1e-15);
        std::cout << std::endl;

        std::cout << "Current (up)" << std::endl;
        sy.print(1e-15);
        std::cout << std::endl;

        std::cout << "Current Divergence" << std::endl;
        divs.print(1e-15);
        std::cout << std::endl;
    }
};

namespace { volatile std::atomic_bool interrupted; }
void interrupt_handler(int signal) { interrupted = true; }

int main(int argc, char **argv) {
    std::signal(SIGINT, interrupt_handler);

    double bias = 1;
    size_t width = 1;
    size_t distance = 10;
    size_t dimension = 20;

    std::cerr << "Usage: " << argv[0] << " [bias [width [dist [sim_size]]]]" << std::endl;
    std::cerr << "       ctrl-c to stop" << std::endl << std::endl;

    if (argc > 1) bias = atof(argv[1]);
    if (argc > 2) width = atoll(argv[2]);
    if (argc > 3) distance = atoll(argv[3]);
    if (argc > 4) dimension = atoll(argv[4]);

    Lattice rw(bias, width, distance, dimension);
    int converged = 0;
    for (int i = 1; !interrupted; i++) {
        rw.evolve();
        if (i%100 == 0) {
            double err = rw.error();
            std::cerr << rw.mfpt() << " (" << err << ")" << std::endl;
            if (err == 0) {
                converged++;
                if (converged > 10)
                    break;
            } else converged = 0;
        }
    }
    rw.print();
    // std::cout << std::endl;
    // rw.print_currents();
    return 0;
}