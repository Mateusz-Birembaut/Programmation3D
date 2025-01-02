#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <immintrin.h> // SIMD

class Vec3;
static inline Vec3 operator + (Vec3 const & a , Vec3 const & b);
static inline Vec3 operator * (float a , Vec3 const & b);

class Vec3 {
private:
    union {
        struct { float x, y, z; };
        __m128 simd; // Utilisation de SIMD pour les opérations vectorielles
    };
public:
    Vec3() : simd(_mm_setzero_ps()) {}
    Vec3(float x, float y, float z) : simd(_mm_set_ps(0.0f, z, y, x)) {}
    Vec3(__m128 simd) : simd(simd) {}

    __m128 getSimd() const { return simd; }

    float & operator [] (unsigned int c) {
        switch(c) {
            case 0: return x;
            case 1: return y;
            default: return z;
        }
    }
    float operator [] (unsigned int c) const {
        switch(c) {
            case 0: return x;
            case 1: return y;
            default: return z;
        }
    }

    Vec3 operator = (Vec3 const & other) {
       simd = other.simd;
       return *this;
    }

    float squareLength() const {
        return x*x + y*y + z*z;
    }

    float length() const { return sqrt(squareLength()); }

    inline float norm() const { return length(); }
    inline float squareNorm() const { return squareLength(); }

    void normalize() {
        float L = length();
        simd = _mm_div_ps(simd, _mm_set1_ps(L));
    }

    static float dot(Vec3 const & a, Vec3 const & b) {
        __m128 mul = _mm_mul_ps(a.simd, b.simd);
        float result[4];
        _mm_storeu_ps(result, mul);
        return result[0] + result[1] + result[2];
    }

    static Vec3 cross(Vec3 const & a, Vec3 const & b) {
        return Vec3(
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        );
    }

    void operator += (Vec3 const & other) {
        simd = _mm_add_ps(simd, other.simd);
    }

    void operator -= (Vec3 const & other) {
        simd = _mm_sub_ps(simd, other.simd);
    }

    void operator *= (float s) {
        simd = _mm_mul_ps(simd, _mm_set1_ps(s));
    }

    void operator /= (float s) {
        simd = _mm_div_ps(simd, _mm_set1_ps(s));
    }

    static Vec3 compProduct(Vec3 const & a, Vec3 const & b) {
        return Vec3(_mm_mul_ps(a.simd, b.simd));
    }

    unsigned int getMaxAbsoluteComponent() const {
        if (fabs(x) > fabs(y)) {
            if (fabs(x) > fabs(z)) {
                return 0;
            }
            return 2;
        }
        if (fabs(y) > fabs(z)) {
            return 1;
        }
        return 2;
    }

    Vec3 getOrthogonal() const {
        unsigned int c1 = getMaxAbsoluteComponent();
        unsigned int c2 = (c1 + 1) % 3;
        Vec3 res(0, 0, 0);
        res[c1] = z;
        res[c2] = -y;
        return res;
    }

};

static inline Vec3 operator + (Vec3 const & a, Vec3 const & b) {
    return Vec3(_mm_add_ps(a.getSimd(), b.getSimd()));
}

static inline Vec3 operator - (Vec3 const & a, Vec3 const & b) {
    return Vec3(_mm_sub_ps(a.getSimd(), b.getSimd()));
}

static inline Vec3 operator * (float a, Vec3 const & b) {
    return Vec3(_mm_mul_ps(_mm_set1_ps(a), b.getSimd()));
}

static inline Vec3 operator * (Vec3 const & b, float a) {
    return Vec3(_mm_mul_ps(_mm_set1_ps(a), b.getSimd()));
}

static inline Vec3 operator / (Vec3 const & a, float b) {
    return Vec3(_mm_div_ps(a.getSimd(), _mm_set1_ps(b)));
}

static inline std::ostream & operator << (std::ostream & s, Vec3 const & p) {
    s << p[0] << " " << p[1] << " " << p[2];
    return s;
}

static inline std::istream & operator >> (std::istream & s, Vec3 & p) {
    s >> p[0] >> p[1] >> p[2];
    return s;
}

static inline Vec3 operator * (Vec3 const & b, Vec3 const & a) {
    return Vec3(_mm_mul_ps(a.getSimd(), b.getSimd()));
}

static inline bool operator == (Vec3 const & b, Vec3 const & a) {
    return (b[0] == a[0] && b[1] == a[1] && b[2] == a[2]);
}


class Mat3 {
public:
    ////////////         CONSTRUCTORS          //////////////
    Mat3() {
        vals[0] = 0;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = 0;
        vals[4] = 0;
        vals[5] = 0;
        vals[6] = 0;
        vals[7] = 0;
        vals[8] = 0;
    }

    Mat3(float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8, float v9) {
        vals[0] = v1;
        vals[1] = v2;
        vals[2] = v3;
        vals[3] = v4;
        vals[4] = v5;
        vals[5] = v6;
        vals[6] = v7;
        vals[7] = v8;
        vals[8] = v9;
    }

    Mat3(const Mat3 &m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                (*this)(i, j) = m(i, j);
    }

    Mat3(const Vec3& col1, const Vec3& col2, const Vec3& col3) {
        vals[0] = col1[0]; vals[1] = col1[1]; vals[2] = col1[2];
        vals[3] = col2[0]; vals[4] = col2[1]; vals[5] = col2[2];
        vals[6] = col3[0]; vals[7] = col3[1]; vals[8] = col3[2];
    }

    // Multiplication de matrice avec un Vec3 : m.p
    //--> application d'un matrice de rotation à un point ou un vecteur
    Vec3 operator*(const Vec3 &p) {
        //Pour acceder a un element de la matrice (*this)(i,j) et du point p[i]
        Vec3 res = Vec3(
                    (*this)(0, 0) * p[0] + (*this)(0, 1) * p[1] + (*this)(0, 2) * p[2],
                (*this)(1, 0) * p[0] + (*this)(1, 1) * p[1] + (*this)(1, 2) * p[2],
                (*this)(2, 0) * p[0] + (*this)(2, 1) * p[1] + (*this)(2, 2) * p[2]);
        return res;
    }

    Mat3 operator*(const Mat3 &m2) { // calcul du produit matriciel m1.m2
        //Pour acceder a un element de la premiere matrice (*this)(i,j) et de la deuxième m2(k,l)
        Mat3 res = Mat3(
                    (*this)(0, 0) * m2(0, 0) + (*this)(0, 1) * m2(1, 0) + (*this)(0, 2) * m2(2, 0),
                    (*this)(0, 0) * m2(0, 1) + (*this)(0, 1) * m2(1, 1) + (*this)(0, 2) * m2(2, 1),
                    (*this)(0, 0) * m2(0, 2) + (*this)(0, 1) * m2(1, 2) + (*this)(0, 2) * m2(2, 2),
                    (*this)(1, 0) * m2(0, 0) + (*this)(1, 1) * m2(1, 0) + (*this)(1, 2) * m2(2, 0),
                    (*this)(1, 0) * m2(0, 1) + (*this)(1, 1) * m2(1, 1) + (*this)(1, 2) * m2(2, 1),
                    (*this)(1, 0) * m2(0, 2) + (*this)(1, 1) * m2(1, 2) + (*this)(1, 2) * m2(2, 2),
                    (*this)(2, 0) * m2(0, 0) + (*this)(2, 1) * m2(1, 0) + (*this)(2, 2) * m2(2, 0),
                    (*this)(2, 0) * m2(0, 1) + (*this)(2, 1) * m2(1, 1) + (*this)(2, 2) * m2(2, 1),
                    (*this)(2, 0) * m2(0, 2) + (*this)(2, 1) * m2(1, 2) + (*this)(2, 2) * m2(2, 2)
                    );
        return res;
    }

    bool isnan() const {
        return std::isnan(vals[0]) || std::isnan(vals[1]) || std::isnan(vals[2])
                || std::isnan(vals[3]) || std::isnan(vals[4]) || std::isnan(vals[5])
                || std::isnan(vals[6]) || std::isnan(vals[7]) || std::isnan(vals[8]);
    }

    void operator=(const Mat3 &m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                (*this)(i, j) = m(i, j);
    }

    void operator+=(const Mat3 &m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                (*this)(i, j) += m(i, j);
    }

    void operator-=(const Mat3 &m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                (*this)(i, j) -= m(i, j);
    }

    void operator/=(double s) {
        for (unsigned int c = 0; c < 9; ++c)
            vals[c] /= s;
    }

    Mat3 operator-(const Mat3 &m2) {
        return Mat3((*this)(0, 0) - m2(0, 0), (*this)(0, 1) - m2(0, 1), (*this)(0, 2) - m2(0, 2),
                    (*this)(1, 0) - m2(1, 0), (*this)(1, 1) - m2(1, 1), (*this)(1, 2) - m2(1, 2),
                    (*this)(2, 0) - m2(2, 0), (*this)(2, 1) - m2(2, 1), (*this)(2, 2) - m2(2, 2));
    }

    Mat3 operator+(const Mat3 &m2) {
        return Mat3((*this)(0, 0) + m2(0, 0), (*this)(0, 1) + m2(0, 1), (*this)(0, 2) + m2(0, 2),
                    (*this)(1, 0) + m2(1, 0), (*this)(1, 1) + m2(1, 1), (*this)(1, 2) + m2(1, 2),
                    (*this)(2, 0) + m2(2, 0), (*this)(2, 1) + m2(2, 1), (*this)(2, 2) + m2(2, 2));
    }

    Mat3 operator/(float s) {
        return Mat3((*this)(0, 0) / s, (*this)(0, 1) / s, (*this)(0, 2) / s, (*this)(1, 0) / s, (*this)(1, 1) / s,
                    (*this)(1, 2) / s, (*this)(2, 0) / s, (*this)(2, 1) / s, (*this)(2, 2) / s);
    }

    Mat3 operator*(float s) {
        return Mat3((*this)(0, 0) * s, (*this)(0, 1) * s, (*this)(0, 2) * s, (*this)(1, 0) * s, (*this)(1, 1) * s,
                    (*this)(1, 2) * s, (*this)(2, 0) * s, (*this)(2, 1) * s, (*this)(2, 2) * s);
    }

    ////////        ACCESS TO COORDINATES      /////////
    float operator()(unsigned int i, unsigned int j) const { return vals[3 * i + j]; }

    float &operator()(unsigned int i, unsigned int j) { return vals[3 * i + j]; }

    ////////        BASICS       /////////
    inline float sqrnorm() {
        return vals[0] * vals[0] + vals[1] * vals[1] + vals[2] * vals[2]
                + vals[3] * vals[3] + vals[4] * vals[4] + vals[5] * vals[5]
                + vals[6] * vals[6] + vals[7] * vals[7] + vals[8] * vals[8];
    }

    inline float norm() { return sqrt(sqrnorm()); }

    inline float determinant() const {
        return vals[0] * (vals[4] * vals[8] - vals[7] * vals[5])
                - vals[1] * (vals[3] * vals[8] - vals[6] * vals[5])
                + vals[2] * (vals[3] * vals[7] - vals[6] * vals[4]);
    }

    inline float trace() const { return vals[0] + vals[4] + vals[8]; }

    ////////        TRANSPOSE       /////////
    inline
    void transpose() {
        float xy = vals[1], xz = vals[2], yz = vals[5];
        vals[1] = vals[3];
        vals[3] = xy;
        vals[2] = vals[6];
        vals[6] = xz;
        vals[5] = vals[7];
        vals[7] = yz;
    }

    Mat3 getTranspose() const {
        return Mat3(vals[0], vals[3], vals[6], vals[1], vals[4], vals[7], vals[2], vals[5], vals[8]);
    }

    Mat3 operator-() const {
        return Mat3(-vals[0], -vals[1], -vals[2], -vals[3], -vals[4], -vals[5], -vals[6], -vals[7], -vals[8]);
    }


private:
    float vals[9];
    // will be noted as :
    // 0 1 2
    // 3 4 5
    // 6 7 8
};


inline static
Mat3 operator*(float s, const Mat3 &m) {
    return Mat3(m(0, 0) * s, m(0, 1) * s, m(0, 2) * s, m(1, 0) * s, m(1, 1) * s, m(1, 2) * s, m(2, 0) * s, m(2, 1) * s,
                m(2, 2) * s);
}


inline static std::ostream &operator<<(std::ostream &s, Mat3 const &m) {
    s << m(0, 0) << " \t" << m(0, 1) << " \t" << m(0, 2) << std::endl << m(1, 0) << " \t" << m(1, 1) << " \t" << m(1, 2)
      << std::endl << m(2, 0) << " \t" << m(2, 1) << " \t" << m(2, 2) << std::endl;
    return s;
}


#endif
