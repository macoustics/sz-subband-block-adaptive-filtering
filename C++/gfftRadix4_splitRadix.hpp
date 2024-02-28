#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

namespace fft{
using uint = unsigned int;
static const double M_PI {3.14159265358979323846};

template <unsigned N, typename C, bool F = false>
inline C Twiddle(const unsigned int& n)
{
    return F ? std::exp(C(0, n * 2.0 * M_PI / ((double)N))) : std::exp(C(0, -n * 2.0 * M_PI / ((double)N)));
};

static constexpr double invSqrt2 { 1.0 / std::sqrt(2) };

template <typename C, bool F>
inline C rot90(C input)
{
    return F ? C(-input.imag(), input.real()) : C(input.imag(), -input.real());
}

template <typename C, bool F>
inline C rot45(C input)
{
    return F ? (C(input.real() - input.imag(), input.real() + input.imag())) * invSqrt2 : (C(input.real() + input.imag(), -input.real() + input.imag())) * invSqrt2;
}

template <typename C, bool F>
inline C rot135(C input)
{
    return F ? (C(-input.real() - input.imag(), input.real() - input.imag())) * invSqrt2 : (C(-input.real() + input.imag(), -input.real() - input.imag())) * invSqrt2;
}

template <unsigned N, typename A, typename C, typename T, size_t M, bool F = false>
class DanielsonLanczos {
    DanielsonLanczos<N / 2, A, C, T, M, F> next2;
    DanielsonLanczos<N / 4, A, C, T, M, F> next4;

private:
    std::array<C, N / 2> m_ScratchEven {};
    std::array<C, N / 4> m_ScratchOdd1 {};
    std::array<C, N / 4> m_ScratchOdd3 {};

public:
    void apply(A& data, const A& twiddle, const size_t& offset, const uint& stride)
    {
        next2.apply(data, twiddle, offset, 2 * stride);
        next4.apply(data, twiddle, offset + stride, 4 * stride);
        next4.apply(data, twiddle, offset - stride, 4 * stride);

        for (unsigned i = 0; i < N / 4; i++) {
            m_ScratchEven[i] = data[(offset + i * 2 * stride) & M];
            m_ScratchEven[i + N / 4] = data[(offset + i * 2 * stride + N / 2 * stride) & M];
            m_ScratchOdd1[i] = data[(offset + i * 4 * stride + stride) & M] * twiddle[i * stride];
            m_ScratchOdd3[i] = data[(offset + i * 4 * stride - stride) & M] * std::conj(twiddle[i * stride]);
        }
        
        for (unsigned i = 0; i < N / 4; i++) {
            // Calculate: data[i*S*2] = u1 + z1 + z3
            data[(offset + i * stride) & M] = m_ScratchEven[i] + m_ScratchOdd1[i] + m_ScratchOdd3[i];

            // Calculate: data[i*S*2 + N/2] = u1 - z1 - z3
            data[(offset + i * stride + N / 2 * stride) & M] = m_ScratchEven[i] - (m_ScratchOdd1[i] + m_ScratchOdd3[i]);

            // Calculate: data[i*S*2 + N/4] = u3 -j*(z1 - z3)
            data[(offset + i * stride + N / 4 * stride) & M] = m_ScratchEven[i + N / 4] + rot90<C, F>(m_ScratchOdd1[i] - m_ScratchOdd3[i]);

            // Calculate: data[i*S*2 + 3*N/4] = u3 + j*(z1 - z3)
            data[(offset + i * stride + 3 * N / 4 * stride) & M] = m_ScratchEven[i + N / 4] - rot90<C, F>(m_ScratchOdd1[i] - m_ScratchOdd3[i]);
        }
    }
};

template <typename A, typename C, typename T, size_t M, bool F>
class DanielsonLanczos<1, A, C, T, M, F> {
public:
    void apply(A& data, const A& twiddle, const size_t& offset, const uint& stride)
    {
        C y { data[0] };
        y = twiddle[0];
        size_t t {};
        t = offset;
        t++;
        uint u {};
        u = stride;
        u++;
    }
};

template <typename A, typename C, typename T, size_t M, bool F>
class DanielsonLanczos<2, A, C, T, M, F> {
public:
    void apply(A& data, const A& twiddle, const size_t& offset, const uint& stride)
    {
        C y0, y1;
        y0 = data[offset & M];
        y1 = twiddle[0];
        y1 = data[(offset + stride) & M];
        // Calculate: data[i] = y0 + y1;
        data[(offset)&M] = y0 + y1;
        // Calculate: data[i+N/2] = y0 - y1;
        data[(offset + stride) & M] = y0 - y1;
    }
};

template <typename A, typename C, typename T, size_t M, bool F>
class DanielsonLanczos<4, A, C, T, M, F> {
public:
    void apply(A& data, const A& twiddle, const size_t& offset, const uint& stride)
    {
        C y0, y1, y2, y3;
        y1 = twiddle[0];
        y0 = data[(offset)&M];
        y1 = data[(offset + stride) & M];
        y2 = data[(offset + 2 * stride) & M];
        y3 = data[(offset + 3 * stride) & M];
        // Calculate: data[i] += y1 + y2 + y3;
        data[offset & M] = y0 + y1 + y2 + y3;
        // Calculate: data[i+N/4] = data[i] - j*y1 -y2 + j*y3;
        data[(offset + stride) & M] = y0 + rot90<C, F>(y1) - y2 - rot90<C, F>(y3);
        // Calculate: data[i+N/2] = data[i] - y1 + y2 - y3;
        data[(offset + 2 * stride) & M] = y0 - y1 + y2 - y3;
        // Calculate: data[i+3*N/4] = data[i] + j*y1 - y2 - j*y3;
        data[(offset + 3 * stride) & M] = y0 - rot90<C, F>(y1) - y2 + rot90<C, F>(y3);
    }
};

template <typename A, typename C, typename T, size_t M, bool F>
class DanielsonLanczos<8, A, C, T, M, F> {
public:
    void apply(A& data, const A& twiddle, const size_t& offset, const uint& stride)
    {
        C y0, y1, y2, y3, y4, y5, y6, y7;
        y1 = twiddle[0];
        y0 = data[offset & M];
        y1 = data[(offset + stride) & M];
        y2 = data[(offset + 2* stride) & M];
        y3 = data[(offset + 3* stride) & M];
        y4 = data[(offset + 4* stride) & M];
        y5 = data[(offset + 5* stride) & M];
        y6 = data[(offset + 6* stride) & M];
        y7 = data[(offset + 7* stride) & M];

        // Calculate: data[i] += y1 + y2 + y3;
        data[offset & M] = y0 + y1 + y2 + y3 + y4 + y5 + y6 + y7;

        // data[i+1] = y0 + w*y1 + w^2*y2 + w^3*y3 + w^4*y4 + w^5*y5 + w^6*y6 + w^7*y7;
        data[(offset + stride) & M] = y0 + rot45<C, F>(y1) + rot90<C, F>(y2) + rot135<C, F>(y3) - y4 - rot45<C, F>(y5) - rot90<C, F>(y6) - rot135<C, F>(y7);

        // data[i+2] = y0 + w^2*y1 + w^4*y2 + w^6*y3 + y4 + w^2*y5 + w^4*y6 + w^6*y7;
        data[(offset + 2 * stride) & M] = y0 + rot90<C, F>(y1) - y2 - rot90<C, F>(y3) + y4 + rot90<C, F>(y5) - y6 - rot90<C, F>(y7);

        // data[i+3] = y0 + w^3*y1 + w^6*y2 + w^9*y3 + w^12*y4 + w^15*y5 + w^18*y6 + w^21*y7;
        // data[i+3] = y0 + w^3*y1 + w^6*y2 + w^1*y3 + w^4*y4 + w^7*y5 + w^2*y6 + w^5*y7;
        data[(offset + 3 * stride) & M] = y0 + rot135<C, F>(y1) - rot90<C, F>(y2) + rot45<C, F>(y3) - y4 - rot135<C, F>(y5) + rot90<C, F>(y6) - rot45<C, F>(y7);

        // data[i+4] = y0 + w^4*y1 + w^8*y2 + w^12*y3 + w^16*y4 + w^20*y5 + w^24*y6 + w^28*y7;
        // data[i+4] = y0 + w^4*y1 + y2 + w^4*y3 + y4 + w^4*y5 + y6 + w^4*y7;
        data[(offset + 4 * stride) & M] = y0 - y1 + y2 - y3 + y4 - y5 + y6 - y7;

        // data[i+5] = y0 + w^5*y1 + w^10*y2 + w^15*y3 + w^20*y4 + w^25*y5 + w^30*y6 + w^35*y7;
        // data[i+5] = y0 + w^5*y1 + w^2*y2 + w^7*y3 + w^4*y4 + w^1*y5 + w^6*y6 + w^3*y7;
        data[(offset + 5 * stride) & M] = y0 - rot45<C, F>(y1) + rot90<C, F>(y2) - rot135<C, F>(y3) - y4 + rot45<C, F>(y5) - rot90<C, F>(y6) + rot135<C, F>(y7);

        // data[i+6] = y0 + w^6*y1 + w^12*y2 + w^18*y3 + w^24*y4 + w^30*y5 + w^36*y6 + w^42*y7;
        // data[i+6] = y0 + w^6*y1 + w^4*y2 + w^2*y3 + y4 + w^6*y5 + w^4*y6 + w^2*y7;
        data[(offset + 6 * stride) & M] = y0 - rot90<C, F>(y1) - y2 + rot90<C, F>(y3) + y4 - rot90<C, F>(y5) - y6 + rot90<C, F>(y7);

        // data[i+7] = y0 + w^7*y1 + w^14*y2 + w^21*y3 + w^28*y4 + w^35*y5 + w^42*y6 + w^49*y7;
        // data[i+7] = y0 + w^7*y1 + w^6*y2 + w^5*y3 + w^4*y4 + w^3*y5 + w^2*y6 + w^1*y7;
        data[(offset + 7 * stride) & M] = y0 - rot135<C, F>(y1) - rot90<C, F>(y2) - rot45<C, F>(y3) - y4 + rot135<C, F>(y5) + rot90<C, F>(y6) + rot45<C, F>(y7);
    }
};


// F is a flag indicating inverse transform
// F = true performs the inverse transform (without normalization)
template <unsigned N, typename A, typename T, bool F>
class gfft {
    using ArrayType = A;
    typedef std::complex<T> C;
    // const static size_t Mask {N};
    DanielsonLanczos<N, A, C, T, N - 1, F> recursion;

public:
    gfft()
        : m_twiddle {}
    {
        populateTwiddle(m_twiddle);
    }

    ~gfft() = default;
    void transform(A& data)
    {
        recursion.apply(data, m_twiddle, 0, 1);
    }

private:
    ArrayType m_twiddle;

    void populateTwiddle(ArrayType& twiddleFactors)
    {
        T pi {std::acos((T)-1)};
        for (size_t i = 0; i < N; i++) {
            twiddleFactors[i] = F ? std::exp(C(0, (T)i * 2.0 * pi / ((double)N))) : std::exp(C(0, -(T)i * 2.0 * pi / ((double)N)));
        }
    }
};
} // end namespace