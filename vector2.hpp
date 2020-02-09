#pragma once

#include <iostream>

#include "config.hpp"

template <typename T>
class Vector2 {
  public:
    T x, y;
    SPH_KERNEL
    Vector2()                       : x(T(0)) , y(T(0))  {}
    SPH_KERNEL
    Vector2(const T _x, const T _y) : x(_x)   , y(_y)    {}
    SPH_KERNEL
    Vector2(const T s)              : x(s)    , y(s)     {}
    SPH_KERNEL
    Vector2(const Vector2& src)     : x(src.x), y(src.y) {}

    SPH_KERNEL
    const Vector2& operator = (const Vector2& rhs) {
      x = rhs.x;
      y = rhs.y;
      return (*this);
    }

    SPH_KERNEL
    const Vector2& operator = (const T s) {
      x = y = s;
      return (*this);
    }

    // addition
    SPH_KERNEL
    Vector2 operator + (const Vector2& rhs) const {
      return Vector2(x + rhs.x, y + rhs.y);
    }

    SPH_KERNEL
    const Vector2& operator += (const Vector2& rhs) {
      (*this) = (*this) + rhs;
      return (*this);
    }

    // subtraction
    SPH_KERNEL
    Vector2 operator - (const Vector2& rhs) const {
      return Vector2(x - rhs.x, y - rhs.y);
    }

    SPH_KERNEL
    const Vector2& operator -= (const Vector2& rhs) {
      (*this) = (*this) - rhs;
      return (*this);
    }

    // product
    SPH_KERNEL
    Vector2 operator * (const T s) const {
      return Vector2(x * s, y * s);
    }

    SPH_KERNEL
    const Vector2& operator *= (const T s) {
      (*this) = (*this) * s;
      return (*this);
    }

    SPH_KERNEL
    friend Vector2 operator * (const T s, const Vector2& v) {
      return (v * s);
    }

    // dividion
    SPH_KERNEL
    Vector2 operator / (const T s) const {
      return Vector2(x / s, y / s);
    }

    SPH_KERNEL
    const Vector2& operator /= (const T s) {
      (*this) = (*this) / s;
      return (*this);
    }

    // sign
    SPH_KERNEL
    const Vector2& operator + () const {
      return (*this);
    }

    SPH_KERNEL
    const Vector2 operator - () const {
      return Vector2(-x, -y);
    }

    // inner product
    SPH_KERNEL
    T operator * (const Vector2& rhs) const {
      return (x * rhs.x) + (y * rhs.y);
    }

    // outer product (retruned value is scholar)
    SPH_KERNEL
    T operator ^ (const Vector2& rhs) const {
      const T z = (x * rhs.y) - (y * rhs.x);
      return z;
    }

    // cast to Vector2<U>
    template <typename U>
    SPH_KERNEL
    operator Vector2<U> () const {
      return Vector2<U>(static_cast<U>(x), static_cast<U>(y));
    }

    // min/max
    SPH_KERNEL
    T getMin() const {
      return x < y ? x : y;
    }

    SPH_KERNEL
    T getMax() const {
      return x > y ? x : y;
    }

    // apply
    template <class F>
    SPH_KERNEL
    Vector2 applyEach(F f) const {
      return Vector2(f(x), f(y));
    }

    template <class F>
    SPH_KERNEL
    friend Vector2 ApplyEach(F f, const Vector2& arg1, const Vector2& arg2) {
      return Vector2(f(arg1.x, arg2.x), f(arg1.y, arg2.y));
    }

    // stream
    friend std::ostream& operator << (std::ostream& c, const Vector2& u) {
      c << u.x << " " << u.y;
      return c;
    }

    friend std::istream& operator >> (std::istream& c, Vector2& u) {
      c >> u.x; c >> u.y;
      return c;
    }

    // index
    const T& operator [] (const int i) const {
      if (0 == i) return x;
      if (1 == i) return y;
      std::cerr << "Vector index = " << i << " is not valid." << std::endl;
      exit(-1);
      return x; // dummy for avoid warning
    }

    T& operator [] (const int i) {
      if (0 == i) return x;
      if (1 == i) return y;
      std::cerr << "Vector index = " << i << " is not valid." << std::endl;
      exit(-1);
      return x; // dummy for avoid warning
    }

    SPH_KERNEL
    T getDistanceSQ(const Vector2& u) const {
      T dx = x - u.x;
      T dy = y - u.y;
      return dx * dx + dy * dy;
    }

    // comparison
    SPH_KERNEL
    bool operator == (const Vector2& u) const {
      return (x == u.x) && (y == u.y);
    }

    SPH_KERNEL
    bool operator != (const Vector2& u) const {
      return (x != u.x) || (y != u.y);
    }

    // orthant
    SPH_KERNEL
    int orthant(const Vector2& origin) const {
	    return (x > origin.x) + ((y > origin.y) << 1);
    }
};

// optimization for division
template <>
SPH_KERNEL
inline Vector2<float> Vector2<float>::operator / (const float s) const {
  const float inv_s = 1.0f / s;
  return Vector2(x * inv_s, y * inv_s);
}

template <>
SPH_KERNEL
inline Vector2<double> Vector2<double>::operator / (const double s) const {
  const double inv_s = 1.0 / s;
  return Vector2(x * inv_s, y * inv_s);
}
