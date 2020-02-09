#pragma once

#include <limits>

#include "config.hpp"
#include "util.hpp"
#include "vector2.hpp"

template <typename T>
class BoundingBox2 {
  public:
    Vector2<T> min, max;
    BoundingBox2() : min(std::numeric_limits<T>::max()), max(std::numeric_limits<T>::min()) {}
    BoundingBox2(const BoundingBox2& src) : min(src.min), max(src.max) {}
    BoundingBox2(const Vector2<T>& vec) : min(vec), max(vec) {}
    BoundingBox2(const Vector2<T>& vec1, const Vector2<T>& vec2) {
      min.x = min_(vec1.x, vec2.x);
      min.y = min_(vec1.y, vec2.y);
      max.x = max_(vec1.x, vec2.x);
      max.y = max_(vec1.y, vec2.y);
    }

    const BoundingBox2& operator = (const BoundingBox2& rhs) {
      min = rhs.min;
      max = rhs.max;
      return (*this);
    }

    const BoundingBox2& merge(const Vector2<T>& vec) {
      max.x = max_(max.x, vec.x);
      max.y = max_(max.y, vec.y);
      min.x = min_(min.x, vec.x);
      min.y = min_(min.y, vec.y);
      return (*this);
    }

    const BoundingBox2& merge(const BoundingBox2& bbox) {
      max.x = max_(max.x, bbox.max.x);
      max.y = max_(max.y, bbox.max.y);
      min.x = min_(min.x, bbox.min.x);
      min.y = min_(min.y, bbox.min.y);
      return (*this);
    }

    const BoundingBox2& expand(const T margin) {
      max.x += margin;
      max.y += margin;
      min.x -= margin;
      min.y -= margin;
      return (*this);
    }

    inline Vector2<T> center() const {
      return (max + min) * T(0.5);
    }

    inline bool intersect(const BoundingBox2& bbox) const {
      return (bbox.min.x <= max.x) && (min.x <= bbox.max.x)
          && (bbox.min.y <= max.y) && (min.y <= bbox.max.y);
    }

    inline BoundingBox2 orthant(const int i) const {
      int px = i & 1;
      int py = (i & (1 << 1)) >> 1;
      T x = px ? max.x : min.x;
      T y = py ? max.y : min.y;
      return BoundingBox2(Vector2<T>(x, y), center());
    }

    // stream
    friend std::ostream& operator << (std::ostream& c, const BoundingBox2& bbox) {
      c << bbox.min << " " << bbox.max;
      return c;
    }
};
