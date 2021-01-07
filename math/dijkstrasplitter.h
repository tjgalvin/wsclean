#ifndef SUBDIVISION_H
#define SUBDIVISION_H

#include <aocommon/uvector.h>

#include <algorithm>
#include <cmath>
#include <queue>

class DijkstraSplitter {
 public:
  DijkstraSplitter(size_t width, size_t height)
      : _width(width), _height(height) {}

  struct Coord {
    Coord() = default;
    Coord(size_t _x, size_t _y) : x(_x), y(_y) {}
    size_t x, y;
  };

  struct Visit {
    float distance;
    Coord to, from;
    bool operator<(const Visit& rhs) const { return distance > rhs.distance; }
  };

  void AddVerticalDivider(const float* image, float* scratch, float* output,
                          size_t x1, size_t x2) const {
    DivideVertically(image, scratch, x1, x2);
    for (size_t y = 0; y != _height; ++y) {
      for (size_t i = y * _width + x1; i != y * _width + x2; ++i)
        output[i] += scratch[i];
    }
  }

  /**
   * Find the shortest vertical path through an image. The path is
   * constrained to lie between vertical lines given by x1 and x2.
   * The output is set to 1 for pixels that are part of the path, and
   * set to 0 otherwise. The reason it's a floating point is because
   * it is also used as scratch.
   */
  void DivideVertically(const float* image, float* output, size_t x1,
                        size_t x2) const {
    using visitset = std::priority_queue<Visit>;
    visitset visits;

    for (size_t x = x1; x != x2; ++x) {
      Visit visit;
      visit.distance = 0.0;
      visit.to = Coord(x, 0);
      visit.from = Coord(x, 0);
      visits.push(visit);
    }
    aocommon::UVector<Coord> path((x2 - x1) * _height);
    fillColumns(output, x1, x2, std::numeric_limits<float>::max());
    Visit visit;
    while (!visits.empty()) {
      visit = visits.top();
      visits.pop();
      size_t x = visit.to.x, y = visit.to.y;
      if (y == _height) break;
      float curDistance = output[x + y * _width];
      float newDistance = visit.distance + std::fabs(image[x + y * _width]);
      if (newDistance < curDistance) {
        output[x + y * _width] = newDistance;
        path[x - x1 + y * (x2 - x1)] = visit.from;
        visit.distance = newDistance;
        visit.from = visit.to;
        if (x > x1) {
          visit.to = Coord(x - 1, y + 1);
          visits.push(visit);
          visit.to = Coord(x - 1, y);
          visits.push(visit);
        }
        visit.to = Coord(x, y + 1);
        visits.push(visit);
        if (x < x2 - 1) {
          visit.to = Coord(x + 1, y + 1);
          visits.push(visit);
          visit.to = Coord(x + 1, y);
          visits.push(visit);
        }
      }
    }
    fillColumns(output, x1, x2, 0.0);
    Coord pCoord = visit.from;
    while (pCoord.y > 0) {
      output[pCoord.x + pCoord.y * _width] = 1.0;
      pCoord = path[pCoord.x - x1 + pCoord.y * (x2 - x1)];
    }
    output[pCoord.x] = 1.0;
  }

  void AddHorizontalDivider(const float* image, float* scratch, float* output,
                            size_t y1, size_t y2) const {
    DivideHorizontally(image, scratch, y1, y2);
    for (size_t i = y1 * _width; i != y2 * _width; ++i) output[i] += scratch[i];
  }

  /**
   * Like DivideVertically, but for horizontal paths. The path is
   * constrained to lie between horizontal lines given by y1 and y2.
   */
  void DivideHorizontally(const float* image, float* output, size_t y1,
                          size_t y2) const {
    using visitset = std::priority_queue<Visit>;
    visitset visits;

    for (size_t y = y1; y != y2; ++y) {
      Visit visit;
      visit.distance = 0.0;
      visit.to = Coord(0, y);
      visit.from = Coord(0, y);
      visits.push(visit);
    }
    aocommon::UVector<Coord> path(_width * (y2 - y1));
    std::fill(output + y1 * _width, output + y2 * _width,
              std::numeric_limits<float>::max());
    Visit visit;
    while (!visits.empty()) {
      visit = visits.top();
      visits.pop();
      size_t x = visit.to.x, y = visit.to.y;
      if (x == _width) break;
      float curDistance = output[x + y * _width];
      float newDistance = visit.distance + std::fabs(image[x + y * _width]);
      if (newDistance < curDistance) {
        output[x + y * _width] = newDistance;
        path[x + (y - y1) * _width] = visit.from;
        visit.distance = newDistance;
        visit.from = visit.to;
        if (y > y1) {
          visit.to = Coord(x + 1, y - 1);
          visits.push(visit);
          visit.to = Coord(x, y - 1);
          visits.push(visit);
        }
        visit.to = Coord(x + 1, y);
        visits.push(visit);
        if (y < y2 - 1) {
          visit.to = Coord(x + 1, y + 1);
          visits.push(visit);
          visit.to = Coord(x, y + 1);
          visits.push(visit);
        }
      }
    }
    std::fill(output + y1 * _width, output + y2 * _width, 0.0);
    Coord pCoord = visit.from;
    while (pCoord.x > 0) {
      output[pCoord.x + pCoord.y * _width] = 1.0;
      pCoord = path[pCoord.x + (pCoord.y - y1) * _width];
    }
    output[pCoord.y * _width] = 1.0;
  }

  /**
   * Mask the space between (typically) two vertical divisions.
   * @param subdivision An image that is the result of earlier calls
   * to DivideVertically().
   * @param subImgX An x-position that is in between the two splits.
   * @param mask A mask image for which pixels will be set to true if
   *   and only if they are part of the area specified by the
   *   two divisions.
   * @param [out] x The left side of the bounding box of the divisions.
   * @param [out] subWidth The bounding width of the two divisions.
   */
  void FloodVerticalArea(const float* subdivision, size_t subImgX, bool* mask,
                         size_t& x, size_t& subWidth) const {
    std::fill(mask, mask + _width * _height, false);
    x = _width;
    size_t x2 = 0;
    for (size_t y = 0; y != _height; ++y) {
      const float* dRow = &subdivision[y * _width];
      bool* maskRow = &mask[y * _width];
      int xIter = subImgX;
      // Move to the left until a border is hit
      while (xIter >= 0 && dRow[xIter] == 0.0) {
        maskRow[xIter] = true;
        --xIter;
      }
      // Continue to the left through the border
      while (xIter >= 0 && dRow[xIter] != 0.0) {
        maskRow[xIter] = true;
        --xIter;
      }
      x = std::min<size_t>(x, xIter + 1);
      xIter = subImgX + 1;
      // Move to the right until a border is hit
      while (size_t(xIter) < _width && dRow[xIter] == 0.0) {
        maskRow[xIter] = true;
        ++xIter;
      }
      x2 = std::max<size_t>(x2, xIter);
    }
    if (x2 < x)
      subWidth = 0;
    else
      subWidth = x2 - x;
  }

  /**
   * Like FloodVerticalArea(), but for horizontal flooding.
   */
  void FloodHorizontalArea(const float* subdivision, size_t subImgY, bool* mask,
                           size_t& y, size_t& subHeight) const {
    std::fill(mask, mask + _width * _height, false);
    y = _height;
    size_t y2 = 0;
    for (size_t x = 0; x != _width; ++x) {
      int yIter = subImgY;
      // Move up until a border is hit
      while (yIter >= 0 && subdivision[yIter * _width + x] == 0.0) {
        mask[yIter * _width + x] = true;
        --yIter;
      }
      // Continue to the left through the border
      while (yIter >= 0 && subdivision[yIter * _width + x] != 0.0) {
        mask[yIter * _width + x] = true;
        --yIter;
      }
      y = std::min<size_t>(y, yIter + 1);
      yIter = subImgY + 1;
      // Move to the right until a border is hit
      while (size_t(yIter) < _height &&
             subdivision[yIter * _width + x] == 0.0) {
        mask[yIter * _width + x] = true;
        ++yIter;
      }
      y2 = std::max<size_t>(y2, yIter);
    }
    if (y2 < y)
      subHeight = 0;
    else
      subHeight = y2 - y;
  }

  /**
   * Combines a horizontally and vertically filled area and extracts a
   * single mask where the areas overlap.
   * @param vMask Mask returned by FloodHorizontalArea(), but trimmed
   * to have the specified width.
   * @param vMaskX x-value returned by FloodHorizontalArea().
   * @param vMaskWidth Width return by FloodHorizontalArea(), and width
   * of vMask.
   * @param hMask Mask returned by FloodVerticalArea().
   * @param [in,out] mask Result
   * @param [out] subX Bounding box x-value
   * @param [out] subY Bounding box y-value
   * @param [out] subWidth Bounding box width
   * @param [out] subHeight Bounding box height
   */
  void GetBoundingMask(const bool* vMask, size_t vMaskX, size_t vMaskWidth,
                       const bool* hMask, bool* mask, size_t& subX,
                       size_t& subY, size_t& subWidth,
                       size_t& subHeight) const {
    subX = vMaskWidth + vMaskX;
    subY = _height;
    size_t subX2 = 0, subY2 = 0;
    for (size_t y = 0; y != _height; ++y) {
      const bool* vMaskRow = &vMask[y * vMaskWidth];
      const bool* hMaskRow = &hMask[y * _width];
      bool* maskRow = &mask[y * _width];
      for (size_t x = 0; x != vMaskWidth; ++x) {
        size_t hx = x + vMaskX;
        if (vMaskRow[x] && hMaskRow[hx]) {
          maskRow[hx] = true;
          subX = std::min(hx, subX);
          subY = std::min(y, subY);
          subX2 = std::max(hx, subX2);
          subY2 = std::max(y, subY2);
        } else
          maskRow[hx] = false;
      }
    }
    if (subX2 < subX) {
      subWidth = 0;
      subHeight = 0;
    } else {
      subWidth = subX2 + 1 - subX;
      subHeight = subY2 + 1 - subY;
    }
    // If dimensions start off even, keep subimages even too
    if (_width % 2 == 0) {
      if (subWidth % 2 != 0) {
        ++subWidth;
        if (subWidth + subX >= _width) --subX;
      }
    }
    if (_height % 2 == 0) {
      if (subHeight % 2 != 0) {
        ++subHeight;
        if (subHeight + subY >= _height) --subY;
      }
    }
  }

 private:
  /**
   * This function sets a rectangular area given by 0 <= y < height and xStart
   * <= x < xEnd.
   */
  void fillColumns(float* output, size_t xStart, size_t xEnd,
                   float newValue) const {
    for (size_t y = 0; y != _height; ++y) {
      std::fill(output + _width * y + xStart, output + _width * y + xEnd,
                newValue);
    }
  }
  size_t _width, _height;
};

#endif
