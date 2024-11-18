#ifndef GRAPHICS_HPP
#define GRAPHICS_HPP

#include "kinematics.hpp"
#include <SFML/Graphics.hpp>
#include <cassert>
#include <cmath>

class Plot
{
 private:
  const float w_;
  const float h_;

  const float x_scale_;
  const float y_scale_;

  sf::VertexArray x_axis{sf::Lines, 2};
  sf::VertexArray y_axis{sf::Lines, 2};

  sf::VertexArray barrier1_line;
  sf::VertexArray barrier2_line;

  std::vector<sf::Vertex> vertices_px;

  template<typename FP1, typename FP2>
  sf::Vector2f to_pixel(FP1 x, FP2 y);

  sf::Vector2f to_pixel(Point p);

 public:
  // creates barrier lines, plot scale is set from barrier dimensions
  Plot(const Barrier& barrier1, const Barrier& barrier2, int w = 600,
       int h = 400);

  void add(const std::vector<Point>& vertices);

  void show();
};
#endif