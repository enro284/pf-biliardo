#include "graphics.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>

template<typename FP1, typename FP2>
sf::Vector2f Plot::to_pixel(FP1 x, FP2 y)
{
  return (sf::Vector2f(static_cast<float>(x) * scale_ + 0.1f * w_,
                       static_cast<float>(-y) * scale_ + 0.5f * h_));
}

sf::Vector2f Plot::to_pixel(Vec2 p)
{
  return to_pixel(p.x_, p.y_);
}

Plot::Plot(const Barrier& barrier1, const Barrier& barrier2, int w, int h)
    : w_{static_cast<float>(w)}
    , h_{static_cast<float>(h)}
    , scale_{std::min(0.8f * w_ / static_cast<float>(barrier1.max()),
                      0.4f * h_ / static_cast<float>(barrier1.pol()(0.)))}
    , barrier1_plot(sf::LineStrip, 2)
    , barrier2_plot(sf::LineStrip, 2)
    , vertices_px(0)
{
  assert(w_ > 0 && h_ > 0);

  x_axis[0] = sf::Vertex({0, h_ / 2.f}, sf::Color::Black);
  x_axis[1] = sf::Vertex({w_, h_ / 2.f}, sf::Color::Black);

  y_axis[0] = sf::Vertex({w_ * 0.1f, 0}, sf::Color::Black);
  y_axis[1] = sf::Vertex({w_ * 0.1f, h_}, sf::Color::Black);

  if (barrier1.pol().deg() == 1) {
    barrier1_plot[0] = {to_pixel(0., barrier1.pol()(0)), sf::Color::Black};
    barrier1_plot[1] = {
        to_pixel(barrier1.max(), barrier1.pol()(barrier1.max())),
        sf::Color::Black};

    barrier2_plot[0] = {to_pixel(0., barrier2.pol()(0)), sf::Color::Black};
    barrier2_plot[1] = {
        to_pixel(barrier2.max(), barrier2.pol()(barrier2.max())),
        sf::Color::Black};
  } else {
    std::size_t div = static_cast<std::size_t>(w_ * 0.8f);

    double x{0.};
    double dx{barrier1.max() / static_cast<double>(div)};

    barrier1_plot.resize(div + 1);
    barrier2_plot.resize(div + 1);

    for (std::size_t i{0}; i <= div; ++i) {
      barrier1_plot[i] = {to_pixel(x, barrier1.pol()(x)), sf::Color::Black};
      barrier2_plot[i] = {to_pixel(x, barrier2.pol()(x)), sf::Color::Black};
      x += dx;
    }
  }
}

void Plot::add(const std::vector<Vec2>& vertices)
{
  vertices_px.resize(vertices.size());
  std::transform(
      vertices.begin(), vertices.end(), vertices_px.begin(),
      [&](const Vec2& p) { return sf::Vertex(to_pixel(p), sf::Color::Red); });
}

void Plot::show()
{
  sf::RenderWindow window(sf::VideoMode(static_cast<unsigned int>(w_),
                                        static_cast<unsigned int>(h_)),
                          "Biliardo triangolare");
  window.setFramerateLimit(60);
  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) {
        window.close();
      } else if (event.type == sf::Event::KeyPressed
                 && event.key.code == sf::Keyboard::Q) {
        window.close();
      }
    }

    window.clear(sf::Color::White);

    window.draw(x_axis);
    window.draw(y_axis);
    window.draw(barrier1_plot);
    window.draw(barrier2_plot);

    window.draw(&vertices_px[0], vertices_px.size(), sf::LineStrip);

    window.display();
  }
}
