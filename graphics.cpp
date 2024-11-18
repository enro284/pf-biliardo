#include "graphics.hpp"
#include "algorithm"

template<typename FP1, typename FP2>
sf::Vector2f Plot::to_pixel(FP1 x, FP2 y)
{
  return (sf::Vector2f(static_cast<float>(x) * x_scale_ + 0.1f * w_,
                       static_cast<float>(-y) * y_scale_ + 0.5f * h_));
}

sf::Vector2f Plot::to_pixel(Point p)
{
  return to_pixel(p.x_, p.y_);
}

Plot::Plot(const Barrier& barrier1, const Barrier& barrier2, int w, int h)
    : w_{static_cast<float>(w)}
    , h_{static_cast<float>(h)}
    , x_scale_{w_ * 0.8f / static_cast<float>(barrier1.max())}
    , y_scale_{h_ * 0.8f / static_cast<float>(barrier1.pol()(0.) * 2.f)}
    , barrier1_line(sf::Lines, 2)
    , barrier2_line(sf::Lines, 2)
    , vertices_px(0)
{
  assert(w_ > 0 && h_ > 0);
  x_axis[0] = sf::Vertex({0, h_ / 2.f}, sf::Color::Black);
  x_axis[1] = sf::Vertex({w_, h_ / 2.f}, sf::Color::Black);

  y_axis[0] = sf::Vertex({w_ * 0.1f, 0}, sf::Color::Black);
  y_axis[1] = sf::Vertex({w_ * 0.1f, h_}, sf::Color::Black);

  barrier1_line[0] =
      sf::Vertex(to_pixel({0., barrier1.pol()(0)}), sf::Color::Black);
  barrier1_line[1] =
      sf::Vertex(to_pixel({barrier1.max(), barrier1.pol()(barrier1.max())}),
                 sf::Color::Black);

  barrier2_line[0] =
      sf::Vertex(to_pixel({0.f, barrier2.pol()(0)}), sf::Color::Black);
  barrier2_line[1] =
      sf::Vertex(to_pixel({barrier2.max(), barrier2.pol()(barrier2.max())}),
                 sf::Color::Black);
}

void Plot::add(const std::vector<Point>& vertices)
{
  vertices_px.resize(vertices.size());
  std::transform(vertices.begin(), vertices.end(), vertices_px.begin(),
                 [this](const Point& p) {
                   return sf::Vertex(to_pixel(p), sf::Color::Red);
                 });
}

void Plot::show()
{
  sf::RenderWindow window(sf::VideoMode(w_, h_), "Biliardo triangolare");
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
    window.draw(barrier1_line);
    window.draw(barrier2_line);

    window.draw(&vertices_px[0], vertices_px.size(), sf::LineStrip);

    window.display();
  }
}
