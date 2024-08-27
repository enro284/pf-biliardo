#ifndef KINEMATICS_HPP
#define KINEMATICS_HPP

#include "mathematics.hpp"
#include <iostream>

struct Point
{
  double x_;
  double y_;

  bool operator==(Point const& rhs) const;
};

class Result
{
  Point p_;
  double theta_;

 public:
  // normal constructors
  Result(double x, double y, double theta);
  Result(Point const& p, double theta);

  // constructs invalid result (x = -1)
  Result();

  bool operator==(Result b) const;

  double get_x() const;
  double get_y() const;
  double get_theta() const;

  bool is_valid() const;
};
std::ostream& operator<<(std::ostream& os, Result const& res);

struct Trajectory
{
  Point p_;
  double m_; 
  bool up_;

  void exit(double x);
  Result result() const;
};

class Barrier
{
  Point max_;
  Pol pol_;

 public:
  Barrier(Pol p, double x_max);

  double max() const;
  Pol pol() const;
};

Point intersect(Trajectory t, Barrier b);

Result simulate_single_particle(Barrier const& barrier_up,
                                Barrier const& barrier_down, Point p0,
                                double m0);

#endif