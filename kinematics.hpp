#ifndef KINEMATICS_HPP
#define KINEMATICS_HPP

#include "mathematics.hpp"
#include <iostream>
#include <vector>

class Result
{
  Vec2 p_;
  double theta_;

 public:
  // see Trajectory.result() instead
  Result(Vec2 const& p, double theta);

  bool operator==(Result b) const;

  double get_x() const;
  double get_y() const;
  double get_theta() const;
};
std::ostream& operator<<(std::ostream& os, Result const& res);

struct Trajectory
{
  Vec2 p_;
  Vec2 v_;

  Trajectory(Vec2 const& p_, double theta);
  void exit(double x); // check for x=0

  Result result() const;
};

class Barrier
{
  Vec2 max_;
  Pol pol_;

 public:
  // generic constructor
  Barrier(Pol const& p, double x_max);
  // constructor for linear barrier
  Barrier(double l, double r1, double r2);

  double max() const; // TODO: x_max
  Pol pol() const;
};

struct Bounce
{
  Vec2 p_;
  Barrier const* b_ptr; // TODO check pointer safety
};

// return all possible collisions (going the right way, in barrier bounds,
// different from current trajectory point)
std::vector<Bounce> intersect(Trajectory const& t, Barrier const* b);

Result simulate_single_particle(Barrier const& barrier_up,
                                Barrier const& barrier_down, Trajectory t,
                                std::vector<Vec2>& points);

Result simulate_single_particle(Barrier const& barrier_up,
                                Barrier const& barrier_down, Trajectory t);

#endif