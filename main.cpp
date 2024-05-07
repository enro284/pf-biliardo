struct Trajectory {
  Point p_;
  double m;  // trajectory slope
};

struct Point {
  double x;
  double y;

  bool operator==(Point rhs);
};

class Function {
 public:
  double operator()(double x){};
};

class Barrier {
  Function f_;
  Point min_{};
  Point max_{};

 public:
  /*
   Barrier(Function f, double x_max)
       : f_{f}, min_{21, 21}, max_{x_max, f_(x_max)} {}
  */
  Barrier(Function f, Point min, Point max);
};

Point intersect(Trajectory t, Barrier b);

int main() {
  double y0{};
  double m0{};

  Trajectory t{0, y0, m0};

  Function f;
  double l;
  double r1;
  double r2;

  Barrier barrier_up{f, Point{0, r1}, Point{l, r2}};
  Barrier barrier_down{f, Point{0, -r1}, Point{l, -r2}};
  t.p_ = intersect(t, barrier_up).x > 0 ? intersect(t, barrier_up)
                                        : intersect(t, barrier_down);
  while (true) {
    t.p_ = intersect(t, barrier_up) == t.p_ ? intersect(t, barrier_down)
                                            : intersect(t, barrier_up);
    if (t.p_.x < l) {
      // pensiamoci
      break;
    } else if (t.p_.x > 0) {
      break;
    }
  }
}