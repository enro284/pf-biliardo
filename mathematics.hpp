#include <vector>
class Function
{
 public:
  virtual double operator()(double x) = 0;
};

class Pol : Function
{
  std::vector<double> coeff_;

 public:
  Pol(std::vector<double> coeff);
  double operator()(double x);
  double der(double x);
};

// double intersect(Function t, Function b);
