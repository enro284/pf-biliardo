#include <vector>
class Function
{
 public:
};

class Pol : Function
{
  std::vector<double> coeff_;

 public:
  double operator()(double x_);
  double der(double x_);
};

double intersect(Function t, Function b);
