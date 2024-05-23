#include <vector>

class Pol
{
  std::vector<double> coeff_;

 public:
  Pol(std::vector<double> coeff);
  double operator()(double x);
  double der(double x);

  int deg() const;
  std::vector<double> coeff() const;

  Pol operator-();
};

double eq_solve(Pol const& pol1, Pol const& pol2);
