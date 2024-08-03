#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <vector>

namespace pf {

struct Statistics
{
  double mean{};
  double sigma{};
  double mean_err{};
  double median{};
  // aggiungere coeff (a)simmetria:
  // https://en.wikipedia.org/wiki/Skewness#Fisher's_moment_coefficient_of_skewness
  // aggiungere coeff appiattimento:
  // https://statorials.org/it/coefficiente-di-appiattimento/
};

class Sample
{
  std::vector<double> entries_;

 public:
  int size() const;

  void add(double x);
  bool remove(double x);

  auto& entries()
  {
    return entries_;
  }
  const auto& entries() const
  {
    return entries_;
  }

  Statistics statistics() const;

  Sample& operator+=(const Sample& r);
};

Sample operator+(const Sample& l, const Sample& r);

} // namespace pf

#endif
