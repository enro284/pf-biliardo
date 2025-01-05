#ifndef STATISTICS_HPP
#define STATISTICS_HPP

struct Statistics
{
  double mean{};
  double std_dev{};
  double skewness{};
  double kurtosis{}; // implemented as excess kurtosis
};

class Sample
{
  int n;
  double s1;
  double s2;
  double s3;
  double s4;

 public:
  Sample();
  
  void add(double x);
  int size();

  Statistics statistics() const;
};
#endif