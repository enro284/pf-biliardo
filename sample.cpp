#include "sample.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

int pf::Sample::size() const
{
  return entries_.size();
}

void pf::Sample::add(double x)
{
  entries_.push_back(x);
}

bool pf::Sample::remove(double x)
{
  for (auto it = entries_.begin(), end = entries_.end(); it != end; ++it) {
    if (*it == x) {
      entries_.erase(it);
      return true;
    }
  }
  return false;
}

pf::Statistics pf::Sample::statistics() const
{
  const int N{size()};

  if (N < 2) {
    throw std::runtime_error{"Not enough entries to run a statistics"};
  }

  double sum_x{};
  double sum_x2{};

  for (auto const& e : entries_) {
    sum_x += e;
    sum_x2 += e * e;
  }

  const double mean     = sum_x / N;
  const double sigma    = std::sqrt((sum_x2 - N * mean * mean) / (N - 1));
  const double mean_err = sigma / std::sqrt(N);

  auto entries = entries_;
  std::sort(entries.begin(), entries.end());

  const double median =
      N % 2 == 0 ? (entries[N / 2 - 1] + entries[N / 2]) * 0.5 : entries[N / 2];

  return {mean, sigma, mean_err, median};
}

pf::Sample& pf::Sample::operator+=(const pf::Sample& r)
{
  entries_.insert(entries_.end(), r.entries().begin(), r.entries().end());
  return *this;
}

pf::Sample pf::operator+(const pf::Sample& l, const pf::Sample& r)
{
  pf::Sample sum{l};
  return sum += r;
}
