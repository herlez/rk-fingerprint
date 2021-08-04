#include <fstream>
#include <iostream>
#include <limits>
#include <streambuf>
#include <string>
#include <vector>

#include "rolling_hash/fingerprinting.hpp"
#include "rolling_hash/rk_prime.hpp"
#include "timer.hpp"

__extension__ typedef unsigned __int128 uint128_t;

inline std::string file_to_string(std::string path) {
  std::ifstream stream(path);
  std::string text((std::istreambuf_iterator<char>(stream)),
                   std::istreambuf_iterator<char>());
  return text;
}

inline std::ostream &operator<<(std::ostream &out,
                                uint128_t value) {
  std::ostream::sentry s(out);
  if (s) {
    char buffer[64];  // 39 should be enough
    char *digit = &(buffer[64]);
    do {
      *(--digit) = "0123456789"[value % 10];
      value /= 10;
    } while (value != 0);
    int len = &(buffer[64]) - digit;
    if (out.rdbuf()->sputn(digit, len) != len) {
      out.setstate(std::ios_base::badbit);
    }
  }
  return out;
}

template <uint128_t prime_exp>
inline void benchmark_a(std::string &text, std::string &path, uint128_t tau, uint128_t base) {
  timer tmr{};
  auto rk = herlez::rolling_hash::rk_prime<decltype(text.cbegin()), prime_exp>(text.cbegin(), tau, base);
  auto fp = rk.get_currect_fp();
  for (size_t i = 0; i < text.size() - tau; ++i) {
    fp = rk.roll();
  }
  size_t time = tmr.get();

  //Correctness test
  size_t last_window_index = text.size() - tau;
  auto rk_test = herlez::rolling_hash::rk_prime<decltype(text.cbegin()), prime_exp>(text.cbegin() + last_window_index, tau, base);
  auto fp_test = rk_test.get_currect_fp();
  if (fp != fp_test) {
    std::cout << "Fingerprints were not calculated correctly.\n"
              << "Current: " << fp << " Correct: " << fp_test << "\n";
  }

  std::cout << "RESULT"
            << " algo=alex" << prime_exp
            << " text=" << path
            << " size=" << text.size()
            << " time=" << time
            << " hash=" << fp
            << " correct=" << std::boolalpha << (fp == fp_test)
            << '\n';
}

template <kr_fingerprinting::MersennePrime p>
inline void benchmark_j(std::string &text, std::string &path, uint128_t tau, uint128_t base) {
  timer tmr{};
  typename kr_fingerprinting::kr_fingerprinter<p>::sliding_window_precompute<true> rk = typename kr_fingerprinting::kr_fingerprinter<p>::sliding_window_precompute<true>(tau, base);
  typename kr_fingerprinting::kr_fingerprinter<p>::uintX_t fp = 0;
  for (size_t i = 0; i < tau; ++i) {
    fp = rk.roll_right(fp, (u_char)0, (u_char)text[i]);
  }
  for (size_t i = 0; i < text.size() - tau; ++i) {
    fp = rk.roll_right(fp, (u_char)text[i], (u_char)text[i + tau]);
  }
  size_t time = tmr.get();

  //Correctness test
  size_t last_window_index = text.size() - tau;
  typename kr_fingerprinting::kr_fingerprinter<p>::sliding_window_precompute<true> rk_test = typename kr_fingerprinting::kr_fingerprinter<p>::sliding_window_precompute<true>(tau, base);
  typename kr_fingerprinting::kr_fingerprinter<p>::uintX_t fp_test = 0;
  for (size_t i = last_window_index; i < text.size(); ++i) {
    fp_test = rk_test.roll_right(fp_test, (u_char)0, (u_char)text[i]);
  }
  if (fp != fp_test) {
    std::cout << "Fingerprints were not calculated correctly.\n"
              << "Current: " << fp << " Correct: " << fp_test << "\n";
  }

  std::cout << "RESULT"
            << " algo=jonas" << p::shift
            << " text=" << path
            << " size=" << text.size()
            << " time=" << time
            << " hash=" << fp
            << " correct=" << std::boolalpha << (fp == fp_test)
            << '\n';
}

int main(int argc, char **argv) {
  //Window length
  size_t tau;
  std::string text;
  std::string path;
  size_t base;
  //Read File
  {
    if (argc != 4) {
      std::cerr << "Format: bench FILEPATH TAU BASE\n";
      return -1;
    }
    path = argv[1];
    text = file_to_string(path);
    if (text.size() == 0) {
      std::cerr << "File " << path << " could not be read\n";
      return -1;
    }
    tau = std::atoi(argv[2]);
    base = std::atoi(argv[3]);
  }

  benchmark_a<61>(text, path, tau, base);
  benchmark_a<89>(text, path, tau, base);
  benchmark_a<107>(text, path, tau, base);

  benchmark_j<kr_fingerprinting::MERSENNE61>(text, path, tau, base);
  benchmark_j<kr_fingerprinting::MERSENNE89>(text, path, tau, base);
  benchmark_j<kr_fingerprinting::MERSENNE107>(text, path, tau, base);
  benchmark_j<kr_fingerprinting::MERSENNE127>(text, path, tau, base);

  //WHY IS THIS SO FAST?
  constexpr uint128_t tau1 = 107;
  {
    timer tmr{};
    auto rk = herlez::rolling_hash::rk_prime<decltype(text.cbegin()), tau1>(text.cbegin(), tau, base);
    auto fp = rk.get_currect_fp();
    for (size_t i = 0; i < text.size() - tau; ++i) {
      fp = rk.roll();
    }
    size_t time = tmr.get();

    //Correctness test
    size_t last_window_index = text.size() - tau;
    auto rk_test = herlez::rolling_hash::rk_prime<decltype(text.cbegin()), tau1>(text.cbegin() + last_window_index, tau, base);
    auto fp_test = rk_test.get_currect_fp();
    if (fp != fp_test) {
      std::cout << "Fingerprints were not calculated correctly.\n"
                << "Current: " << fp << " Correct: " << fp_test << "\n";
    }

    std::cout << "RESULT"
              << " algo=alex" << tau1
              << " text=" << path
              << " size=" << text.size()
              << " time=" << time
              << " hash=" << fp
              << " correct=" << std::boolalpha << (fp == fp_test)
              << '\n';
  }
  return 0;
}
