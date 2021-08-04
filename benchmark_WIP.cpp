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

std::string file_to_string(std::string path) {
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

template <uint128_t prime_exp, uint128_t tau>
void benchmark_a(std::string &text, uint128_t base, std::string &path) {
  timer tmr{};
  auto rk = herlez::rolling_hash::rk_prime<decltype(text.cbegin()), tau, prime_exp>(text.cbegin(), text.cend(), base);
  uint128_t fp = rk.get_currect_fp();
  for (size_t i = 0; i < text.size() - tau; ++i) {
    fp = rk.roll();
  }
  size_t time = tmr.get();

  //Correctness test
  size_t last_window_index = text.size() - tau;
  auto rk_test = herlez::rolling_hash::rk_prime<decltype(text.cbegin()), tau, prime_exp>(text.cbegin() + last_window_index, text.cend(), base);
  uint128_t fp_test = rk_test.get_currect_fp();
  if (fp != fp_test) {
    std::cout << "Fingerprints were not calculated correctly.\n"
              << "Current: " << fp << " Correct: " << fp_test << "\n";
  }

  std::cout << "RESULT"
            << " algo=alex107"
            << " text=" << path
            << " size=" << text.size()
            << " time=" << time
            << " hash=" << fp
            << " correct=" << std::boolalpha << (fp == fp_test)
            << '\n';
}

/*
//Benchmark Jonas61
template <kr_fingerprinter::kr_fingerprinting::MersennePrime p, uint128_t tau>
void benchmark_j(std::string &text, uint128_t base, std::string &path) {
    timer tmr{};
    auto rk = kr_fingerprinting::kr_fingerprinter<p>::sliding_window_precompute<true>(tau, base);
    uint64_t fp = 0;
    for (size_t i = 0; i < tau; ++i) {
      fp = rk.roll_right(fp, 0, text[i]);
    }
    for (size_t i = 0; i < text.size() - tau; ++i) {
      fp = rk.roll_right(fp, text[i], text[i + tau]);
    }
    size_t time = tmr.get();

    //Correctness test
    size_t last_window_index = text.size() - tau;
    auto rk_test = kr_fingerprinting::kr_fingerprinter<p>::sliding_window_precompute<true>(tau, base);
    uint64_t fp_test = 0;
    for (size_t i = last_window_index; i < text.size(); ++i) {
      fp_test = rk_test.roll_right(fp_test, 0, text[i]);
    }
    if (fp != fp_test) {
      std::cout << "Fingerprints were not calculated correctly.\n"
                << "Current: " << fp << " Correct: " << fp_test << "\n";
    }

    std::cout << "RESULT"
              << " algo=jonas61"
              << " text=" << path
              << " size=" << text.size()
              << " time=" << time
              << " hash=" << fp
              << " correct=" << std::boolalpha << (fp == fp_test)
              << '\n';
  }
*/
  int main(int argc, char **argv) {
    //Window length
    static constexpr size_t tau = 256;
    size_t base;
    std::string text;
    std::string path;
    //Read File
    {
      if (argc != 3) {
        std::cerr << "Format: bench FILEPATH BASE\n";
        return -1;
      }
      path = argv[1];
      text = file_to_string(path);
      if (text.size() == 0) {
        std::cerr << "File " << path << " could not be read\n";
        return -1;
      }
      base = std::atoi(argv[2]);
    }

    benchmark_a<61, 1024>(text, base, path);
    benchmark_a<89, 1024>(text, base, path);
    benchmark_a<107, 1024>(text, base, path);

    //Benchmark Jonas61
    {
      timer tmr{};
      auto rk = kr_fingerprinting::kr_fingerprinter<kr_fingerprinting::MERSENNE61>::sliding_window_precompute<true>(tau, base);
      uint64_t fp = 0;
      for (size_t i = 0; i < tau; ++i) {
        fp = rk.roll_right(fp, 0, text[i]);
      }
      for (size_t i = 0; i < text.size() - tau; ++i) {
        fp = rk.roll_right(fp, text[i], text[i + tau]);
      }
      size_t time = tmr.get();

      //Correctness test
      size_t last_window_index = text.size() - tau;
      auto rk_test = kr_fingerprinting::kr_fingerprinter<kr_fingerprinting::MERSENNE61>::sliding_window_precompute<true>(tau, base);
      uint64_t fp_test = 0;
      for (size_t i = last_window_index; i < text.size(); ++i) {
        fp_test = rk_test.roll_right(fp_test, 0, text[i]);
      }
      if (fp != fp_test) {
        std::cout << "Fingerprints were not calculated correctly.\n"
                  << "Current: " << fp << " Correct: " << fp_test << "\n";
      }

      std::cout << "RESULT"
                << " algo=jonas61"
                << " text=" << path
                << " size=" << text.size()
                << " time=" << time
                << " hash=" << fp
                << " correct=" << std::boolalpha << (fp == fp_test)
                << '\n';
    }

    //Benchmark Jonas107
    {
      timer tmr{};
      auto rk = kr_fingerprinting::kr_fingerprinter<kr_fingerprinting::MERSENNE107>::sliding_window_precompute<true>(tau, base);
      uint128_t fp = 0;
      for (size_t i = 0; i < tau; ++i) {
        fp = rk.roll_right(fp, 0, text[i]);
      }
      for (size_t i = 0; i < text.size() - tau; ++i) {
        fp = rk.roll_right(fp, text[i], text[i + tau]);
      }
      size_t time = tmr.get();

      //Correctness test
      size_t last_window_index = text.size() - tau;
      auto rk_test = kr_fingerprinting::kr_fingerprinter<kr_fingerprinting::MERSENNE107>::sliding_window_precompute<true>(tau, base);
      uint128_t fp_test = 0;
      for (size_t i = last_window_index; i < text.size(); ++i) {
        fp_test = rk_test.roll_right(fp_test, 0, text[i]);
      }
      if (fp != fp_test) {
        std::cout << "Fingerprints were not calculated correctly.\n"
                  << "Current: " << fp << " Correct: " << fp_test << "\n";
      }

      std::cout << "RESULT"
                << " algo=jonas107"
                << " text=" << path
                << " size=" << text.size()
                << " time=" << time
                << " hash=" << fp
                << " correct=" << std::boolalpha << (fp == fp_test)
                << '\n';
    }

    //Benchmark Jonas127
    {
      timer tmr{};
      auto rk = kr_fingerprinting::kr_fingerprinter<kr_fingerprinting::MERSENNE127>::sliding_window(tau, base);
      uint128_t fp = 0;
      for (size_t i = 0; i < tau; ++i) {
        fp = rk.roll_right(fp, 0, text[i]);
      }
      for (size_t i = 0; i < text.size() - tau; ++i) {
        fp = rk.roll_right(fp, text[i], text[i + tau]);
      }
      size_t time = tmr.get();

      //Correctness test
      size_t last_window_index = text.size() - tau;
      auto rk_test = kr_fingerprinting::kr_fingerprinter<kr_fingerprinting::MERSENNE127>::sliding_window(tau, base);
      uint128_t fp_test = 0;
      for (size_t i = last_window_index; i < text.size(); ++i) {
        fp_test = rk_test.roll_right(fp_test, 0, text[i]);
      }
      if (fp != fp_test) {
        std::cout << "Fingerprints were not calculated correctly.\n"
                  << "Current: " << fp << " Correct: " << fp_test << "\n";
      }

      std::cout << "RESULT"
                << " algo=jonas127"
                << " text=" << path
                << " size=" << text.size()
                << " time=" << time
                << " hash=" << fp
                << " correct=" << std::boolalpha << (fp == fp_test)
                << '\n';
    }

   
    return 0;
  }
