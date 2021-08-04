#pragma once
#include <assert.h>

#include <chrono>
#include <cmath>
#include <iterator>

namespace herlez::rolling_hash {
__extension__ typedef unsigned __int128 uint128_t;

template <typename t_it, size_t t_prime_exp = 107>
class rk_prime {
 public:
  rk_prime(t_it fp_begin, uint128_t tau, uint128_t base)
      : m_tau(tau),
        m_fp_begin(fp_begin),
        m_fp_end(fp_begin + m_tau),
        m_cur_fp(0),
        m_base(base) {
    //prime should be mersenne and (prime*base + prime) should not overflow
    //static_assert(t_prime_exp == 107 || t_prime_exp == 61 || t_prime_exp == 89);
    assert(t_prime_exp + std::log2(m_base) < 126);
    fillPowerTable();
    //Calculate first window
    for (size_t i = 0; i < m_tau; ++i) {
      m_cur_fp *= m_base;
      m_cur_fp += (unsigned char)m_fp_begin[i];
      m_cur_fp = mod_m_prime(m_cur_fp);
    }
  }

  inline uint128_t roll([[maybe_unused]] size_t i = 1) {
    m_cur_fp *= m_base;
    uint128_t border_char_influence = m_char_influence[(unsigned char)(*m_fp_begin)][(unsigned char)(*m_fp_end)];
    m_cur_fp += border_char_influence;
    m_cur_fp = mod_m_prime(m_cur_fp);

    std::advance(m_fp_begin, 1);
    std::advance(m_fp_end, 1);
    return m_cur_fp;
  }

  inline uint128_t get_currect_fp() const {
    return m_cur_fp;
  }

 private:
  static constexpr uint128_t m_prime = (uint128_t{1} << t_prime_exp) - 1;
  uint128_t m_tau;
  t_it m_fp_begin;
  t_it m_fp_end;
  uint128_t m_cur_fp;

  uint128_t m_base;
  uint128_t m_char_influence[256][256];

  inline uint128_t mod_m_prime(uint128_t num) const {
    //Does only work for 2^127 - 1
    // uint128_t const z = (num + 1) >> t_prime_exp;
    // return (num + z) & m_prime;


    // uint128_t const v = num + 1;
    // uint64_t const z = ((v >> t_prime_exp) + v) >> t_prime_exp;
    // return (num + z) & m_prime;

    num = (num & m_prime) + (num >> t_prime_exp);
    return (num >= m_prime) ? (num - m_prime) : num;
  }

  // To compute (a * b) % mod
  inline uint128_t mulmod(uint128_t a, uint128_t b) const {
    uint128_t res = 0;  // Initialize result
    a = a % m_prime;
    while (b > 0) {
      // If b is odd, add 'a' to result
      if (b % 2 == 1)
        res = (res + a) % m_prime;

      // Multiply 'a' with 2
      a = (a * 2) % m_prime;

      // Divide b by 2
      b /= 2;
    }
    // Return result
    return res % m_prime;
  }

  inline uint128_t calculatePowerModulo() const {
    assert(__builtin_popcount(m_tau) == 1);
    uint128_t x = m_base;
    for (unsigned int i = 0; i < std::log2(m_tau); i++) {
      x = mulmod(x, x);
      //x = mod_m_prime(x*x); //overflows
    }
    return x;
  }

  inline void fillPowerTable() {
    //std::chrono::system_clock::time_point begin_ = std::chrono::system_clock::now();
    const uint128_t m_two_pow_tau_mod_q = calculatePowerModulo();
    for (size_t i = 0; i < 256; ++i) {
      m_char_influence[i][0] = m_prime - (mod_m_prime(m_two_pow_tau_mod_q * i));
      for (size_t j = 1; j < 256; ++j) {
        m_char_influence[i][j] = mod_m_prime(m_char_influence[i][j - 1] + 1);
      }
    }
    //std::chrono::system_clock::time_point const end = std::chrono::system_clock::now();
    //std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - begin_).count() << '\n';
  }
};
}  // namespace herlez::rolling_hash