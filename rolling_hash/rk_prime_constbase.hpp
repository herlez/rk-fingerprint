#pragma once
#include <array>
#include <cmath>
#include <iterator>

namespace herlez::rolling_hash {
typedef size_t t_fp;
__extension__ typedef unsigned __int128 uint128_t;
static constexpr uint128_t m_prime = (1ULL << 61) - 1;


template <typename t_it, size_t t_tau = 1024>
class rk_prime {
  static constexpr uint64_t calculatePowerModulo(const uint64_t power,
                                                 const uint128_t kPrime) {
    uint128_t x = m_base;
    for (unsigned int i = 0; i < power; i++) {
      x = (x * x) % kPrime;
    }
    return static_cast<uint64_t>(x);
  }

  std::array<uint64_t, 256 * 256> calculatePowerTable() {
    std::array<uint64_t, 256 * 256> char_influence;
    for (size_t i = 0; i < 256; ++i) {
      for (size_t j = 0; j < 256; ++j) {
        const uint128_t first_char_influence = i * m_two_pow_tau_mod_q;
        const uint128_t last_char_influence = j; // * m_base ^ 0 = 1;
        const uint128_t char_pair_influence = m_prime - (first_char_influence > last_char_influence ? mod_m_prime(first_char_influence - last_char_influence)
                                        : mod_m_prime(m_prime - (last_char_influence - first_char_influence)));
        char_influence[256 * i + j] = char_pair_influence;
      }
    }
    return char_influence;
  }

 public:
  rk_prime(t_it fp_begin, t_it text_end)
      : m_fp_begin(fp_begin),
        m_fp_end(fp_begin + t_tau),
        m_text_end(text_end),
        m_cur_fp(0) {
    m_char_influence = calculatePowerTable();
    //Calculate first window
    for (uint64_t i = 0; i < t_tau; ++i) {
      m_cur_fp *= m_base;
      m_cur_fp += (unsigned char)m_fp_begin[i];
      m_cur_fp = mod_m_prime(m_cur_fp);
    }
  }

  t_fp roll_safe(size_t i = 1) {
    m_cur_fp *= m_base;
    m_cur_fp += (unsigned char)*m_fp_end;
    m_cur_fp = mod_m_prime(m_cur_fp);

    uint128_t first_char_influence = m_char_influence[*m_fp_begin];
    m_cur_fp = first_char_influence < m_cur_fp ? m_cur_fp - first_char_influence : m_prime - (first_char_influence - m_cur_fp);

    std::advance(m_fp_begin, 1);
    if (m_fp_end != m_text_end) { //safety check
      std::advance(m_fp_end, 1);
    }
    return m_cur_fp;
  }

  t_fp roll([[maybe_unused]] size_t i = 1) {
    
    // m_cur_fp *= m_base;
    // m_cur_fp += (unsigned char)*m_fp_end;
    // m_cur_fp = mod_m_prime(m_cur_fp);
    // uint128_t first_char_influence = m_char_influence[256 * *m_fp_begin];
    // m_cur_fp = first_char_influence < m_cur_fp ? m_cur_fp - first_char_influence : m_prime - (first_char_influence - m_cur_fp);
    
    m_cur_fp *= m_base;
    uint128_t border_char_influence = m_char_influence[256 * *m_fp_begin + *m_fp_end];
    m_cur_fp += border_char_influence;
    m_cur_fp = mod_m_prime(m_cur_fp);

    std::advance(m_fp_begin, 1);
    std::advance(m_fp_end, 1);
    return m_cur_fp;
  }

  t_fp get_currect_fp() {
    return static_cast<t_fp>(m_cur_fp);
  }

  inline uint128_t mod_m_prime(uint128_t num) const {
    num = (num & m_prime) + (num >> 61);
    num = (num > m_prime) ? (num - m_prime) : num;
    return num;
  }

 private:
  t_it m_fp_begin;
  t_it m_fp_end;
  t_it m_text_end;
  uint128_t m_cur_fp;

  static constexpr uint64_t m_base = 256;
  static constexpr uint64_t m_two_pow_tau_mod_q = calculatePowerModulo(std::log2(t_tau), m_prime);
  std::array<uint64_t, 256 * 256> m_char_influence;
};
}  // namespace herlez::rolling_hash