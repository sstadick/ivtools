#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

class Interval {
public:
  int start;
  int stop;
  int val;

  inline bool overlaps(const int o_start, const int o_stop) {
    return start < o_stop && stop > o_start;
  }
};

// True if lhs is less than the rhs
bool iv_sorter(const Interval &lhs, const Interval &rhs) {
  if (lhs.start < rhs.start) {
    return true;
  } else if (lhs.start == rhs.start) {
    if (lhs.stop < rhs.stop) {
      return true;
    }
  }
  return false;
}

constexpr int infty = std::numeric_limits<int>::infinity();
typedef int int8_sim_t __attribute__((vector_size(8 * sizeof(int))));
static int8_sim_t *int8_sim_alloc(std::size_t n) {
  void *tmp = 0;
  if (posix_memalign(&tmp, sizeof(int8_sim_t), sizeof(int8_sim_t) * n)) {
    throw std::bad_alloc();
  }
  return (int8_sim_t *)tmp;
}
class Lapper {
public:
  std::vector<Interval> intervals;
  int max_len;
  int8_sim_t *lapper_starts;
  int8_sim_t *lapper_stops;
  int lapper_total_vecs;
  const int elem_per_vec = 8;

  Lapper(std::vector<Interval> ivs) {
    // Sort and set the intervals
    std::sort(ivs.begin(), ivs.end(), &iv_sorter);
    int max_len = 0;
    for (auto const &iv : ivs) {
      int len = iv.stop - iv.start;
      max_len = len > max_len ? len : max_len;
    }
    intervals = ivs;

    // Generate and fill the starts and stops
    // calculate the vec sizes
    int raw_vecs = intervals.size() / elem_per_vec;
    int mod = intervals.size() % elem_per_vec;
    int total_vecs = raw_vecs + (mod > 0 ? 1 : 0);

    std::cerr << "Total ves: " << total_vecs << std::endl;

    // starts - remainder
    int8_sim_t *starts = int8_sim_alloc(total_vecs);
    // stops - remainder
    int8_sim_t *stops = int8_sim_alloc(total_vecs);

    // fill the vectors
    for (int i = 0; i < total_vecs; ++i) {
      for (int j = 0; j < elem_per_vec; ++j) {
        if (i * j < intervals.size()) {
          starts[i][j] = intervals[i * elem_per_vec + j].start;
          stops[i][j] = intervals[i * elem_per_vec + j].stop;
        } else {
          starts[i][j] = infty;
          stops[i][j] = infty;
        }
      }
    }
    lapper_total_vecs = total_vecs;
    lapper_starts = starts;
    lapper_stops = stops;
  }

  ~Lapper() {
    std::free(lapper_starts);
    std::free(lapper_stops);
  }

  std::vector<std::shared_ptr<Interval>> find(const int f_start,
                                              const int f_stop) {
    int moded_start = f_start - max_len;
    int offset = lower_bound(moded_start > 0 ? moded_start : 0);
    std::vector<std::shared_ptr<Interval>> result;

    // search
    // Do a lower bounds search first, then start at the nearest full vector
    std::cerr << "Offset is: " << offset << std::endl;
    std::cerr << "Vec to start at is: " << offset / elem_per_vec << std::endl;
    bool break_early = false;
    for (int i = offset / elem_per_vec; i < lapper_total_vecs; ++i) {
      std::cerr << "i is: " << i << std::endl;
      int8_sim_t s_starts = lapper_starts[i];
      int8_sim_t s_stops = lapper_stops[i];

      // q(11, 15) vs i(10, 12)

      int8_sim_t start_v_stop = s_starts < f_stop; // 10 < 15 -> true
      int8_sim_t stop_v_start = s_stops > f_start; // 12 > 11 -> true

      int8_sim_t start_and_stop =
          start_v_stop && stop_v_start; // any both true are overlaps
      int8_sim_t start_gte_stop =
          s_starts >= f_stop; // check for break condition

      for (int b = 0; b < elem_per_vec; ++b) {
        std::cerr << b << " (" << s_starts[b] << ", " << s_stops[b] << ") -> "
                  << start_and_stop[b] << std::endl;
      }

      for (int j = 0; j < elem_per_vec; ++j) {
        if (start_and_stop[j]) {
          result.push_back(
              std::make_shared<Interval>(intervals[i * elem_per_vec + j]));
        }
        if (start_gte_stop[j]) {
          break_early = true;
          break;
        }
      }
      if (break_early) {
        break;
      }
    }
    return result;
  }

private:
  int lower_bound(int s) {
    int size = intervals.size();
    int low = 0;

    while (size > 0) {
      int half = size / 2;
      int other_half = size - half;
      int probe = low + half;
      int other_low = low + other_half;
      auto v = &intervals[probe];
      size = half;
      low = v->start < s ? other_low : low;
    }
    return low;
  }
};

int main() {
  std::vector<Interval> ivs = {
      Interval{.start = 0, .stop = 5, .val = 1},
      Interval{.start = 1, .stop = 5, .val = 2},
      Interval{.start = 4, .stop = 5, .val = 2},
      Interval{.start = 5, .stop = 9, .val = 4},
      Interval{.start = 0, .stop = 9, .val = 5},
      Interval{.start = 0, .stop = 5, .val = 6},
      Interval{.start = 0, .stop = 5, .val = 7},
      Interval{.start = 0, .stop = 2, .val = 8},
      Interval{.start = 0, .stop = 8, .val = 9},
      Interval{.start = 1, .stop = 9, .val = 10},
      Interval{.start = 3, .stop = 5, .val = 11},
      Interval{.start = 10, .stop = 12, .val = 12},
  };
  Lapper lapper(ivs);
  auto found = lapper.find(11, 15);
  for (const auto f : found) {
    std::cout << f->start << " " << f->stop << " " << f->val << std::endl;
  }
  return 1;
}
