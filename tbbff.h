#include <vector>
#include "tbb/pipeline.h"

namespace tbbff {

template<typename T>
using vec = std::vector<T>;

template<typename T, typename Iter>
tbb::filter_t<void, vec<T>* >
input_buffer(Iter &begin, Iter &end, size_t buffer_size) {
  return tbb::make_filter<void, vec<T>* >(
    tbb::filter::serial_out_of_order,
    [&begin, &end, buffer_size](tbb::flow_control &fc) -> vec<T>* {
      vec<T> *buffer = new vec<T>();
      buffer->reserve(buffer_size);
      if(begin == end) {
        fc.stop();
        return buffer;
      }
      size_t num_added = 0;
      while(begin != end && num_added < buffer_size) {
        ++num_added;
        buffer->push_back(*begin++);
      }
      return buffer;
    }
  );
}

template<typename T>
tbb::filter_t<vec<T>*, vec<T>* >
keep_if(std::function<bool(T &)> pred) {
  return tbb::make_filter<vec<T>*, vec<T>* >(
    tbb::filter::parallel,
    [&pred](vec<T>* buffer) -> vec<T>* {
      buffer->erase(std::remove_if(buffer->begin(), buffer->end(), std::not1(pred)),
                   buffer->end());
     return buffer;
    }
  );
}

template<typename T, typename Iter>
tbb::filter_t<vec<T>*, void> output_buffer(Iter &out) {
  return tbb::make_filter<vec<T>*, void>(
    tbb::filter::serial_out_of_order,
    [&out](vec<T>* buffer) -> void {
      std::move(buffer->begin(), buffer->end(), out);
      buffer->clear();
      delete buffer;
    }
  );
}

} // namespace tbbff


