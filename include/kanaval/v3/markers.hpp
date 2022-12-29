#ifndef KANAVAL_MARKERS_V3_HPP
#define KANAVAL_MARKERS_V3_HPP

#include "H5Cpp.h"
#include <vector>

namespace kanaval {

namespace v3 {

namespace markers {
    
inline const std::vector<std::string> effects { "lfc", "delta_detected", "cohen", "auc" };

inline bool check_common_parameters(const H5::Group& phandle) {
    auto lfc_threshold = utils::load_float_scalar<>(phandle, "lfc_threshold");
    if (lfc_threshold < 0) {
        throw std::runtime_error("'lfc_threshold' must be non-negative");
    }
    return utils::load_integer_scalar<>(phandle, "compute_auc");
}

}

}

}

#endif 
