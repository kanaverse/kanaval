#ifndef KANAVAL_CELL_FILTERING_V2_HPP
#define KANAVAL_CELL_FILTERING_V2_HPP

#include "H5Cpp.h"
#include <vector>
#include "../utils.hpp"
#include "misc.hpp"

namespace kanaval {

namespace v2 {

inline int validate_cell_filtering(const H5::H5File& handle, int num_cells, int num_modalities, int version) {
    if (version < 2000000) { // didn't exist before v2.
        return -1;
    }

    auto qhandle = utils::check_and_open_group(handle, "cell_filtering");
    utils::check_and_open_group(qhandle, "parameters");

    int remaining = -1;
    try {
        auto rhandle = utils::check_and_open_group(qhandle, "results");
        if (num_modalities > 1) {
            remaining = quality_control::check_discard_vector(rhandle, num_cells, false);
        }
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to retrieve results from 'cell_filtering'");
    }

    return remaining;
}

}

}

#endif
