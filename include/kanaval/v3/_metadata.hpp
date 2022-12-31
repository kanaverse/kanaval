#ifndef KANAVAL_METADATA_V3_HPP
#define KANAVAL_METADATA_V3_HPP

#include "H5Cpp.h"
#include "../utils.hpp"
#include <stdexcept>

namespace kanaval {

namespace v3 {

void validate_metadata(const H5::H5File& handle, int version) {
    try {    
        auto mhandle = utils::check_and_open_group(handle, "_metadata");
        auto v = utils::load_integer_scalar(mhandle, "format_version");
        if (v != version) {
            throw std::runtime_error("'format_version' is not consistent with kana file version");
        }

        utils::load_string(mhandle, "application_name");
        utils::load_string(mhandle, "application_version");
    } catch (std::exception& e) {
        throw utils::combine_errors(e, "failed to check the '_metadata'");
    }
}

}

}

#endif
