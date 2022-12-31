#ifndef KANAVAL__METADATA_V3_HPP
#define KANAVAL__METADATA_V3_HPP

#include "H5Cpp.h"
#include "../utils.hpp"
#include <stdexcept>

namespace kanaval {

namespace v3 {

// Here we're using a double underscore for consistency,
// given that everything else has validate_<step_name>.
inline void validate__metadata(const H5::H5File& handle, int version) {
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
