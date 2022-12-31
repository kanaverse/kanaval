#ifndef KANAVAL_VALIDATE_HPP
#define KANAVAL_VALIDATE_HPP

#include "v2/_validate.hpp"
#include "v3/_validate.hpp"

/**
 * @file kanaval.hpp
 *
 * @brief Validate the embedded state file of a kana file.
 */

/** 
 * @namespace kanaval
 * @brief kana file validation 
 */
namespace kanaval{ 

/**
 * Validate the analysis state HDF5 file embedded inside a `*.kana` file.
 * An error is raised if an invalid structure is detected in any step.
 *
 * @param handle Open handle to a HDF5 file.
 * @param embedded Whether the data files are embedded.
 * @param version Version of the kana file.
 */
inline void validate(const H5::H5File& handle, bool embedded, int version) {
    if (version < 3000000) {
        v2::validate(handle, embedded, version);
    } else {
        v3::validate(handle, embedded, version);
    }
}

}

#endif
