#ifndef KANAVAL_UTILS_HPP
#define KANAVAL_UTILS_HPP
#include "H5Cpp.h"

namespace kanaval {

namespace utils {

template<class Object>
H5::Group check_and_open_group(const Object& handle, const std::string& name) {
    if (!handle.exists(name) || handle.childObjType(name) != H5O_TYPE_GROUP) {
        throw std::runtime_error("'" + name + "' group does not exist");
    }
    return handle.openGroup(name);
}

template<class Object>
H5::DataSet check_and_open_dataset(const Object& handle, const std::string& name) {
    if (!handle.exists(name) || handle.childObjType(name) != H5O_TYPE_DATASET) {
        throw std::runtime_error("'" + name + "' dataset does not exist");
    }
    return handle.openDataSet(name);
}

template<class Object>
H5::DataSet check_and_open_dataset(const Object& handle, const std::string& name, H5T_class_t expected_type) {
    auto dhandle = check_and_open_dataset(handle, name);    
    if (dhandle.getTypeClass() != expected_type) {
        std::string expected = "other";
        if (expected_type == H5T_INTEGER) {
            expected = "integer";
        } else if (expected_type == H5T_STRING) {
            expected = "string";
        } else if (expected_type == H5T_FLOAT) {
            expected = "float";
        }
        throw std::runtime_error("'" + name + "' dataset should be of type " + expected);
    }
    return dhandle;
}

template<class Object>
H5::DataSet check_and_open_dataset(const Object& handle, const std::string& name, H5T_class_t expected_type, const std::vector<size_t>& expected_dims) {
    auto dhandle = check_and_open_dataset(handle, name, expected_type);
    auto dspace = dhandle.getSpace();

    size_t ndims = dspace.getSimpleExtentNdims();
    if (ndims != expected_dims.size()) {
        throw std::runtime_error("'" + name + "' dataset does not have the expected dimensions");
    }

    if (ndims) {
        std::vector<hsize_t> observed(ndims);
        dspace.getSimpleExtentDims(observed.data());
        for (size_t i = 0; i < ndims; ++i) {
            if (observed[i] != expected_dims[i]) {
                throw std::runtime_error("'" + name + "' dataset does not have the expected dimensions");
            }
        }
    }

    return dhandle;
}

template<class Object>
H5::DataSet check_and_open_scalar(const Object& handle, const std::string& name, H5T_class_t expected_type) {
    auto dhandle = check_and_open_dataset(handle, name, expected_type);
    auto dspace = dhandle.getSpace();

    size_t ndims = dspace.getSimpleExtentNdims();
    if (ndims) {
        throw std::runtime_error("'" + name + "' dataset should be a scalar");
    }
    return dhandle;
}

template<typename T = int, class Object>
T load_integer_scalar(const Object&handle, const std::string& name) {
    auto dhandle = check_and_open_scalar(handle, name, H5T_INTEGER);

    // TODO: support more int types.
    T output;
    dhandle.read(&output, H5::PredType::NATIVE_INT);

    return output;    
}

template<typename T = double, class Object>
T load_float_scalar(const Object&handle, const std::string& name) {
    auto dhandle = check_and_open_scalar(handle, name, H5T_FLOAT);

    // TODO: support more types.
    T output;
    dhandle.read(&output, H5::PredType::NATIVE_DOUBLE);

    return output;    
}

inline std::runtime_error combine_errors(const std::exception& e, const std::string& msg) {
    return std::runtime_error(msg + "\n  - " + std::string(e.what()));
}

}

}

#endif
