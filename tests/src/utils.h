#ifndef UTILS_H
#define UTILS_H

#include "H5Cpp.h"

static constexpr int latest = 2001000;

template<class Object>
void quick_write_dataset(Object& handle, std::string name, int val) {
    H5::DataSpace space;
    auto dhandle = handle.createDataSet(name, H5::PredType::NATIVE_INT, space);
    dhandle.write(&val, H5::PredType::NATIVE_INT);
    return;
}

template<class Object>
void quick_write_dataset(Object& handle, std::string name, double val) {
    H5::DataSpace space;
    auto dhandle = handle.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
    dhandle.write(&val, H5::PredType::NATIVE_DOUBLE);
    return;
}

template<class Object>
void quick_write_dataset(Object& handle, std::string name, std::string val) {
    H5::DataSpace space;
    H5::StrType stype(0, H5T_VARIABLE);
    auto dhandle = handle.createDataSet(name, stype, space);
    dhandle.write(val, stype);
    return;
}

inline H5::DataSpace create_space(hsize_t n) {
    H5::DataSpace space;
    space.setExtentSimple(1, &n);
    return space;
}

template<class Object>
void quick_write_dataset(Object& handle, std::string name, std::vector<int> val) {
    auto space = create_space(val.size());
    auto dhandle = handle.createDataSet(name, H5::PredType::NATIVE_INT, space);
    dhandle.write(val.data(), H5::PredType::NATIVE_INT);
    return;
}

template<class Object>
void quick_write_dataset(Object& handle, std::string name, std::vector<double> val) {
    auto space = create_space(val.size());
    auto dhandle = handle.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
    dhandle.write(val.data(), H5::PredType::NATIVE_DOUBLE);
    return;
}

template<class Object>
void quick_write_dataset(Object& handle, std::string name, std::vector<std::string> val) {
    size_t maxlen = 0;
    for (auto v : val) {
        if (v.size() > maxlen) {
            maxlen = v.size();
        }
    }

    std::vector<char> buffer(maxlen * val.size());
    auto bIt = buffer.begin();
    for (size_t i = 0; i < val.size(); ++i, bIt += maxlen) {
        std::copy(val[i].begin(), val[i].end(), bIt);
    }

    H5::StrType stype(H5::PredType::C_S1, maxlen);
    auto space = create_space(val.size());
    auto dhandle = handle.createDataSet(name, stype, space);
    dhandle.write(buffer.data(), stype);
    return;
}


template<class Function>
void quick_throw(Function fun, std::string msg) {
    EXPECT_ANY_THROW({
        try {
            fun();
            std::cout << "failed to throw '" << msg << "'" << std::endl;
        } catch (std::exception& e) {
            std::string found(e.what());
            if (found.find(msg) == std::string::npos) {
                std::cout << "error '" << e.what() << "' does not match '" << msg << "'" << std::endl;
                EXPECT_FALSE(true);
            }
            throw e;
        }
    });
}

#endif
