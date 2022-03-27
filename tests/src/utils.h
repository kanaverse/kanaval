#ifndef UTILS_H
#define UTILS_H

#include "H5Cpp.h"

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

#endif
