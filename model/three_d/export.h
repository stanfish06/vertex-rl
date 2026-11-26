#pragma once
#include "solver.h"
#include <filesystem>
#include <map>

namespace vertex3d {
class Collection {
    std::filesystem::path directory_;
    std::map<int,double> frames_;
public:
    explicit Collection(std::filesystem::path directory) : directory_(std::move(directory)) {}
    void write(const Solver& solver, int frame);
};
}
