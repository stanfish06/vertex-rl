#pragma once
#include "solver.h"
#include <filesystem>
#include <json-c/json.h>
#include <memory>
#include <string>

namespace vertex3d {
using Json = std::unique_ptr<json_object, decltype(&json_object_put)>;
using Path = std::filesystem::path;
Json parse(const std::string& text);
Json read_json(const Path& path);
json_object* field(json_object* object, const char* key, bool required = true);
void keys(json_object* object, std::initializer_list<const char*> allowed);
std::string string(json_object* value);
double number(json_object* value);
int integer(json_object* value, int minimum = 0);
unsigned seed_value(json_object* value);
Vec3 vector3(json_object* value);
double setting(json_object* object, const char* key, double fallback);
void zone_overrides(std::vector<Region>& regions, json_object* array);
void controls(std::vector<Region>& regions, json_object* array);
std::vector<Action> cell_actions(json_object* array);
Experiment experiment(json_object* object, const Path& base);
Json state(const Solver& solver, bool training);
}
