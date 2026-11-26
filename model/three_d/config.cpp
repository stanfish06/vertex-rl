#include "config.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <set>
#include <stdexcept>

namespace vertex3d {

Json parse(const std::string& text) {
    std::unique_ptr<json_tokener, decltype(&json_tokener_free)> tok(json_tokener_new(), json_tokener_free);
    json_tokener_set_flags(tok.get(), JSON_TOKENER_STRICT);
    Json result(json_tokener_parse_ex(tok.get(), text.c_str(), text.size()+1), json_object_put);
    if (!result || json_tokener_get_error(tok.get()) != json_tokener_success) {
        throw std::invalid_argument("invalid JSON");
    }
    return result;
}

Json read_json(const Path& path) {
    std::ifstream stream(path);
    if (!stream) throw std::invalid_argument("cannot open " + path.string());
    std::ostringstream text;
    text << stream.rdbuf();
    return parse(text.str());
}

void keys(json_object* object, std::initializer_list<const char*> allowed) {
    if (!json_object_is_type(object, json_type_object)) throw std::invalid_argument("expected an object");
    json_object_object_foreach(object, key, value) {
        (void)value;
        if (std::none_of(allowed.begin(), allowed.end(), [&](auto name) { return std::strcmp(key, name) == 0; })) {
            throw std::invalid_argument("unknown field " + std::string(key));
        }
    }
}

json_object* field(json_object* object, const char* key, bool required) {
    json_object* value = nullptr;
    if (!json_object_object_get_ex(object, key, &value)) {
        if (required) throw std::invalid_argument("missing field " + std::string(key));
    } else if (!value) throw std::invalid_argument("null field " + std::string(key));
    return value;
}

std::string string(json_object* value) {
    if (!json_object_is_type(value, json_type_string)) throw std::invalid_argument("expected a string");
    return std::string(json_object_get_string(value), json_object_get_string_len(value));
}

double number(json_object* value) {
    if (!json_object_is_type(value,json_type_int) && !json_object_is_type(value,json_type_double)) {
        throw std::invalid_argument("expected a number");
    }
    double x = json_object_get_double(value);
    if (!std::isfinite(x)) throw std::invalid_argument("expected a finite number");
    return x;
}

int integer(json_object* value, int minimum) {
    if (!json_object_is_type(value,json_type_int)) throw std::invalid_argument("expected an integer");
    auto x = json_object_get_int64(value);
    if (x < minimum || x > std::numeric_limits<int>::max()) throw std::invalid_argument("integer out of range");
    return static_cast<int>(x);
}

unsigned seed_value(json_object* value) {
    if (!json_object_is_type(value,json_type_int)) throw std::invalid_argument("seed must be an unsigned integer");
    auto x = json_object_get_int64(value);
    if (x < 0 || x > UINT32_MAX) throw std::invalid_argument("seed out of range");
    return static_cast<unsigned>(x);
}

static bool flag(json_object* object, const char* key, bool fallback) {
    auto* value = field(object, key, false);
    if (!value) return fallback;
    if (!json_object_is_type(value,json_type_boolean)) throw std::invalid_argument(std::string(key)+" must be boolean");
    return json_object_get_boolean(value);
}

double setting(json_object* object, const char* key, double fallback) {
    auto* value = field(object, key, false);
    return value ? number(value) : fallback;
}

static size_t count(json_object* array) {
    if (!json_object_is_type(array,json_type_array)) throw std::invalid_argument("expected an array");
    return json_object_array_length(array);
}

Vec3 vector3(json_object* value) {
    if (count(value) != 3) throw std::invalid_argument("expected three coordinates");
    return {number(json_object_array_get_idx(value,0)), number(json_object_array_get_idx(value,1)),
            number(json_object_array_get_idx(value,2))};
}

static Mesh mesh(json_object* value, const Path& base) {
    auto kind = string(field(value,"shape"));
    auto vec = [&](const char* key) { return vector3(field(value,key)); };
    auto num = [&](const char* key) { return number(field(value,key)); };
    if (kind == "box") {
        keys(value,{"shape","min","max"});
        return box(vec("min"), vec("max"));
    }
    if (kind == "sphere") {
        keys(value,{"shape","center","radius","max_surface_error"});
        return sphere(vec("center"), num("radius"), num("max_surface_error"));
    }
    if (kind == "ellipsoid") {
        keys(value,{"shape","center","radii","max_surface_error"});
        return ellipsoid(vec("center"), vec("radii"), num("max_surface_error"));
    }
    if (kind == "cylinder") {
        keys(value,{"shape","center","radius","height","max_surface_error"});
        return cylinder(vec("center"), num("radius"), num("height"), num("max_surface_error"));
    }
    if (kind == "import") {
        keys(value,{"shape","path"});
        return import_wall((base/string(field(value,"path"))).string());
    }
    if (kind == "mesh") {
        keys(value,{"shape","vertices","faces"});
        Mesh result;
        auto* points = field(value,"vertices");
        auto* faces = field(value,"faces");
        for (size_t i=0; i<count(points); ++i) result.vertices.push_back(vector3(json_object_array_get_idx(points,i)));
        for (size_t i=0; i<count(faces); ++i) {
            auto* face = json_object_array_get_idx(faces,i);
            if (count(face) != 3) throw std::invalid_argument("expected a triangle");
            result.faces.push_back({integer(json_object_array_get_idx(face,0)),integer(json_object_array_get_idx(face,1)),
                                    integer(json_object_array_get_idx(face,2))});
        }
        return result;
    }
    throw std::invalid_argument("unknown mesh shape " + kind);
}

static Zone zone(json_object* spec, int id, const Path& base) {
    auto kind = string(field(spec,"kind"));
    double inf = std::numeric_limits<double>::infinity();
    if (kind == "all") {
        keys(spec,{"kind"});
        return Zone{id, [](Vec3) { return true; }};
    }
    if (kind == "layer") {
        keys(spec,{"kind","origin","axis","lower","upper"});
        return layer(id, vector3(field(spec,"origin")), vector3(field(spec,"axis")),
                     setting(spec,"lower",-inf), setting(spec,"upper",inf));
    }
    if (kind == "radial") {
        keys(spec,{"kind","center","lower","upper"});
        return radial_band(id, vector3(field(spec,"center")), setting(spec,"lower",0), setting(spec,"upper",inf));
    }
    if (kind == "surface") {
        keys(spec,{"kind","mesh","inside"});
        return surface_zone(id,mesh(field(spec,"mesh"),base),flag(spec,"inside",true));
    }
    if (kind == "interface") {
        keys(spec,{"kind","path","positive_side"});
        return interface_zone(id,(base/string(field(spec,"path"))).string(),flag(spec,"positive_side",true));
    }
    throw std::invalid_argument("unknown zone rule " + kind);
}

static void zone_settings(Region& r, json_object* z) {
    r.allow_growth = flag(z,"allow_growth",r.allow_growth);
    r.allow_stop = flag(z,"allow_stop",r.allow_stop);
    r.allow_death = flag(z,"allow_death",r.allow_death);
    r.allow_division = flag(z,"allow_division",r.allow_division);
    r.growth_rate = setting(z,"growth_rate",r.growth_rate);
    r.target_size = setting(z,"target_size",r.target_size);
    r.death_rate = setting(z,"death_rate",r.death_rate);
    r.division_volume = setting(z,"division_volume",r.division_volume);
    r.division_rate = setting(z,"division_rate",r.division_rate);
}

static Region& region_by_id(std::vector<Region>& regions, int id) {
    auto r = std::find_if(regions.begin(), regions.end(), [id](const Region& r){ return r.zone.id == id; });
    if (r == regions.end()) throw std::invalid_argument("unknown zone " + std::to_string(id));
    return *r;
}

void zone_overrides(std::vector<Region>& regions, json_object* array) {
    if (regions.empty()) regions.push_back({Zone{0, [](Vec3) { return true; }}});
    std::set<int> seen;
    for (size_t i=0; i<count(array); ++i) {
        auto* z = json_object_array_get_idx(array,i);
        keys(z,{"id","allow_growth","allow_stop","allow_death","allow_division","growth_rate","target_size",
                "division_rate","division_volume","death_rate"});
        int id = integer(field(z,"id"));
        if (!seen.insert(id).second) throw std::invalid_argument("duplicate zone");
        zone_settings(region_by_id(regions,id),z);
    }
}

void controls(std::vector<Region>& regions, json_object* array) {
    std::set<std::pair<int,std::string>> seen;
    for (size_t i=0; i<count(array); ++i) {
        auto* control = json_object_array_get_idx(array,i);
        keys(control,{"zone","parameter","value"});
        int id = integer(field(control,"zone"));
        auto& r = region_by_id(regions,id);
        auto name = string(field(control,"parameter"));
        if (!seen.emplace(id,name).second) throw std::invalid_argument("duplicate control " + name);
        double value = number(field(control,"value"));
        if (value < 0 || (name == "target_size" && value == 0)) throw std::invalid_argument("invalid " + name);
        if (name == "growth_rate") r.growth_rate = value;
        else if (name == "target_size") r.target_size = value;
        else if (name == "division_rate") r.division_rate = value;
        else if (name == "death_rate") r.death_rate = value;
        else throw std::invalid_argument("unknown control " + name);
    }
}

std::vector<Action> cell_actions(json_object* array) {
    std::vector<Action> result;
    for (size_t i=0; i<count(array); ++i) {
        auto* a = json_object_array_get_idx(array,i);
        int id = integer(field(a,"id"),1);
        auto name = string(field(a,"action"));
        if (name == "divide") {
            keys(a,{"id","action","orientation"});
            Vec3 normal = vector3(field(a,"orientation"));
            if (!(norm(normal)>0) || !std::isfinite(norm(normal))) throw std::invalid_argument("invalid plane normal");
            result.push_back({id,Fate::divide,normal});
        } else {
            keys(a,{"id","action"});
            if (name == "grow") result.push_back({id,Fate::grow});
            else if (name == "stop_growing") result.push_back({id,Fate::stop});
            else if (name == "die") result.push_back({id,Fate::die});
            else throw std::invalid_argument("unknown fate " + name);
        }
    }
    return result;
}

Experiment experiment(json_object* object, const Path& base) {
    keys(object,{"domain","cells","zones","mechanics","seed","time_limit","output_interval","divisions"});
    Experiment result;
    result.wall = mesh(field(object,"domain"),base);
    auto* cells = field(object,"cells");
    for (size_t i=0; i<count(cells); ++i) {
        auto* spec = json_object_array_get_idx(cells,i);
        keys(spec,{"mesh","target_volume","growing"});
        Cell c;
        c.mesh = mesh(field(spec,"mesh"),base);
        if (auto* target = field(spec,"target_volume",false)) {
            c.target_volume = number(target);
            if (c.target_volume <= 0) throw std::invalid_argument("target_volume must be positive");
        }
        c.growing = flag(spec,"growing",true);
        result.cells.push_back(std::move(c));
    }
    if (auto* zones = field(object,"zones",false)) {
        for (size_t i=0; i<count(zones); ++i) {
            auto* z = json_object_array_get_idx(zones,i);
            keys(z,{"id","rule","allow_growth","allow_division","growth_rate","division_volume","division_rate",
                    "allow_stop","allow_death","death_rate","target_size"});
            Region r{zone(field(z,"rule"), integer(field(z,"id")), base)};
            zone_settings(r,z);
            result.regions.push_back(std::move(r));
        }
    }
    if (auto* obj = field(object,"mechanics",false)) {
        keys(obj,{"volume_stiffness","surface_tension","pressure","contact_range","bin_width",
                  "max_displacement","max_step","min_step","growth_half_size","remesh_level","contact_barrier","smoothing"});
        auto& m = result.material;
        m.volume_stiffness = setting(obj,"volume_stiffness",m.volume_stiffness);
        m.surface_tension = setting(obj,"surface_tension",m.surface_tension);
        m.pressure = setting(obj,"pressure",m.pressure);
        m.contact_range = setting(obj,"contact_range",m.contact_range);
        m.bin_width = setting(obj,"bin_width",m.bin_width);
        m.max_displacement = setting(obj,"max_displacement",m.max_displacement);
        m.max_step = setting(obj,"max_step",m.max_step);
        m.growth_half_size = setting(obj,"growth_half_size",m.growth_half_size);
        m.remesh_level = static_cast<int>(setting(obj,"remesh_level",m.remesh_level));
        m.contact_barrier = flag(obj,"contact_barrier",m.contact_barrier);
        m.smoothing = setting(obj,"smoothing",m.smoothing);
        m.min_step = setting(obj,"min_step",m.min_step);
    }
    return result;
}

static json_object* point(Vec3 p) {
    auto* a=json_object_new_array_ext(3);
    for (double x : {p.x,p.y,p.z}) json_object_array_add(a,json_object_new_double(x));
    return a;
}

Json state(const Solver& solver, bool training) {
    Json result(json_object_new_object(),json_object_put);
    json_object_object_add(result.get(),"time",json_object_new_double(solver.time()));
    auto* cells = json_object_new_array();
    json_object_object_add(result.get(),"cells",cells);
    auto* m = json_object_new_object();
    json_object_object_add(result.get(),"measurements",m);
    std::vector<int> counts(solver.regions().size());
    double total = 0;
    for (const auto& c : solver.cells()) {
        auto* obj = json_object_new_object();
        json_object_array_add(cells,obj);
        json_object_object_add(obj,"id",json_object_new_int(c.id));
        json_object_object_add(obj,"zone",json_object_new_int(c.zone));
        if (training) json_object_object_add(obj,"age",json_object_new_double(solver.time()-c.born));
        json_object_object_add(obj,"growing",json_object_new_boolean(c.growing));
        json_object_object_add(obj,"center",point(c.geometry.centroid));
        json_object_object_add(obj,"volume",json_object_new_double(c.geometry.volume));
        json_object_object_add(obj,"surface_area",json_object_new_double(c.geometry.area));
        json_object_object_add(obj,"target_volume",json_object_new_double(c.target_volume));
        auto* vertices = json_object_new_array();
        auto* faces = json_object_new_array();
        json_object_object_add(obj,"vertices",vertices);
        json_object_object_add(obj,"faces",faces);
        for (auto p : c.mesh.vertices) json_object_array_add(vertices,point(p));
        for (auto f : c.mesh.faces) {
            auto* triangle = json_object_new_array_ext(3);
            for (int j : f) json_object_array_add(triangle,json_object_new_int(j));
            json_object_array_add(faces,triangle);
        }
        total += c.geometry.volume;
        for (size_t i=0;i<counts.size();++i) if (solver.regions()[i].zone.id==c.zone) ++counts[i];
    }
    json_object_object_add(m,"cell_count",json_object_new_int64(solver.cells().size()));
    json_object_object_add(m,"total_volume",json_object_new_double(total));
    auto* zcounts = json_object_new_array();
    json_object_object_add(m,"zone_counts",zcounts);
    for (int n : counts) json_object_array_add(zcounts,json_object_new_int(n));
    json_object_object_add(result.get(),"accepted_steps",json_object_new_uint64(solver.accepted_steps()));
    json_object_object_add(result.get(),"rejected_trials",json_object_new_uint64(solver.retries()));
    return result;
}
}
