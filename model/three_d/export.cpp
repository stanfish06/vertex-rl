#include "export.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace vertex3d {
namespace {
template<class F>
void write_file(const std::filesystem::path& destination, F write) {
    auto temp = destination;
    temp += ".tmp";
    try {
        std::ofstream out;
        out.exceptions(std::ios::badbit | std::ios::failbit);
        out.open(temp);
        out << std::setprecision(17);
        write(out);
        out.close();
        std::filesystem::rename(temp,destination);
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(temp,ignored);
        throw;
    }
}

void frame_xml(std::ostream& out, const std::vector<Cell>& cells) {
    size_t points=0, faces=0;
    for (const auto& cell : cells) {
        points += cell.mesh.vertices.size();
        faces += cell.mesh.faces.size();
    }
    out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<PolyData><Piece NumberOfPoints=\"" << points << "\" NumberOfPolys=\"" << faces
        << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\">\n"
        << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& c : cells) for (auto p : c.mesh.vertices) out << p.x << ' ' << p.y << ' ' << p.z << '\n';
    out << "</DataArray></Points>\n<Polys><DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    size_t base = 0;
    for (const auto& c : cells) {
        for (auto f : c.mesh.faces) out << base+f[0] << ' ' << base+f[1] << ' ' << base+f[2] << '\n';
        base += c.mesh.vertices.size();
    }
    out << "</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i=1; i<=faces; ++i) out << 3*i << ' ';
    out << "\n</DataArray></Polys>\n<CellData Scalars=\"cell_id\">\n";
    auto array = [&](const char* name, const char* type, auto get) {
        out << "<DataArray type=\"" << type << "\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (const auto& c : cells) {
            for (size_t i=0; i<c.mesh.faces.size(); ++i) out << get(c) << ' ';
            out << '\n';
        }
        out << "</DataArray>\n";
    };
    array("cell_id","Int32",[](const Cell& c){return c.id;});
    array("zone_id","Int32",[](const Cell& c){return c.zone;});
    array("volume","Float64",[](const Cell& c){return c.geometry.volume;});
    out << "</CellData></Piece></PolyData></VTKFile>\n";
}
}

void Collection::write(const Solver& solver, int frame) {
    if (frame < 1) throw std::invalid_argument("frame ID must be positive");
    auto found = frames_.find(frame);
    if (found != frames_.end() && found->second != solver.time()) throw std::invalid_argument("frame ID already used");
    std::filesystem::create_directories(directory_);
    std::string filename = "frame_" + std::to_string(frame) + ".vtp";
    write_file(directory_/filename,[&](auto& out){frame_xml(out,solver.cells());});
    auto next = frames_;
    next[frame] = solver.time();
    std::vector<std::pair<double,int>> ordered;
    for (auto [id,time] : next) ordered.emplace_back(time,id);
    std::sort(ordered.begin(),ordered.end());
    write_file(directory_/"tissue.pvd",[&](auto& out) {
        out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n";
        for (auto [time,id] : ordered) {
            out << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"frame_" << id << ".vtp\"/>\n";
        }
        out << "</Collection></VTKFile>\n";
    });
    frames_.swap(next);
}
}
