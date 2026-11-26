#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

struct BoundaryVertex {
    float x, y, shoulder_x, shoulder_y;
};

struct Cell2D {
    int id, zone, dead, growing;
    float x, y, area, target_area;
    std::vector<BoundaryVertex> vertices;
};

struct Tissue2D {
    double time = 0;
    std::vector<Cell2D> cells;
};

Tissue2D measure_tissue(double time) {
    Tissue2D tissue{time, {}};
    for (int i = 1; i <= Nc; ++i) if (exist[3][i]) {
        Cell2D cell{i, get_subdomain(i), c_dead[i], c_growth[i],
                    c_center[i][1], c_center[i][2], c_area[i], c_A0[i], {}};
        cell.vertices.reserve(c_vertices[i][0]);
        for (int j = 1; j <= c_vertices[i][0]; ++j) {
            int id = c_vertices[i][j];
            cell.vertices.push_back({v[id][1], v[id][2], vS[id][1], vS[id][2]});
        }
        tissue.cells.push_back(std::move(cell));
    }
    return tissue;
}

class OutputFile {
    std::string path;
    std::unique_ptr<FILE, decltype(&fclose)> stream;
public:
    OutputFile(const std::string& filename, const char* mode = "w")
        : path("./output/" + filename), stream(fopen(path.c_str(), mode), fclose) {
        if (!stream) throw std::runtime_error("cannot open " + path);
    }
    FILE* get() const { return stream.get(); }
    void finish() {
        bool failed = ferror(stream.get());
        if (fclose(stream.release()) != 0 || failed) throw std::runtime_error("cannot write " + path);
    }
};

void export_tissue(const Tissue2D& tissue, int frame) {
    const std::string suffix = "_" + std::to_string(frame) + ".txt";
    OutputFile xs("X" + suffix), ys("Y" + suffix), sx("XS" + suffix), sy("YS" + suffix);
    OutputFile cx("XC" + suffix), cy("YC" + suffix), flags("C_dead" + suffix);
    OutputFile snapshot("snapshot" + suffix);
    fprintf(snapshot.get(), "cell_id,v_x,v_y,vs_x,vs_y,c_x,c_y,A,A0,dead\n");
    for (const Cell2D& cell : tissue.cells) {
        std::array<FILE*, 4> polygons{xs.get(), ys.get(), sx.get(), sy.get()};
        for (FILE* file : polygons) fprintf(file, "%zu ", cell.vertices.size());
        for (const auto& point : cell.vertices) {
            fprintf(xs.get(), "%.16g ", point.x);
            fprintf(ys.get(), "%.16g ", point.y);
            fprintf(sx.get(), "%.16g ", point.shoulder_x);
            fprintf(sy.get(), "%.16g ", point.shoulder_y);
            // The dev snapshot format stores the two shoulder columns in this order.
            fprintf(snapshot.get(), "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%d\n",
                    cell.id, point.x, point.y, point.shoulder_y, point.shoulder_x, cell.x,
                    cell.y, cell.area, cell.target_area, cell.dead);
        }
        for (FILE* file : polygons) {
            for (size_t j = cell.vertices.size(); j <= 15; ++j) fprintf(file, "0 ");
            fprintf(file, "\n");
        }
        fprintf(cx.get(), "%d %.16g \n", cell.id, cell.x);
        fprintf(cy.get(), "%d %.16g \n", cell.id, cell.y);
        fprintf(flags.get(), "%d %d \n", cell.id, cell.dead);
    }
    for (OutputFile* file : {&xs, &ys, &sx, &sy, &cx, &cy, &flags, &snapshot}) file->finish();
}

void export_scheduled_frame() {
    ++kount;
    export_tissue(measure_tissue(Time), kount);
    if (reference_output) {
        OutputFile index("frame_index.csv", kount == 1 ? "w" : "a");
        if (kount == 1) fprintf(index.get(), "frame_id,scheduled_time,actual_time\n");
        fprintf(index.get(), "%d,%.9g,%.9g\n", kount, TimeCount, Time);
        index.finish();
    }
}
