// Closed, individually owned triangle meshes in simulation length units.
#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace vertex3d {

struct Vec3 {
    double x=0, y=0, z=0;
    double& operator[](std::size_t i) { return i==0 ? x : i==1 ? y : z; }
    const double& operator[](std::size_t i) const { return i==0 ? x : i==1 ? y : z; }
};
inline Vec3 operator+(Vec3 a,Vec3 b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
inline Vec3 operator-(Vec3 a,Vec3 b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
inline Vec3 operator-(Vec3 a) { return {-a.x,-a.y,-a.z}; }
inline Vec3 operator*(Vec3 a,double b) { return {a.x*b,a.y*b,a.z*b}; }
inline Vec3 operator*(double a,Vec3 b) { return b*a; }
inline Vec3 operator/(Vec3 a,double b) { return {a.x/b,a.y/b,a.z/b}; }
inline Vec3& operator+=(Vec3& a,Vec3 b) { return a=a+b; }
inline Vec3& operator-=(Vec3& a,Vec3 b) { return a=a-b; }
inline double dot(Vec3 a,Vec3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Vec3 cross(Vec3 a,Vec3 b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
inline double norm(Vec3 a) { return std::sqrt(dot(a,a)); }

using Face=std::array<int,3>;
using Edge=std::pair<int,int>;
struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<Face> faces;
    // Faces traversing the edge first in ascending, then descending vertex order.
    std::map<Edge,std::array<int,2>> edges;
};
struct Measures { double volume, area; Vec3 centroid; };
Measures measure(const Mesh&);
void validate(Mesh&);  // checks a single embedded, outward-wound, closed manifold
void validate_geometry(const Mesh&);  // same connectivity, new vertex positions
std::vector<Vec3> internal_forces(const Mesh&,double target_volume,double stiffness,double tension);
// First daughter is on the negative side. Invalid cuts throw and preserve the input.
std::pair<Mesh,Mesh> split(const Mesh&,Vec3 plane_point,Vec3 normal);

Mesh box(Vec3 lo,Vec3 hi);
Mesh sphere(Vec3 center,double radius,double max_error);
Mesh ellipsoid(Vec3 center,Vec3 radii,double max_error);
Mesh cylinder(Vec3 center,double radius,double height,double max_error);  // axis +z
Mesh import_wall(const std::string& path);  // OBJ or STL, validated as supplied
int locate_point(const Mesh& closed,Vec3 point);  // exterior -1, boundary 0, interior +1
bool contains_point(const Mesh& closed,Vec3 point);  // includes boundary
bool contains(const Mesh& wall,const Mesh& cell);    // wholly in the strict interior

struct Zone { int id; std::function<bool(Vec3)> matches; };
Zone layer(int id,Vec3 origin,Vec3 axis,double lower,double upper);
Zone radial_band(int id,Vec3 center,double lower,double upper);
Zone surface_zone(int id,Mesh closed_surface,bool inside);
Zone interface_zone(int id,const std::string& path,bool positive_side);
int classify(Vec3 centroid,const std::vector<Zone>& ordered);
} // namespace vertex3d
