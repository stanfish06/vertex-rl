#include "mesh.h"
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace vertex3d {
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double roundoff = 128 * std::numeric_limits<double>::epsilon();

void require(bool ok, const char* reason) { if (!ok) throw std::invalid_argument(reason); }
bool finite(Vec3 v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }
bool positive(double v) { return std::isfinite(v) && v > 0; }
Vec3 unit(Vec3 v) {
    double d = norm(v);
    require(positive(d), "expected a finite nonzero direction");
    return v/d;
}
Edge edge(int a,int b) { return std::minmax(a,b); }

struct Bounds { Vec3 lo, hi; };
Bounds bounds(const std::vector<Vec3>& points) {
    require(!points.empty(), "empty surface");
    Bounds b{points[0],points[0]};
    for (auto v: points) {
        require(finite(v), "nonfinite vertex");
        for (int k=0;k<3;++k) {
            b.lo[k]=std::min(b.lo[k],v[k]);
            b.hi[k]=std::max(b.hi[k],v[k]);
        }
    }
    return b;
}
double tolerance(const Mesh& m) {
    const auto b=bounds(m.vertices);
    return roundoff*norm(b.hi-b.lo);
}
using Triangle=std::array<Vec3,3>;
Triangle triangle(const Mesh& m,Face f) {
    return {m.vertices[f[0]],m.vertices[f[1]],m.vertices[f[2]]};
}
Bounds tri_bounds(const Triangle& t) {
    Bounds b{t[0],t[0]};
    for (int k=0;k<3;++k) {
        b.lo[k]=std::min({t[0][k],t[1][k],t[2][k]});
        b.hi[k]=std::max({t[0][k],t[1][k],t[2][k]});
    }
    return b;
}
bool overlap(Bounds a,Bounds b,double eps) {
    for (int k=0;k<3;++k)
        if (a.hi[k]+eps<b.lo[k] || b.hi[k]+eps<a.lo[k]) return false;
    return true;
}

// Inside the triangle's three edge halfspaces on its plane.
bool in_triangle(Vec3 p,const Triangle& t,Vec3 n,double eps) {
    for (int i=0;i<3;++i) {
        auto e=t[(i+1)%3]-t[i];
        if (dot(cross(e,p-t[i]),n)<-eps*norm(e)) return false;
    }
    return true;
}

// Report all extremes of the segment/triangle intersection, including coplanar overlap.
template<class Report>
bool segment_hits(Vec3 a,Vec3 b,const Triangle& t,Vec3 n,double eps,Report report) {
    double da=dot(a-t[0],n), db=dot(b-t[0],n);
    if (std::abs(da)<=eps && std::abs(db)<=eps) {
        double lo=0,hi=1;
        for (int i=0;i<3;++i) {
            auto e=t[(i+1)%3]-t[i];
            double u=dot(cross(e,a-t[i]),n), v=dot(cross(e,b-t[i]),n);
            double error=eps*norm(e);
            if (std::abs(u)<=error) u=0;
            if (std::abs(v)<=error) v=0;
            if (u<0 && v<0) return false;
            if (u<0) lo=std::max(lo,u/(u-v));
            if (v<0) hi=std::min(hi,u/(u-v));
        }
        return lo<=hi && (report(a+(b-a)*lo) || report(a+(b-a)*hi));
    }
    if ((da>eps && db>eps) || (da< -eps && db< -eps)) return false;
    double u=std::clamp(da/(da-db),0.,1.);
    auto p=a+(b-a)*u;
    return in_triangle(p,t,n,eps) && report(p);
}
template<class Report>
bool intersections(const Triangle& a,const Triangle& b,double eps,Report report) {
    auto na=unit(cross(a[1]-a[0],a[2]-a[0]));
    auto nb=unit(cross(b[1]-b[0],b[2]-b[0]));
    for (int k=0;k<3;++k) {
        if (segment_hits(a[k],a[(k+1)%3],b,nb,eps,report) ||
            segment_hits(b[k],b[(k+1)%3],a,na,eps,report)) return true;
    }
    return false;
}

// Pairs sharing an index may meet only on that vertex or shared edge.
void embedded(const Mesh& m,double eps) {
    std::vector<Bounds> boxes;
    for (auto f: m.faces) boxes.push_back(tri_bounds(triangle(m,f)));
    std::vector<int> order(m.faces.size());
    std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](int a,int b){ return boxes[a].lo.x<boxes[b].lo.x; });
    for (size_t i=0;i<order.size();++i) {
        int a=order[i];
        for (size_t j=i+1;j<order.size();++j) {
            int b=order[j];
            if (boxes[b].lo.x>boxes[a].hi.x+eps) break;
            if (!overlap(boxes[a],boxes[b],eps)) continue;
            auto ta=triangle(m,m.faces[a]),tb=triangle(m,m.faces[b]);
            auto origin=ta[0];
            for (auto& v:ta) v-=origin;
            for (auto& v:tb) v-=origin;
            std::vector<Vec3> shared;
            for (int k: m.faces[a])
                for (int l: m.faces[b]) if (k==l) shared.push_back(m.vertices[k]-origin);
            auto unexpected=[&](Vec3 p) {
                if (shared.empty()) return true;
                auto d=p-shared[0];
                if (shared.size()==1) return norm(d)>4*eps;
                auto e=shared[1]-shared[0];
                return norm(d-e*std::clamp(dot(d,e)/dot(e,e),0.,1.))>4*eps;
            };
            require(!intersections(ta,tb,eps,unexpected),
                    "self-intersecting surface");
        }
    }
}

void topology(Mesh& m,bool closed) {
    double eps=tolerance(m);
    std::map<Edge,std::array<int,2>> edges;
    std::vector<std::vector<int>> incident(m.vertices.size());
    require(!m.faces.empty(), "empty face list");
    for (int i=0;i<int(m.faces.size());++i) {
        Face f=m.faces[i];
        for (int v:f) {
            require(v>=0 && v<int(m.vertices.size()), "invalid face index");
            incident[v].push_back(i);
        }
        auto t=triangle(m,f);
        auto u=t[1]-t[0], v=t[2]-t[0];
        require(norm(cross(u,v))>eps*std::max(norm(u),norm(v)), "degenerate triangle");
        for (int j=0;j<3;++j) {
            int a=f[j], b=f[(j+1)%3];
            auto& sides=edges.try_emplace(edge(a,b),std::array<int,2>{-1,-1}).first->second;
            int& side=sides[a>b];
            require(side<0, "edge has inconsistent winding or more than two faces");
            side=i;
        }
    }
    for (const auto& entry: edges)
        if (closed) require(entry.second[0]>=0 && entry.second[1]>=0,"open mesh boundary");

    // A manifold vertex has one connected fan. Edge incidence alone misses pinched vertices.
    std::vector<int> visit(m.faces.size(),-1);
    for (int v=0;v<int(incident.size());++v) {
        require(!incident[v].empty(),"unused vertex");
        std::vector<int> pending{incident[v][0]};
        visit[pending[0]]=v;
        size_t reached=0;
        while (!pending.empty()) {
            int i=pending.back(); pending.pop_back(); ++reached;
            for (int j: m.faces[i]) if (j!=v)
                for (int n:edges.at(edge(j,v))) if (n>=0 && visit[n]!=v) {
                    visit[n]=v; pending.push_back(n);
                }
        }
        require(reached==incident[v].size(), "nonmanifold vertex");
    }
    std::vector<bool> visited(m.faces.size());
    std::vector<int> pending{0};
    visited[0]=true;
    size_t count=0;
    while (!pending.empty()) {
        int i=pending.back(); pending.pop_back(); ++count;
        auto f=m.faces[i];
        for (int k=0;k<3;++k)
            for (int n:edges.at(edge(f[k],f[(k+1)%3]))) if (n>=0 && !visited[n]) {
                visited[n]=true; pending.push_back(n);
            }
    }
    require(count==m.faces.size(),"disconnected surface");
    if (closed) require(positive(measure(m).volume), "surface requires positive outward volume");
    m.edges=std::move(edges);
}
} // namespace

Measures measure(const Mesh& mesh) {
    require(!mesh.vertices.empty(),"empty mesh");
    Vec3 origin=mesh.vertices[0];
    std::array<long double,3> moment{};
    long double volume6=0,area=0;
    for (auto f:mesh.faces) {
        auto t=triangle(mesh,f);
        Vec3 a=t[0]-origin, b=t[1]-origin, c=t[2]-origin;
        double dv6=dot(a,cross(b,c));
        volume6+=dv6;
        area+=norm(cross(b-a,c-a))*0.5;
        for (int k=0;k<3;++k)
            moment[k]+=(static_cast<long double>(a[k])+b[k]+c[k])*dv6;
    }
    double volume=double(volume6/6);
    require(volume!=0 && std::isfinite(volume) && std::isfinite(area), "invalid mesh measure");
    Vec3 center=origin;
    for (int k=0;k<3;++k) center[k]+=double(moment[k]/(4*volume6));
    require(finite(center), "invalid volume centroid");
    return {volume,double(area),center};
}
void validate(Mesh& m) {
    topology(m,true);
    embedded(m,tolerance(m));
}
std::vector<Vec3> internal_forces(const Mesh& m,double target,double stiffness,double tension) {
    require(positive(target) && std::isfinite(stiffness) && stiffness>=0 &&
            std::isfinite(tension) && tension>=0, "invalid mechanics coefficient or target volume");
    double pressure=-stiffness*(measure(m).volume-target);
    auto origin=m.vertices[0];
    std::vector<Vec3> result(m.vertices.size());
    for (auto f:m.faces) {
        auto t=triangle(m,f);
        Vec3 n=unit(cross(t[1]-t[0],t[2]-t[0]));
        for (int k=0;k<3;++k) {
            Vec3 b=t[(k+1)%3],c=t[(k+2)%3];
            result[f[k]]+=cross(b-origin,c-origin)*(pressure/6)-cross(b-c,n)*(tension*0.5);
        }
    }
    return result;
}

namespace {
void triangulate(const std::vector<Vec3>& vertices,std::vector<int> polygon,std::vector<Face>& faces) {
    require(polygon.size()>=3,"polygon has fewer than three vertices");
    Vec3 p=vertices.at(polygon[0]), sum{};
    double extent=0;
    for (size_t i=0;i<polygon.size();++i) {
        auto a=vertices.at(polygon[i])-p, b=vertices.at(polygon[(i+1)%polygon.size()])-p;
        sum+=cross(a,b);
        extent=std::max(extent,norm(a));
    }
    Vec3 n=unit(sum);
    double eps=extent*roundoff;
    for (int i:polygon)
        require(std::abs(dot(vertices[i]-p,n))<=eps, "nonplanar polygon");
    while (polygon.size()>3) {
        bool clipped=false;
        for (size_t j=0;j<polygon.size();++j) {
            int a=polygon[(j+polygon.size()-1)%polygon.size()], b=polygon[j], c=polygon[(j+1)%polygon.size()];
            Triangle t{vertices[a],vertices[b],vertices[c]};
            if (dot(cross(t[1]-t[0],t[2]-t[0]),n)<=eps*extent) continue;
            bool blocked=false;
            for (int v:polygon) {
                if (v!=a && v!=b && v!=c && in_triangle(vertices[v],t,n,eps)) {
                    blocked=true; break;
                }
            }
            if (blocked) continue;
            faces.push_back({a,b,c});
            polygon.erase(polygon.begin()+j);
            clipped=true;
            break;
        }
        require(clipped,"cannot triangulate polygon boundary");
    }
    faces.push_back({polygon[0],polygon[1],polygon[2]});
}

Mesh compact(const std::vector<Vec3>& vertices,std::vector<Face> faces) {
    Mesh m;
    std::map<int,int> ids;
    for (auto& f:faces) for (int& v:f) {
        auto entry=ids.try_emplace(v,int(m.vertices.size()));
        if (entry.second) m.vertices.push_back(vertices[v]);
        v=entry.first->second;
    }
    m.faces=std::move(faces);
    return m;
}
} // namespace

std::pair<Mesh,Mesh> split(const Mesh& parent,Vec3 point,Vec3 normal) {
    auto n=unit(normal);
    require(finite(point),"nonfinite division point");
    // Clip and cap in coordinates local to the cell. Rounding a small
    // cut point into large world coordinates before triangulation can spoil planarity.
    auto vertices=parent.vertices;
    Vec3 origin=parent.vertices.at(0);
    double offset=dot(n,point-origin);
    double eps=tolerance(parent);
    std::vector<double> distance;
    for (auto& p:vertices) {
        double input_error=std::numeric_limits<double>::epsilon()*
            (std::abs(n.x*p.x)+std::abs(n.y*p.y)+std::abs(n.z*p.z));
        p-=origin;
        double d=dot(n,p)-offset;
        distance.push_back(std::abs(d)<=eps+input_error ? 0. : d);
    }
    auto extremes=std::minmax_element(distance.begin(),distance.end());
    require(*extremes.first<0 && *extremes.second>0,"division plane does not cut the cell");
    std::map<Edge,int> intersections;
    auto intersection=[&](int a,int b) {
        auto e=edge(a,b);
        auto entry=intersections.try_emplace(e,int(vertices.size()));
        if (entry.second) {
            double fraction=distance[a]/(distance[a]-distance[b]);
            vertices.push_back(vertices[a]+(vertices[b]-vertices[a])*fraction);
        }
        return entry.first->second;
    };
    std::array<std::vector<Face>,2> sides;
    for (auto f:parent.faces) {
        require(distance[f[0]]!=0 || distance[f[1]]!=0 || distance[f[2]]!=0,
                "division plane coincides with a surface face");
        for (int side=0;side<2;++side) {
            std::vector<int> polygon;
            for (int k=0;k<3;++k) {
                int a=f[k], b=f[(k+1)%3];
                if ((side==0 && distance[a]<=0) || (side==1 && distance[a]>=0)) polygon.push_back(a);
                if ((distance[a]<0 && distance[b]>0) || (distance[a]>0 && distance[b]<0))
                    polygon.push_back(intersection(a,b));
            }
            for (size_t j=2;j<polygon.size();++j)
                sides[side].push_back({polygon[0],polygon[j-1],polygon[j]});
        }
    }

    // Unpaired directed edges are precisely the boundary of the clipped negative half.
    std::set<Edge> boundary;
    for (auto f:sides[0]) for (int k=0;k<3;++k) {
        int a=f[k],b=f[(k+1)%3];
        if (!boundary.erase({b,a}))
            require(boundary.emplace(a,b).second,"ambiguous cut boundary");
    }
    std::map<int,int> next;
    for (auto e:boundary)
        require(next.emplace(e.second,e.first).second,"cut boundary branches");
    require(!next.empty(),"division plane has no cap");
    while (!next.empty()) {
        int start=next.begin()->first, i=start;
        std::vector<int> loop;
        do {
            loop.push_back(i);
            auto it=next.find(i);
            require(it!=next.end(),"cut boundary is open");
            i=it->second;
            next.erase(it);
        } while (i!=start);
        std::vector<Face> cap;
        triangulate(vertices,loop,cap);
        auto first=cap.front();
        require(dot(cross(vertices[first[1]]-vertices[first[0]],vertices[first[2]]-vertices[first[0]]),n)>0,
                "division section has a hole and cannot form a simple cap");
        for (auto f:cap) {
            sides[0].push_back(f);
            sides[1].push_back({f[2],f[1],f[0]});
        }
    }
    std::copy(parent.vertices.begin(),parent.vertices.end(),vertices.begin());
    for (size_t i=parent.vertices.size();i<vertices.size();++i) vertices[i]+=origin;
    auto negative=compact(vertices,std::move(sides[0])), positive=compact(vertices,std::move(sides[1]));
    validate(negative);
    validate(positive);
    auto volume=measure(parent).volume;
    require(std::abs(measure(negative).volume+measure(positive).volume-volume)<=1e-10*volume,
            "daughter volumes do not conserve parent volume");
    return {std::move(negative),std::move(positive)};
}

Mesh box(Vec3 lo,Vec3 hi) {
    require(finite(lo) && finite(hi) && lo.x<hi.x && lo.y<hi.y && lo.z<hi.z,"invalid box bounds");
    Mesh m;
    m.vertices={{lo.x,lo.y,lo.z},{hi.x,lo.y,lo.z},{hi.x,hi.y,lo.z},{lo.x,hi.y,lo.z},
                {lo.x,lo.y,hi.z},{hi.x,lo.y,hi.z},{hi.x,hi.y,hi.z},{lo.x,hi.y,hi.z}};
    m.faces={{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},
             {1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
    topology(m,true);
    return m;
}
Mesh ellipsoid(Vec3 center,Vec3 radii,double max_error) {
    require(finite(center) && positive(radii.x) && positive(radii.y) && positive(radii.z) &&
            positive(max_error),"invalid ellipsoid radii, center or surface error");
    Mesh m;
    m.vertices={{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0,0,1},{0,0,-1}};
    for (int i=0;i<4;++i) {
        m.faces.push_back({i,(i+1)%4,4});
        m.faces.push_back({(i+1)%4,i,5});
    }
    double largest=std::max({radii.x,radii.y,radii.z});
    for (;;) {
        // Unit-sphere radial projection differs by <= 1 - face-plane distance
        // throughout each facet. The affine ellipsoid map stretches by at most largest.
        bool sufficient=true;
        for (auto f:m.faces) {
            auto t=triangle(m,f);
            double h=dot(t[0],unit(cross(t[1]-t[0],t[2]-t[0])));
            if (largest*(1-h)>max_error) { sufficient=false; break; }
        }
        if (sufficient) break;
        require(m.faces.size()<=2000000/4,"requested ellipsoid surface exceeds two million facets");
        std::map<Edge,int> cache;
        auto midpoint=[&](int a,int b) {
            auto entry=cache.try_emplace(edge(a,b),int(m.vertices.size()));
            if (entry.second) m.vertices.push_back(unit(m.vertices[a]+m.vertices[b]));
            return entry.first->second;
        };
        std::vector<Face> refined;
        refined.reserve(4*m.faces.size());
        for (auto f:m.faces) {
            int a=midpoint(f[0],f[1]),b=midpoint(f[1],f[2]),c=midpoint(f[2],f[0]);
            refined.insert(refined.end(),{{f[0],a,c},{f[1],b,a},{f[2],c,b},{a,b,c}});
        }
        m.faces=std::move(refined);
    }
    for (auto& p:m.vertices) p=center+Vec3{p.x*radii.x,p.y*radii.y,p.z*radii.z};
    topology(m,true);
    return m;
}
Mesh sphere(Vec3 center,double radius,double max_error) {
    return ellipsoid(center,{radius,radius,radius},max_error);
}
Mesh cylinder(Vec3 center,double radius,double height,double max_error) {
    require(finite(center) && positive(radius) && positive(height) && positive(max_error),
            "invalid cylinder dimensions, center or surface error");
    // Sagitta = 2r sin²(pi/2n), including the interior of the planar side quads.
    double angle=std::asin(std::sqrt(std::min(1.,max_error/radius*0.5)));
    double required=std::max(3.,std::ceil(pi/(2*angle)));
    require(required<=500000,"requested cylinder surface exceeds two million facets");
    int count=int(required);
    Mesh m;
    for (int j=0;j<2;++j) {
        for (int i=0;i<count;++i) {
            double theta=2*pi*i/count;
            m.vertices.push_back(center+Vec3{radius*std::cos(theta),radius*std::sin(theta),(j-0.5)*height});
        }
    }
    int bottom=int(m.vertices.size()),top=bottom+1;
    m.vertices.push_back(center+Vec3{0,0,-height/2});
    m.vertices.push_back(center+Vec3{0,0,height/2});
    for (int i=0;i<count;++i) {
        int j=(i+1)%count;
        m.faces.insert(m.faces.end(),{{i,j,j+count},{i,j+count,i+count},{bottom,j,i},{top,i+count,j+count}});
    }
    topology(m,true);
    return m;
}

namespace {
Mesh read_obj(std::istream& file) {
    Mesh m;
    std::string line;
    while (std::getline(file,line)) {
        line.resize(line.find('#')==std::string::npos ? line.size() : line.find('#'));
        std::istringstream row(line);
        std::string kind;
        row>>kind;
        if (kind=="v") {
            Vec3 p;
            require(bool(row>>p.x>>p.y>>p.z) && finite(p),"invalid OBJ vertex");
            m.vertices.push_back(p);
        } else if (kind=="f") {
            std::vector<int> face;
            std::string word;
            while (row>>word) {
                auto token=word.substr(0,word.find('/'));
                size_t consumed;
                int index=std::stoi(token,&consumed);
                require(consumed==token.size() && index!=0,"invalid OBJ face index");
                index=index<0 ? int(m.vertices.size())+index : index-1;
                require(index>=0 && index<int(m.vertices.size()),"OBJ face index out of range");
                face.push_back(index);
            }
            triangulate(m.vertices,face,m.faces);
        }
    }
    return m;
}
uint32_t little32(const char* p) {
    const auto* u=reinterpret_cast<const unsigned char*>(p);
    return uint32_t(u[0]) | uint32_t(u[1])<<8 | uint32_t(u[2])<<16 | uint32_t(u[3])<<24;
}
double stl_float(const char* p) {
    uint32_t bits=little32(p);
    float value;
    static_assert(sizeof(value)==sizeof(bits),"STL requires 32-bit floats");
    std::memcpy(&value,&bits,sizeof(value));
    return value;
}
Mesh read_stl(std::istream& file) {
    Mesh m;
    std::map<std::tuple<double,double,double>,int> indices;
    auto vertex=[&](Vec3 p) {
        require(finite(p),"nonfinite STL vertex");
        auto entry=indices.try_emplace(std::make_tuple(p.x,p.y,p.z),int(m.vertices.size()));
        if (entry.second) m.vertices.push_back(p);
        return entry.first->second;
    };
    file.seekg(0,std::ios::end);
    auto size=file.tellg();
    file.seekg(0);
    char header[84];
    if (file.read(header,84)) {
        uint32_t count=little32(header+80);
        if (size==std::streamoff(84)+std::streamoff(count)*50) {
            for (uint32_t i=0;i<count;++i) {
                char data[50];
                require(bool(file.read(data,sizeof data)),"truncated binary STL");
                Face f;
                for (int j=0;j<3;++j) {
                    const char* p=data+12+12*j;
                    f[j]=vertex({stl_float(p),stl_float(p+4),stl_float(p+8)});
                }
                m.faces.push_back(f);
            }
            return m;
        }
    }
    file.clear();
    file.seekg(0);
    std::string word;
    auto expect=[&](const char* expected) {
        require(bool(file>>word) && word==expected,"invalid ASCII STL facet");
    };
    expect("solid");
    std::getline(file,word);
    while (file>>word && word=="facet") {
        expect("normal");
        Vec3 n;
        require(bool(file>>n.x>>n.y>>n.z),"invalid STL facet normal");
        expect("outer"); expect("loop");
        Face f;
        for (int& i:f) {
            expect("vertex");
            Vec3 p;
            require(bool(file>>p.x>>p.y>>p.z),"invalid STL coordinates");
            i=vertex(p);
        }
        expect("endloop"); expect("endfacet");
        m.faces.push_back(f);
    }
    require(word=="endsolid","missing STL endsolid");
    return m;
}
Mesh import_surface(const std::string& path) {
    auto ext=path.substr(path.find_last_of('.')+1);
    std::transform(ext.begin(),ext.end(),ext.begin(),[](unsigned char c){ return std::tolower(c); });
    std::ifstream file(path,std::ios::binary);
    require(file.good(),"cannot open mesh file");
    if (ext=="obj") return read_obj(file);
    if (ext=="stl") return read_stl(file);
    throw std::invalid_argument("mesh file must have an OBJ or STL extension");
}

// Surface points count as members of both adjacent zones.
int location(const Mesh& m,Vec3 p,double eps) {
    double solid_angle=0;
    for (auto f:m.faces) {
        auto t=triangle(m,f);
        auto n=unit(cross(t[1]-t[0],t[2]-t[0]));
        if (std::abs(dot(p-t[0],n))<=eps && in_triangle(p,t,n,eps)) return 0;
        auto a=t[0]-p,b=t[1]-p,c=t[2]-p;
        auto na=norm(a),nb=norm(b),nc=norm(c);
        solid_angle+=2*std::atan2(dot(a,cross(b,c)),na*nb*nc+dot(a,b)*nc+dot(b,c)*na+dot(c,a)*nb);
    }
    return std::abs(solid_angle)>2*pi ? 1 : -1;
}

Vec3 closest_on_triangle(Vec3 p,const Triangle& t,Vec3 normal) {
    Vec3 projected=p-normal*dot(p-t[0],normal);
    if (in_triangle(projected,t,normal,0)) return projected;
    double best=std::numeric_limits<double>::infinity();
    Vec3 result{};
    for (int j=0;j<3;++j) {
        Vec3 a=t[j],e=t[(j+1)%3]-a;
        Vec3 q=a+e*std::clamp(dot(p-a,e)/dot(e,e),0.,1.);
        double d=dot(p-q,p-q);
        if (d<best) { best=d; result=q; }
    }
    return result;
}
bool in_interval(double v,double lo,double hi,double scale) {
    double eps=roundoff*std::max({1.,std::abs(v),scale});
    return lo-eps<=v && v<=hi+eps;
}
} // namespace

Mesh import_wall(const std::string& path) {
    Mesh m=import_surface(path);
    validate(m);
    return m;
}
void validate_geometry(const Mesh& m) {
    double eps=tolerance(m);
    for (auto f:m.faces) {
        auto t=triangle(m,f);
        auto a=t[1]-t[0],b=t[2]-t[0];
        require(norm(cross(a,b))>eps*std::max(norm(a),norm(b)),"degenerate triangle");
    }
    require(positive(measure(m).volume),"surface requires positive outward volume");
    embedded(m,eps);
}
int locate_point(const Mesh& m,Vec3 p) {
    require(finite(p),"nonfinite query point");
    return location(m,p,tolerance(m));
}
bool contains_point(const Mesh& m,Vec3 p) {
    return locate_point(m,p)>=0;
}
bool contains(const Mesh& wall,const Mesh& cell) {
    double eps=std::max(tolerance(wall),tolerance(cell));
    if (location(wall,cell.vertices[0],eps)!=1) return false;
    std::vector<Bounds> wall_boxes;
    for (auto f:wall.faces) wall_boxes.push_back(tri_bounds(triangle(wall,f)));
    for (auto f:cell.faces) {
        auto a=triangle(cell,f);
        auto b=tri_bounds(a);
        for (size_t j=0;j<wall.faces.size();++j) {
            if (overlap(b,wall_boxes[j],eps) &&
                intersections(a,triangle(wall,wall.faces[j]),eps,[](Vec3){return true;})) return false;
        }
    }
    return true;
}
Zone layer(int id,Vec3 origin,Vec3 axis,double lower,double upper) {
    require(finite(origin) && lower<=upper,"invalid layer limits or origin");
    Vec3 n=unit(axis);
    return {id,[=](Vec3 p){ auto d=p-origin; return in_interval(dot(d,n),lower,upper,norm(d)); }};
}
Zone radial_band(int id,Vec3 center,double lower,double upper) {
    require(finite(center) && lower>=0 && lower<=upper,"invalid radial band limits or center");
    return {id,[=](Vec3 p){ return in_interval(norm(p-center),lower,upper,0); }};
}
Zone surface_zone(int id,Mesh m,bool inside) {
    validate(m);
    double eps=tolerance(m);
    return {id,[m=std::move(m),inside,eps](Vec3 p) {
        int side=location(m,p,eps);
        return side==0 || (side==1)==inside;
    }};
}
Zone interface_zone(int id,const std::string& path,bool positive_side) {
    Mesh m=import_surface(path);
    topology(m,false);
    double eps=tolerance(m);
    embedded(m,eps);
    return {id,[m=std::move(m),positive_side,eps](Vec3 p) {
        double best=std::numeric_limits<double>::infinity(),signed_distance=0;
        for (auto f:m.faces) {
            auto t=triangle(m,f);
            Vec3 n=unit(cross(t[1]-t[0],t[2]-t[0]));
            Vec3 d=p-closest_on_triangle(p,t,n);
            double distance=norm(d);
            if (distance<best-eps) { best=distance; signed_distance=dot(d,n); }
            else if (std::abs(distance-best)<=eps) signed_distance+=dot(d,n);
        }
        return positive_side ? signed_distance>=-eps : signed_distance<=eps;
    }};
}
int classify(Vec3 centroid,const std::vector<Zone>& ordered) {
    require(finite(centroid),"nonfinite zone centroid");
    for (const auto& zone:ordered) if (zone.matches(centroid)) return zone.id;
    throw std::invalid_argument("no configured zone contains the centroid");
}
} // namespace vertex3d
