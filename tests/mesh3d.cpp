// Pass a scratch directory for OBJ/STL import fixtures.
#include "mesh.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

using namespace vertex3d;
namespace fs=std::filesystem;
namespace {
void check(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
bool near(double a,double b,double eps=1e-10) { return std::abs(a-b)<eps; }
template<class Operation>
void rejects(Operation op) {
    try { op(); }
    catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("expected invalid geometry to be rejected");
}
void gradient() {
    auto m=box({0,0,0},{1,1,1});
    m.vertices[5]+={0.12,0.23,0.07};
    m.vertices[3]+={-0.01,0.04,-0.1};
    validate_geometry(m);
    for (bool volume:{false,true}) {
        auto energy=[&](const Mesh& a) {
            auto s=measure(a);
            return volume ? 0.5*2.3*std::pow(s.volume-1.4,2) : 0.7*s.area;
        };
        auto f=internal_forces(m,1.4,volume ? 2.3 : 0,volume ? 0 : 0.7);
        for (size_t i=0;i<f.size();++i) for (int k=0;k<3;++k) {
            auto lo=m,hi=m;
            lo.vertices[i][k]-=1e-5; hi.vertices[i][k]+=1e-5;
            check(near(f[i][k],-(energy(hi)-energy(lo))/2e-5,1e-8),"internal energy derivative mismatch");
        }
    }
}
void geometry() {
    auto m=box({0,0,0},{1,1,1});
    auto s=measure(m);
    check(near(s.volume,1) && near(s.area,6) && norm(s.centroid-Vec3{.5,.5,.5})<1e-12,"wrong cube measure");
    check(m.edges.size()==18,"missing edge-face adjacency");
    for (auto [edge,adjacent]:m.edges) for (int side:{0,1}) {
        auto f=m.faces[adjacent[side]];
        bool found=false;
        for (int j=0;j<3;++j)
            if (f[j]==(side ? edge.second : edge.first) && f[(j+1)%3]==(side ? edge.first : edge.second)) found=true;
        check(found,"wrong traversal in face adjacency");
    }
    auto cuts=split(m,s.centroid,{1,0,0});
    for (const Mesh* side:{&cuts.first,&cuts.second}) {
        check(near(measure(*side).volume,.5) && near(measure(*side).area,4),"bisected cube measure");
        validate_geometry(*side);
        double cap_area=0;
        for (auto f:side->faces) {
            auto a=side->vertices[f[0]],b=side->vertices[f[1]],c=side->vertices[f[2]];
            if (a.x==.5 && b.x==.5 && c.x==.5) cap_area+=cross(b-a,c-a).x/2;
        }
        check(near(cap_area,side==&cuts.first ? 1 : -1),"wrong signed daughter cap area");
    }
    for (Vec3 n:std::vector<Vec3>{{1,1,0},{0,1,1},{1,2,3},{1,1,1},{-.4,.7,.2}}) {
        auto daughters=split(m,s.centroid,n);
        check(near(measure(daughters.first).volume,.5) && near(measure(daughters.second).volume,.5),"oblique cut");
    }
    rejects([&]{ split(m,{0,0,0},{1,0,0}); });
    rejects([&]{ split(m,s.centroid,{}); });
    check(m.faces.size()==12 && m.vertices.size()==8 && measure(m).volume==s.volume,"rejected cut modified parent");

    auto bad=m;
    bad.faces.pop_back(); rejects([&]{validate(bad);});
    bad=m; std::swap(bad.faces[0][0],bad.faces[0][1]); rejects([&]{validate(bad);});
    bad=m; bad.vertices[0]=bad.vertices[1]; rejects([&]{validate(bad);});
    bad=m; bad.vertices[6]={.4,.4,-.3}; rejects([&]{validate(bad);});
    bad=m; bad.vertices[0].x=std::numeric_limits<double>::infinity(); rejects([&]{validate(bad);});
    bad=m; for (auto& f:bad.faces) std::swap(f[0],f[1]); rejects([&]{validate(bad);});
}
void zones(const fs::path& tmp) {
    auto m=box({-1,-1,-1},{1,1,1});
    auto left=layer(7,{},{1,0,0},-100,0),right=layer(9,{},{1,0,0},0,100);
    for (auto zones:std::vector<std::vector<Zone>>{{left,right},{right,left}}) {
        check(classify(measure(m).centroid,zones)==zones[0].id,"ordered interface tie");
        for (double dx:{-.25,.25}) {
            auto moved=m;
            for (auto& p:moved.vertices) p.x+=dx;
            check(classify(measure(moved).centroid,zones)==(dx<0 ? 7 : 9),"zone off interface");
        }
        auto refined=m;
        auto f=m.faces[6];
        auto p=(m.vertices[f[0]]+m.vertices[f[1]]+m.vertices[f[2]])/3;
        refined.vertices.push_back(p);
        refined.faces[6]={f[0],f[1],8};
        refined.faces.push_back({f[1],f[2],8});
        refined.faces.push_back({f[2],f[0],8});
        validate(refined);
        check(norm(measure(refined).centroid-measure(m).centroid)<1e-12,"refined centroid");
        check(classify(measure(refined).centroid,zones)==zones[0].id,"refined mesh changed zone");
    }
    auto interface=tmp/"interface.obj";
    std::ofstream(interface)<<"v 0 -2 -2\nv 0 2 -2\nv 0 2 2\nv 0 -2 2\nf 1 2 3 4\n";
    auto positive=interface_zone(9,interface.string(),true),negative=interface_zone(7,interface.string(),false);
    check(classify({.5,0,0},{negative,positive})==9 && classify({-.5,0,0},{positive,negative})==7 &&
          classify({0,1,0},{negative,positive})==7,"oriented imported interface");
    auto inner=surface_zone(1,m,true),outer=surface_zone(2,m,false);
    check(classify({1,0,0},{outer,inner})==2 && classify({0,0,0},{outer,inner})==1,"closed imported zones");
    check(radial_band(5,{2,1,0},1,2).matches({2,2,0}),"radial band");
    double inf=std::numeric_limits<double>::infinity();
    auto example=box({-.5,-.5,-.5},{.5,.5,.5});
    check(classify(measure(example).centroid,
        {layer(10,{},{1,0,0},-inf,0),layer(20,{},{1,0,0},0,inf)})==10,"example cube zone at reset");
}
void translated_cuts() {
    std::mt19937 rng(94);
    std::uniform_real_distribution<double> component(-1,1);
    for (double scale:{1e-3,1.,1e3}) for (double shift:{0.,10000.}) {
        Vec3 center{shift,0,0},half{scale/2,scale/2,scale/2};
        auto cube=box(center-half,center+half);
        for (int i=0;i<30;++i) {
            Vec3 n{component(rng),component(rng),component(rng)};
            if (i==0) n={1,1,0};
            if (i==1) n={1,1,1};
            auto halves=split(cube,center,n);
            for (const auto* part:{&halves.first,&halves.second})
                check(near(measure(*part).volume/measure(cube).volume,.5,1e-6),"translated half volume");
        }
    }
}
void binary_word(std::ostream& out,uint32_t value) {
    for (int i=0;i<4;++i) out.put(char(value>>(8*i)));
}
void imports(const fs::path& tmp) {
    auto m=box({0,0,0},{1,1,1});
    auto obj=tmp/"wall.obj",ascii=tmp/"ascii.stl",binary=tmp/"binary.stl";
    std::ofstream o(obj),a(ascii),b(binary,std::ios::binary);
    for (auto p:m.vertices) o<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    a<<"solid wall\n";
    b<<std::string(80,' '); binary_word(b,m.faces.size());
    for (auto f:m.faces) {
        o<<"f "<<f[0]+1<<" "<<f[1]+1<<" "<<f[2]+1<<"\n";
        a<<"facet normal 0 0 0\nouter loop\n"; b<<std::string(12,0);
        for (int k:f) {
            auto p=m.vertices[k]; a<<"vertex "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
            for (int j=0;j<3;++j) {
                float value=float(p[j]); uint32_t bits; std::memcpy(&bits,&value,4); binary_word(b,bits);
            }
        }
        a<<"endloop\nendfacet\n"; b<<std::string(2,0);
    }
    a<<"endsolid wall\n"; a.close(); b.close(); o.close();
    for (auto file:{obj,ascii,binary}) {
        auto wall=import_wall(file.string());
        check(near(measure(wall).volume,1) && near(measure(wall).area,6),"wrong imported wall");
        check(contains(wall,box({.2,.2,.2},{.8,.8,.8})) && !contains(wall,box({.5,.5,.5},{1.2,1.2,1.2})),
              "imported wall containment");
        check(contains_point(wall,{1,.5,.5}) && !contains_point(wall,{1.1,.5,.5}),"boundary point location");
    }
    auto missing=tmp/"open.obj";
    std::ofstream(missing)<<"v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
    rejects([&]{import_wall(missing.string());});

    // All eight seed corners are interior to this U-shaped wall. Its triangles
    // bridge the two arms and cross the notch, so vertex-only containment fails.
    auto concave=tmp/"concave.obj";
    std::vector<Vec3> rim{{0,0,0},{3,0,0},{3,3,0},{2,3,0},{2,1,0},{1,1,0},{1,3,0},{0,3,0}};
    std::ofstream u(concave);
    for (int z:{0,1}) for (auto p:rim) u<<"v "<<p.x<<" "<<p.y<<" "<<z<<"\n";
    u<<"f 8 7 6 5 4 3 2 1\nf 9 10 11 12 13 14 15 16\n";
    for (int i=0;i<8;++i) {
        int j=(i+1)%8;
        u<<"f "<<i+1<<" "<<j+1<<" "<<j+9<<" "<<i+9<<"\n";
    }
    u.close();
    auto wall=import_wall(concave.string()),bridge=box({.5,1.5,.25},{2.5,2.5,.75});
    for (auto p:bridge.vertices) check(locate_point(wall,p)==1,"expected bridge corner inside wall");
    check(!contains(wall,bridge),"bridging cell spans a concave wall");
    auto parts=split(wall,{0,0,.5},{0,0,1});
    check(near(measure(parts.first).volume,3.5) && near(measure(parts.second).volume,3.5),
          "nonconvex cap triangulation");

    Mesh tube;
    for (double z:{0.,1.}) for (double r:{2.,1.})
        for (Vec3 p:std::vector<Vec3>{{-r,-r,z},{r,-r,z},{r,r,z},{-r,r,z}}) tube.vertices.push_back(p);
    auto quad=[&](int a,int b,int c,int d) {
        tube.faces.push_back({a,b,c}); tube.faces.push_back({a,c,d});
    };
    for (int i=0;i<4;++i) {
        int j=(i+1)%4;
        quad(i,j,j+8,i+8); quad(j+4,i+4,i+12,j+12);
        quad(i+8,j+8,j+12,i+12); quad(i+4,j+4,j,i);
    }
    validate(tube);
    rejects([&]{split(tube,{0,0,.5},{0,0,1});});
    check(near(measure(tube).volume,12) && tube.faces.size()==32,"rejected annular cut modified parent");
}
void curves() {
    double error=.05;
    for (Vec3 radii:std::vector<Vec3>{{10,10,10},{5,10,2}}) {
        auto m=ellipsoid({},radii,error);
        for (auto f:m.faces) {
            Vec3 q[3];
            for (int i=0;i<3;++i) {
                auto p=m.vertices[f[i]];
                q[i]={p.x/radii.x,p.y/radii.y,p.z/radii.z};
            }
            auto n=cross(q[1]-q[0],q[2]-q[0]);
            double h=dot(n,q[0])/norm(n);
            check(10*(1-h)<=error,"facet interior exceeds global deviation bound");
        }
    }
    auto c=cylinder({},10,4,error);
    double deepest=0;
    for (auto f:c.faces) {
        auto a=c.vertices[f[0]],b=c.vertices[f[1]],d=c.vertices[f[2]];
        if (a.z==b.z && b.z==d.z) continue;
        Vec3 chord=a.z==b.z ? (a+b)/2 : (b+d)/2;
        deepest=std::max(deepest,10-std::hypot(chord.x,chord.y));
    }
    check(deepest>0 && deepest<=error,"cylinder facet sagitta");
    rejects([]{sphere({},1,0);}); rejects([]{cylinder({},1,1,-.1);});
}
}
int main(int argc,char** argv) {
    try {
        check(argc==2,"pass a scratch directory");
        fs::create_directories(argv[1]);
        gradient(); geometry(); zones(argv[1]); imports(argv[1]); curves(); translated_cuts();
        std::cout<<"Mesh mechanics, division, zones, curved bounds and imports passed\n";
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<"\n"; return 1;
    }
}
