#include "../model/three_d/contact.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace vertex3d;

namespace {
void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}
bool near(double actual,double expected,double tol=1e-10) {
    return std::abs(actual-expected)<=tol;
}
Mesh shift(Mesh mesh,Vec3 displacement) {
    for (auto& p:mesh.vertices) p+=displacement;
    return mesh;
}
// Unit-square patches with independent point sets. A shift in x gives a
// rectangular overlap with exactly known area.
Mesh patch(int n,double z,double x_offset,bool reverse) {
    Mesh m;
    for (int j=0;j<=n;++j) for (int i=0;i<=n;++i)
        m.vertices.push_back({x_offset+double(i)/n,double(j)/n,z});
    for (int j=0;j<n;++j) for (int i=0;i<n;++i) {
        int a=(n+1)*j+i,b=a+1,c=a+n+1,d=c+1;
        m.faces.push_back({a,b,d}); m.faces.push_back({a,d,c});
    }
    if (reverse) for (auto& f:m.faces) std::swap(f[0],f[1]);
    return m;
}
Vec3 sum(const std::vector<Contribution>& forces) {
    Vec3 s{}; for (const auto& entry:forces) s+=entry.force; return s;
}
void grid_tests() {
    Mesh a=box({-2,-1,-1},{0,1,1}), b=box({.5,-1,-1},{2,1,1});
    Mesh far=shift(b,{10,0,0});
    MeshList cells{&far,&a,&b};
    std::vector<CellPair> expected{{1,2}};
    for (double width:{0.125,0.25,1.,4.,100.})
        require(contact_pairs(cells,.25,width)==expected,"grid lost pair touching at bin boundary");
    require(expected==exhaustive_pairs(cells,.25),"unexpected exhaustive pairs");
    auto daughters=split(a,{-1,0,0},{1,0,0});
    cells={&daughters.first,&daughters.second,&b,&far};
    for (double width:{.25,1.,3.})
        require(contact_pairs(cells,.3,width)==exhaustive_pairs(cells,.3),"stale post-division bounds");
}
void pressure_tests() {
    constexpr double range=.01,P=4;
    for (double gap:{range,range/2,0.}) {
        for (auto ns:{std::pair{1,1},std::pair{1,4},std::pair{4,4}}) {
            auto a=patch(ns.first,0,0,false), b=patch(ns.second,gap,0,true);
            auto f=contact_forces(a,b,range,P), brute=contact_forces(a,b,range,P,false);
            require(f.contacts==brute.contacts && f.first.size()==brute.first.size(),"face filter changed contacts");
            for (std::size_t k=0;k<f.first.size();++k)
                require(f.first[k].vertex==brute.first[k].vertex &&
                        norm(f.first[k].force-brute.first[k].force)==0,"face filter changed forces");
            Vec3 left=sum(f.first), right=sum(f.second);
            require(near(left.z,-P*(1-gap/range)) && near(norm(left+right),0),"patch pressure or balance");
            require(near(left.x,0) && near(left.y,0),"spurious tangential pressure");
        }
    }
    // A single pair of triangles contributes all three vertex pairs and all
    // three parallel edge-interior pairs. Adding the adjacent triangle cannot
    // duplicate either point or midpoint on their shared diagonal.
    Mesh a=patch(1,0,0,false),b=patch(1,.005,0,true);
    a.faces.resize(1); b.faces.resize(1);
    auto single=contact_forces(a,b,range,P);
    require(single.contacts==6,"triangle pair lost or duplicated a distinct feature contact");
    a=patch(1,0,0,false); b=patch(1,.005,0,true);
    require(contact_forces(a,b,range,P).contacts==9,"shared-edge duplicates counted twice");

    double previous=2;
    for (int n:{2,4,8,16}) {
        auto left=patch(n,0,0,false),right=patch(n,range/2,.5,true);
        double measured=-sum(contact_forces(left,right,range,P).first).z;
        double expected=P*.5*.5;
        double error=std::abs(measured-expected);
        std::cout<<"partial overlap n="<<n<<" force="<<measured<<" expected="<<expected<<" error="<<error<<"\n";
        require(error>0 && error<previous,"partial-overlap corner budget failed to converge");
        previous=error;
    }
    // Crossing middle edges of narrow skew triangles have no nearby vertices.
    Mesh h{{{-1,0,0},{1,0,0},{1,-1,0}},{{0,1,2}},{}};
    Mesh v{{{0,-2,.005},{0,2,.005},{0,2,1}},{{0,1,2}},{}};
    auto force=contact_forces(h,v,range,P);
    require(force.contacts>=1 && norm(sum(force.first))>0,"missed edge-interior approach");
    require(near(norm(sum(force.first)+sum(force.second)),0),"barycentric edge reactions unbalanced");
    // One point projects to the lower triangle at barycentrics (1/2,1/4,1/4).
    // Only that point is within range, so its whole reaction goes to those vertices.
    Mesh base{{{0,0,0},{1,0,0},{0,1,0}},{{0,1,2}},{}};
    Mesh tip{{{.25,.25,.005},{.25,.25,1},{.3,.25,1}},{{0,1,2}},{}};
    auto point_contact=contact_forces(base,tip,range,P);
    require(point_contact.contacts==1,"unique point contact expected");
    double total=sum(point_contact.first).z;
    require(point_contact.first.size()==3 && near(point_contact.first[0].force.z,total*.5) &&
            near(point_contact.first[1].force.z,total*.25) && near(point_contact.first[2].force.z,total*.25),
            "incorrect barycentric force weights");
}
void crossing_tests() {
    Mesh a=box({-2,-.5,-.5},{-1,.5,.5});
    Mesh b=box({1,-.5,-.5},{2,.5,.5});
    auto ae=shift(a,{4,0,0}), be=shift(b,{-4,0,0});
    require(!bodies_intersect(a,b) && !bodies_intersect(ae,be),"exchange endpoints must be disjoint");
    require(swept_crossing(a,ae,b,be),"missed moving cells exchanging sides");
    require(!swept_crossing(a,shift(a,{.5,0,0}),b,shift(b,{-.5,0,0})),"reduced trial rejected");

    Mesh face{{{0,-1,-1},{0,1,-1},{0,0,1}},{{0,1,2}},{}};
    Mesh probe{{{-1,0,0},{-1,-.1,0},{-1,0,.1}},{{0,1,2}},{}};
    require(swept_crossing(face,face,probe,shift(probe,{2,0,0})),"missed mid-trial vertex-face crossing");
    Mesh long_edge{{{-2,0,0},{2,0,0},{2,-.1,0}},{{0,1,2}},{}};
    Mesh other_edge{{{0,-2,-1},{0,2,-1},{.1,2,-1}},{{0,1,2}},{}};
    require(swept_crossing(long_edge,long_edge,other_edge,shift(other_edge,{0,0,2})),
            "missed mid-trial edge-edge crossing without in-range vertices");
    Mesh coplanar=shift(probe,{1,-4,0});
    require(swept_crossing(face,face,coplanar,shift(coplanar,{0,8,0})),"missed coplanar crossing");

    // A moving vertex grazes the interior of a deforming triangle at t=3/8.
    // z_face=t*x and x_point=t, so their signed difference is -(t-3/8)^2.
    // The determinant has no endpoint sign change.
    constexpr double t0=.375;
    Mesh tilt{{{0,0,0},{1,0,0},{0,1,0}},{{0,1,2}},{}}, tilted=tilt;
    tilted.vertices[1].z=1;
    Mesh grazing{{{0,.2,-t0*t0},{0,.21,-t0*t0-.1},{0,.19,-t0*t0-.1}},{{0,1,2}},{}};
    Mesh grazing_end=shift(grazing,{1,0,2*t0});
    require(swept_crossing(tilt,tilted,grazing,grazing_end),"missed repeated-root grazing contact");
    Mesh clear=shift(grazing,{0,0,-.001}),clear_end=shift(grazing_end,{0,0,-.001});
    require(!swept_crossing(tilt,tilted,clear,clear_end),"nearby disjoint deforming trial rejected");

    Mesh cube=box({0,0,0},{1,1,1});
    auto [left,right]=split(cube,{.5,0,0},{1,0,0});
    require(near(measure(left).volume,.5) && near(measure(right).volume,.5),"invalid daughter fixture");
    require(bodies_intersect(left,right),"initial cap contact missing");
    require(!swept_crossing(left,shift(left,{-.1,0,0}),right,shift(right,{.1,0,0})),
            "initially touching daughters cannot separate");
    require(swept_crossing(left,shift(left,{.1,0,0}),right,shift(right,{-.1,0,0})),
            "missed initial cap penetration");
    require(!swept_crossing(left,left,right,right),"resting contact rejected");
    auto reaction=contact_forces(left,right,.01,2);
    require(sum(reaction.first).x<0 && sum(reaction.second).x>0,"cap pressure direction");
    require(near(norm(sum(reaction.first)+sum(reaction.second)),0),"unbalanced daughter reactions");
    auto [oblique_a,oblique_b]=split(cube,{.5,.5,.5},{2,.8,.35});
    Vec3 opening{.01,.004,.00175};
    require(!swept_crossing(oblique_a,shift(oblique_a,-opening),oblique_b,shift(oblique_b,opening)),
            "initial oblique cap opening rejected");
    require(swept_crossing(oblique_a,shift(oblique_a,opening),oblique_b,shift(oblique_b,-opening)),
            "initial oblique cap penetration missed");
    Mesh wall=box({-3,-3,-3},{3,3,3}), through=shift(a,{10,0,0}), safe=shift(a,{1,0,0});
    require(swept_crossing({&a},{&through},wall),"tunneled through stationary wall");
    require(!swept_crossing({&a},{&safe},wall),"valid interior trial rejected");
    Mesh seed=box({-1,-1,-1},{1,1,1}), inner=box({-.2,-.2,-.2},{.2,.2,.2});
    require(bodies_intersect(seed,inner) && !surfaces_intersect(seed,inner),"missed wholly nested bodies");
}
}
int main() {
    grid_tests(); pressure_tests(); crossing_tests();
    std::cout<<"3D grid, contact and continuous crossing fixtures passed\n";
}
