#include "contact.h"

#include <algorithm>
#include <limits>
#include <set>
#include <stdexcept>

namespace vertex3d {
namespace {
constexpr double eps = 64 * std::numeric_limits<double>::epsilon();
using Triangle = std::array<Vec3,3>;
using Weights = std::array<double,3>;

struct Bounds {
    Vec3 lo, hi;
    explicit Bounds(Vec3 p) : lo(p), hi(p) {}
    void include(Vec3 p) {
        for (int d=0; d<3; ++d) {
            lo[d]=std::min(lo[d],p[d]);
            hi[d]=std::max(hi[d],p[d]);
        }
    }
    void expand(double r) {
        Vec3 pad{r,r,r}; lo-=pad; hi+=pad;
    }
    bool overlaps(const Bounds& b) const {
        for (int d=0; d<3; ++d)
            if (hi[d]<b.lo[d] || b.hi[d]<lo[d]) return false;
        return true;
    }
};

Bounds bounds(const Mesh& m) {
    Bounds b(m.vertices.at(0));
    for (auto p:m.vertices) b.include(p);
    return b;
}
Triangle triangle(const Mesh& m,const Face& f) {
    return {m.vertices[f[0]],m.vertices[f[1]],m.vertices[f[2]]};
}
Bounds bounds(const Triangle& t) {
    Bounds b(t[0]); b.include(t[1]); b.include(t[2]); return b;
}
std::vector<Bounds> expanded_boxes(const MeshList& cells,double range) {
    if (!(std::isfinite(range) && range>=0)) throw std::invalid_argument("invalid contact range");
    std::vector<Bounds> boxes;
    for (const Mesh* cell:cells) {
        boxes.push_back(bounds(*cell)); boxes.back().expand(range);
    }
    return boxes;
}
std::vector<CellPair> scan_boxes(const std::vector<Bounds>& boxes) {
    std::vector<CellPair> pairs;
    for (int i=0; i<int(boxes.size()); ++i)
        for (int j=i+1; j<int(boxes.size()); ++j)
            if (boxes[i].overlaps(boxes[j])) pairs.emplace_back(i,j);
    return pairs;
}
Vec3 weighted(const Triangle& t,Weights w) {
    return t[0]+(t[1]-t[0])*w[1]+(t[2]-t[0])*w[2];
}
Weights vertex(int i) { Weights w{}; w[i]=1; return w; }
Weights edge(int i,double u) { Weights w{}; w[i]=1-u; w[(i+1)%3]=u; return w; }
Vec3 normal(const Triangle& t) { return cross(t[1]-t[0],t[2]-t[0]); }
double span(const Triangle& t) {
    return std::max({norm(t[1]-t[0]),norm(t[2]-t[0]),norm(t[2]-t[1])});
}

// Orthogonal projection if it falls in the face, otherwise the closest edge point.
Weights closest(Vec3 p,const Triangle& t) {
    Vec3 u=t[1]-t[0], v=t[2]-t[0], q=p-t[0], n=cross(u,v);
    double n2=dot(n,n);
    if (n2>0) {
        double s=dot(cross(q,v),n)/n2, r=dot(cross(u,q),n)/n2;
        if (s>=0 && r>=0 && s+r<=1) return {1-s-r,s,r};
    }
    double best=std::numeric_limits<double>::infinity();
    Weights weights=vertex(0);
    for (int i=0;i<3;++i) {
        Vec3 d=t[(i+1)%3]-t[i];
        double l2=dot(d,d);
        double k=l2 ? std::clamp(dot(p-t[i],d)/l2,0.,1.) : 0.;
        Vec3 delta=p-(t[i]+k*d);
        double distance2=dot(delta,delta);
        if (distance2<best) { best=distance2; weights=edge(i,k); }
    }
    return weights;
}

// For parallel segments the overlap of their projections has a continuum of
// closest pairs. Its midpoint represents the edge-interior interaction.
bool edge_interiors(Vec3 p,Vec3 q,Vec3 r,Vec3 s,double& u,double& v) {
    Vec3 a=q-p,b=s-r,h=r-p,n=cross(a,b);
    double aa=dot(a,a), bb=dot(b,b), nn=dot(n,n);
    if (!aa || !bb) return false;
    if (nn>eps*eps*aa*bb) {
        u=dot(cross(h,b),n)/nn;
        v=dot(cross(h,a),n)/nn;
    } else {
        double from=dot(h,a)/aa,to=dot(s-p,a)/aa;
        double lo=std::max(0.,std::min(from,to)), hi=std::min(1.,std::max(from,to));
        if (lo>=hi) return false;
        u=(lo+hi)/2;
        v=dot(p+u*a-r,b)/bb;
    }
    return u>eps && u<1-eps && v>eps && v<1-eps;
}
struct Location {
    // Positive supports sorted by vertex ID. The count distinguishes face,
    // edge and vertex. A shared edge therefore has the same key on both faces.
    Face ids{-1,-1,-1};
    Weights weights{};
};
Location location(const Face& f,Weights w) {
    std::vector<std::pair<int,double>> active;
    double total=0;
    for (int i=0;i<3;++i) if (w[i]>eps) {
        active.emplace_back(f[i],w[i]); total+=w[i];
    }
    std::sort(active.begin(),active.end());
    Location loc;
    for (std::size_t i=0;i<active.size();++i) {
        loc.ids[i]=active[i].first; loc.weights[i]=active[i].second/total;
    }
    return loc;
}
struct Contact {
    Location a,b;
    Vec3 delta{},zero_direction{};
    double normal_alignment=2;
    std::set<int> supports_a,supports_b;
};
using Contacts=std::map<std::pair<Face,Face>,Contact>;
struct Facet {
    Face ids;
    Triangle t;
    Bounds box;
    Vec3 unit;
    double area;
};
std::vector<Facet> facets(const Mesh& m,double range) {
    std::vector<Facet> fs;
    for (const auto& face:m.faces) {
        Triangle t=triangle(m,face); Vec3 n=normal(t);
        double length=norm(n);
        if (!(length>0)) throw std::invalid_argument("degenerate contact face");
        fs.push_back({face,t,bounds(t),n/length,length/2});
        fs.back().box.expand(range);
    }
    return fs;
}
void collect(Contacts& contacts,const Facet& a,const Facet& b,int ia,int ib,double range) {
    auto add=[&](Weights wa,Weights wb) {
        Vec3 delta=weighted(a.t,wa)-weighted(b.t,wb);
        if (norm(delta)>=range) return;
        Location la=location(a.ids,wa), lb=location(b.ids,wb);
        // A vertex projection has a unique point on its opposite feature. Two
        // edge interiors similarly have one closest pair (or the midpoint above).
        auto [it,inserted]=contacts.try_emplace({la.ids,lb.ids});
        Contact& c=it->second;
        if (inserted) { c.a=la; c.b=lb; c.delta=delta; }
        c.supports_a.insert(ia); c.supports_b.insert(ib);
        double alignment=dot(a.unit,b.unit);
        if (alignment<c.normal_alignment) {
            c.normal_alignment=alignment;
            Vec3 reaction=b.unit-a.unit;
            c.zero_direction=norm(reaction)>eps ? reaction/norm(reaction) : -a.unit;
        }
    };
    for (int i=0;i<3;++i) {
        add(vertex(i),closest(a.t[i],b.t));
        add(closest(b.t[i],a.t),vertex(i));
    }
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
        double u,v;
        if (edge_interiors(a.t[i],a.t[(i+1)%3],b.t[j],b.t[(j+1)%3],u,v))
            add(edge(i,u),edge(j,v));
    }
}
void distribute(std::vector<Contribution>& result,const Location& loc,Vec3 force) {
    for (int k=0;k<3 && loc.ids[k]>=0;++k)
        result.push_back({loc.ids[k],loc.weights[k]*force});
}
}

std::vector<CellPair> exhaustive_pairs(const MeshList& cells,double range) {
    return scan_boxes(expanded_boxes(cells,range));
}
std::vector<CellPair> contact_pairs(const MeshList& cells,double range,double bin_width) {
    if (!(std::isfinite(bin_width) && bin_width>0)) throw std::invalid_argument("invalid bin width");
    auto boxes=expanded_boxes(cells,range);
    using Bin=std::array<long long,3>;
    std::map<Bin,std::vector<int>> bins;
    std::set<CellPair> pairs;
    for (int i=0;i<int(boxes.size());++i) {
        Bin lo{},hi{};
        long double bin_count=1;
        for (int d=0;d<3;++d) {
            auto bin=[&](double x){return std::floor((long double)x/bin_width);};
            long double a=bin(boxes[i].lo[d]),b=bin(boxes[i].hi[d]);
            // With an extreme width, checking all real pairs is cheaper than
            // allocating millions of empty spatial bins or overflowing indices.
            if (a<=std::numeric_limits<long long>::min() || b>=std::numeric_limits<long long>::max())
                return scan_boxes(boxes);
            bin_count*=b-a+1;
            lo[d]=static_cast<long long>(a); hi[d]=static_cast<long long>(b);
        }
        if (bin_count>1000000) return scan_boxes(boxes);
        for (auto x=lo[0];x<=hi[0];++x) for (auto y=lo[1];y<=hi[1];++y)
            for (auto z=lo[2];z<=hi[2];++z) {
                auto& occupants=bins[{x,y,z}];
                for (int j:occupants) if (boxes[j].overlaps(boxes[i])) pairs.emplace(j,i);
                occupants.push_back(i);
            }
    }
    return {pairs.begin(),pairs.end()};
}
double contact_law(double gap,double range,bool barrier) {
    if (!(gap<range)) return 0;
    if (!barrier) return 1-gap/range;
    return range/std::max(gap,1e-4*range)-1;
}
PairForces contact_forces(const Mesh& a,const Mesh& b,double range,double pressure,bool box_filter,bool barrier) {
    if (!(std::isfinite(range) && range>0 && std::isfinite(pressure) && pressure>=0))
        throw std::invalid_argument("invalid contact parameters");
    auto fa=facets(a,range), fb=facets(b,range);
    Contacts contacts;
    for (int i=0;i<int(fa.size());++i) for (int j=0;j<int(fb.size());++j)
        if (!box_filter || fa[i].box.overlaps(fb[j].box)) collect(contacts,fa[i],fb[j],i,j,range);
    std::vector<int> counts_a(fa.size()),counts_b(fb.size());
    for (const auto& item:contacts) {
        for (int i:item.second.supports_a) ++counts_a[i];
        for (int i:item.second.supports_b) ++counts_b[i];
    }
    PairForces result; result.contacts=contacts.size();
    for (const auto& item:contacts) {
        const Contact& c=item.second;
        double area=0;
        for (int i:c.supports_a) area+=fa[i].area/counts_a[i];
        for (int i:c.supports_b) area+=fb[i].area/counts_b[i];
        double gap=norm(c.delta);
        Vec3 direction=gap>eps*std::min(span(fa[*c.supports_a.begin()].t),span(fb[*c.supports_b.begin()].t))
                     ? c.delta/gap : c.zero_direction;
        Vec3 force=direction*(0.5*area*pressure*contact_law(gap,range,barrier));
        distribute(result.first,c.a,force);
        distribute(result.second,c.b,-force);
    }
    return result;
}

namespace {
bool on_triangle(Vec3 p,const Triangle& t,double tol) {
    return norm(p-weighted(t,closest(p,t)))<=tol;
}
bool segments_touch(Vec3 p,Vec3 q,Vec3 r,Vec3 s,double tol) {
    auto near=[&](Vec3 x,Vec3 a,Vec3 b) {
        Vec3 d=b-a; double dd=dot(d,d);
        double u=dd ? std::clamp(dot(x-a,d)/dd,0.,1.) : 0.;
        return norm((x-a)-u*d)<=tol;
    };
    if (near(p,r,s)||near(q,r,s)||near(r,p,q)||near(s,p,q)) return true;
    double u,v;
    return edge_interiors(p,q,r,s,u,v) && norm((p-r)+u*(q-p)-v*(s-r))<=tol;
}
bool edge_cuts_face(Vec3 p,Vec3 q,const Triangle& t,double tol) {
    Vec3 n=normal(t); double h=dot(p-t[0],n), k=dot(q-t[0],n);
    if ((h<0 && k>0)||(h>0 && k<0))
        return on_triangle(p+(q-p)*(h/(h-k)),t,tol);
    return false;
}
bool triangles_touch(const Triangle& a,const Triangle& b) {
    if (!bounds(a).overlaps(bounds(b))) return false;
    double tol=eps*std::max(span(a),span(b));
    for (int i=0;i<3;++i) {
        int next=(i+1)%3;
        if (on_triangle(a[i],b,tol)||on_triangle(b[i],a,tol)) return true;
        if (edge_cuts_face(a[i],a[next],b,tol)||edge_cuts_face(b[i],b[next],a,tol)) return true;
        for (int j=0;j<3;++j)
            if (segments_touch(a[i],a[next],b[j],b[(j+1)%3],tol)) return true;
    }
    return false;
}

// Isolate all real roots in [0,1]. Between consecutive derivative roots a
// cubic is monotone. Checking the extrema also retains repeated (grazing) roots.
using Polynomial=std::array<long double,4>;
long double evaluate(const Polynomial& c,int degree,long double x) {
    long double y=c[degree];
    while (degree) y=y*x+c[--degree];
    return y;
}
void roots(const Polynomial& c,int degree,std::vector<double>& result) {
    long double scale=0;
    for (int i=0;i<=degree;++i) scale+=std::abs(c[i]);
    if (!scale) return;
    while (degree && std::abs(c[degree])<=eps*scale) --degree;
    if (!degree) return;
    if (degree==1) {
        long double t=-c[0]/c[1];
        if (t>=0 && t<=1) result.push_back(double(t));
        return;
    }
    Polynomial derivative{};
    for (int i=1;i<=degree;++i) derivative[i-1]=i*c[i];
    std::vector<double> knots{0,1};
    roots(derivative,degree-1,knots);
    std::sort(knots.begin(),knots.end());
    for (double x:knots) if (std::abs(evaluate(c,degree,x))<=eps*scale) result.push_back(x);
    for (std::size_t i=1;i<knots.size();++i) {
        long double lo=knots[i-1],hi=knots[i];
        long double f=evaluate(c,degree,lo), g=evaluate(c,degree,hi);
        if ((f>=0 && g>=0)||(f<=0 && g<=0)) continue;
        for (int k=0;k<64;++k) {
            long double mid=(lo+hi)/2, m=evaluate(c,degree,mid);
            if ((f<0 && m<0)||(f>0 && m>0)) lo=mid;
            else hi=mid;
        }
        result.push_back(double((lo+hi)/2));
    }
}
struct Path {
    Vec3 p,v;
    Vec3 at(double t) const { return p+v*t; }
};
Path difference(Path a,Path b) { return {a.p-b.p,a.v-b.v}; }
using Four=std::array<Path,4>;
using MovingTriangle=std::array<Path,3>;
MovingTriangle paths(const Mesh& m,const Mesh& end,const Face& f) {
    MovingTriangle t;
    for (int i=0;i<3;++i) {
        int j=f[i];
        t[i]={m.vertices[j],end.vertices[j]-m.vertices[j]};
    }
    return t;
}
Bounds swept_bounds(const MovingTriangle& t) {
    Bounds b(t[0].p);
    for (const auto& p:t) { b.include(p.p); b.include(p.at(1)); }
    return b;
}
Bounds swept_bounds(const Mesh& m,const Mesh& end) {
    Bounds b=bounds(m); for (auto p:end.vertices) b.include(p); return b;
}
struct SweptFace {
    MovingTriangle path;
    Bounds bounds;
};
std::vector<SweptFace> swept_faces(const Mesh& a,const Mesh& end) {
    std::vector<SweptFace> faces;
    for (const Face& f:a.faces) {
        auto t=paths(a,end,f); faces.push_back({t,swept_bounds(t)});
    }
    return faces;
}
void boundary_roots(Path u,Path v,std::vector<double>& times) {
    // Collinearity of a moving point and edge is quadratic in each projection.
    Vec3 c0=cross(u.p,v.p),c1=cross(u.v,v.p)+cross(u.p,v.v),c2=cross(u.v,v.v);
    for (int axis=0;axis<3;++axis) roots({c0[axis],c1[axis],c2[axis],0},2,times);
}
enum class FeatureTest {vertex_face,edge_edge};

struct TrialEvents {
    double first_positive=1;
    bool initial_touch=false;
};
bool crossing_feature(const Four& p,FeatureTest test,double tol,TrialEvents& events) {
    Path u=difference(p[1],p[0]),v=difference(p[2],p[0]),w=difference(p[3],p[0]);
    Vec3 n0=cross(u.p,v.p), n1=cross(u.p,v.v)+cross(u.v,v.p), n2=cross(u.v,v.v);
    Polynomial c{dot(n0,w.p), dot(n0,w.v)+dot(n1,w.p),
                 dot(n1,w.v)+dot(n2,w.p),dot(n2,w.v)};
    std::vector<double> times{0,1};
    roots(c,3,times);
    double magnitude=(norm(u.p)+norm(u.v))*(norm(v.p)+norm(v.v))*(norm(w.p)+norm(w.v));
    bool coplanar=true;
    for (auto coefficient:c) if (std::abs(coefficient)>eps*magnitude) coplanar=false;
    if (coplanar) {
        // When the coplanarity determinant vanishes for an interval, the
        // containment predicate can only change at a collinearity event.
        // Endpoint coincidences cover trajectories confined to a single line.
        for (int i=0;i<4;++i) for (int j=i+1;j<4;++j) {
            Path d=difference(p[j],p[i]);
            for (int k=j+1;k<4;++k) boundary_roots(d,difference(p[k],p[i]),times);
            for (int axis=0;axis<3;++axis) roots({d.p[axis],d.v[axis],0,0},1,times);
        }
    }
    std::sort(times.begin(),times.end());
    times.erase(std::unique(times.begin(),times.end()),times.end());
    // Evaluate in the frame of the first moving vertex. A common translation
    // cannot change a feature contact and should not cost coordinate precision.
    auto touches=[&](double t) {
        if (test==FeatureTest::vertex_face)
            return on_triangle(w.at(t),{Vec3{},u.at(t),v.at(t)},tol);
        return segments_touch({},u.at(t),v.at(t),w.at(t),tol);
    };
    bool initial=touches(0);
    events.initial_touch|=initial;
    bool free=!initial;
    for (std::size_t i=1;i<times.size();++i) {
        double t=times[i];
        events.first_positive=std::min(events.first_positive,t);
        // Mid-interval tests distinguish a persistent initial contact from a
        // separation followed by a return, including wholly coplanar motion.
        for (double probe:{(times[i-1]+t)/2,t}) {
            if (touches(probe)) { if (free) return true; }
            else free=true;
        }
    }
    return false;
}
bool initial_patch_enters(const MovingTriangle& ma,const MovingTriangle& mb,bool enclosing) {
    Triangle a{ma[0].p,ma[1].p,ma[2].p},b{mb[0].p,mb[1].p,mb[2].p};
    if (!bounds(a).overlaps(bounds(b))) return false;
    Vec3 na=normal(a),nb=normal(b);
    double al=norm(na),bl=norm(nb);
    if (!al || !bl) return false;
    na=na/al; nb=nb/bl;
    if (norm(na-(enclosing ? nb : -nb))>8*eps) return false;
    double length=std::max(span(a),span(b));
    if (std::abs(dot(a[0]-b[0],nb))>eps*length) return false;
    // Clip the coincident patch to both triangles. At each point of positive
    // area the two surface velocities are affine in the barycentric weights,
    // so the largest inward normal speed occurs at a clipped polygon vertex.
    std::vector<Vec3> polygon(a.begin(),a.end());
    for (int i=0;i<3 && !polygon.empty();++i) {
        Vec3 edge_vector=b[(i+1)%3]-b[i];
        Vec3 in=cross(nb,edge_vector);
        std::vector<Vec3> clipped;
        Vec3 p=polygon.back();
        double h=dot(p-b[i],in);
        for (Vec3 q:polygon) {
            double k=dot(q-b[i],in);
            if ((h<0 && k>=0)||(k<0 && h>=0)) clipped.push_back(p+(q-p)*(h/(h-k)));
            if (k>=0) clipped.push_back(q);
            p=q; h=k;
        }
        polygon=std::move(clipped);
    }
    if (polygon.size()<3) return false;
    Vec3 twice_area{};
    for (std::size_t i=2;i<polygon.size();++i)
        twice_area+=cross(polygon[i-1]-polygon[0],polygon[i]-polygon[0]);
    if (norm(twice_area)<=eps*std::min(al,bl)) return false;
    Triangle va{ma[0].v,ma[1].v,ma[2].v},vb{mb[0].v,mb[1].v,mb[2].v};
    for (Vec3 p:polygon) {
        Vec3 relative=weighted(va,closest(p,a))-weighted(vb,closest(p,b));
        if (dot(na,relative)>eps*std::max(length,norm(relative))) return true;
    }
    return false;
}
bool crossing_faces(const MovingTriangle& a,const MovingTriangle& b,TrialEvents& events) {
    Triangle ta{a[0].p,a[1].p,a[2].p},tb{b[0].p,b[1].p,b[2].p};
    double tol=eps*std::max(span(ta),span(tb));
    for (int i=0;i<3;++i) {
        if (crossing_feature({a[0],a[1],a[2],b[i]},FeatureTest::vertex_face,tol,events) ||
            crossing_feature({b[0],b[1],b[2],a[i]},FeatureTest::vertex_face,tol,events))
            return true;
        for (int j=0;j<3;++j)
            if (crossing_feature({a[i],a[(i+1)%3],b[j],b[(j+1)%3]},FeatureTest::edge_edge,tol,events))
                return true;
    }
    return false;
}
Mesh interpolated(const Mesh& a,const Mesh& b,double t) {
    Mesh m=a;
    for (std::size_t i=0;i<m.vertices.size();++i) m.vertices[i]+=(b.vertices[i]-a.vertices[i])*t;
    return m;
}
bool sample_intrusion(const Mesh& a,const Mesh& b,bool enclosing) {
    auto forbidden=[&](Vec3 p) {
        int side=locate_point(b,p);
        return enclosing ? side<0 : side>0;
    };
    for (Vec3 p:a.vertices) if (forbidden(p)) return true;
    // After coincident daughter caps enter each other their rim vertices may
    // still lie on side faces, while the interiors of the caps are embedded.
    for (const Face& f:a.faces) {
        Triangle t=triangle(a,f);
        if (forbidden(t[0]+((t[1]-t[0])+(t[2]-t[0]))/3)) return true;
    }
    return false;
}
bool transverse_intersection(const Mesh& a,const Mesh& b) {
    // An edge piercing a face strictly inside both features also detects a
    // partial cap intrusion whose vertices and triangle centroids remain outside.
    for (const Face& af:a.faces) {
        Triangle ta=triangle(a,af); Bounds ab=bounds(ta);
        for (const Face& bf:b.faces) {
            Triangle tb=triangle(b,bf);
            if (!ab.overlaps(bounds(tb))) continue;
            Vec3 n=normal(tb);
            double tol=eps*norm(n)*std::max(span(ta),span(tb));
            for (int i=0;i<3;++i) {
                Vec3 p=ta[i],q=ta[(i+1)%3];
                double h=dot(p-tb[0],n), k=dot(q-tb[0],n);
                if (!((h>tol && k<-tol)||(h<-tol && k>tol))) continue;
                Weights w=closest(p+(q-p)*(h/(h-k)),tb);
                if (*std::min_element(w.begin(),w.end())>eps) return true;
            }
        }
    }
    return false;
}
bool crosses(const Mesh& a,const Mesh& ae,const Mesh& b,const Mesh& be,bool enclosing) {
    if (a.vertices.size()!=ae.vertices.size() || b.vertices.size()!=be.vertices.size() ||
        a.faces!=ae.faces || b.faces!=be.faces) throw std::invalid_argument("trial mesh topology changed");
    if (!swept_bounds(a,ae).overlaps(swept_bounds(b,be))) return false;
    auto fa=swept_faces(a,ae),fb=swept_faces(b,be);
    TrialEvents events;
    for (const auto& x:fa) for (const auto& y:fb)
        if (x.bounds.overlaps(y.bounds) &&
            (initial_patch_enters(x.path,y.path,enclosing) || crossing_faces(x.path,y.path,events))) return true;
    if (events.initial_touch) {
        // Determine on which side an initial contact departs before the first
        // possible subsequent event. Coplanarity alone has no sign at birth.
        auto trial_a=interpolated(a,ae,events.first_positive/2);
        auto trial_b=interpolated(b,be,events.first_positive/2);
        if (transverse_intersection(trial_a,trial_b) || transverse_intersection(trial_b,trial_a)) return true;
        if (sample_intrusion(trial_a,trial_b,enclosing)) return true;
        if (!enclosing && sample_intrusion(trial_b,trial_a,false)) return true;
    }
    return false;
}
}

bool surfaces_intersect(const Mesh& a,const Mesh& b) {
    if (!bounds(a).overlaps(bounds(b))) return false;
    for (const auto& x:a.faces) for (const auto& y:b.faces)
        if (triangles_touch(triangle(a,x),triangle(b,y))) return true;
    return false;
}
bool bodies_intersect(const Mesh& a,const Mesh& b) {
    if (!bounds(a).overlaps(bounds(b))) return false;
    return surfaces_intersect(a,b) || contains_point(a,b.vertices.at(0)) || contains_point(b,a.vertices.at(0));
}
bool swept_crossing(const Mesh& a,const Mesh& ae,const Mesh& b,const Mesh& be) {
    return crosses(a,ae,b,be,false);
}
bool swept_crossing(const MeshList& cells,const MeshList& trial,const Mesh& wall) {
    if (cells.size()!=trial.size()) throw std::invalid_argument("trial cell list size changed");
    std::vector<Bounds> boxes;
    for (std::size_t i=0;i<cells.size();++i) boxes.push_back(swept_bounds(*cells[i],*trial[i]));
    for (auto [i,j]:scan_boxes(boxes))
        if (crosses(*cells[i],*trial[i],*cells[j],*trial[j],false)) return true;
    for (std::size_t i=0;i<cells.size();++i)
        if (crosses(*cells[i],*trial[i],wall,wall,true)) return true;
    return false;
}
}
