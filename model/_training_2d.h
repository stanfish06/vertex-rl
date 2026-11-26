// Cell-policy state belongs to one episode in the native worker process.
#pragma once
#include <map>
using PlanarPoint = std::array<double,2>;

using Ring = std::vector<PlanarPoint>;

double orient(PlanarPoint a, PlanarPoint b, PlanarPoint p) {
    return (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0]);
}

double polygon_area(const Ring& ring) {
    double twice=0;
    for (size_t i=2; i<ring.size(); ++i) twice+=orient(ring[0],ring[i-1],ring[i]);
    return twice/2;
}

bool simple_ring(const Ring& ring) {
    if (ring.size()<3 || !(std::abs(polygon_area(ring))>0)) return false;
    for (size_t i=0; i<ring.size(); ++i) {
        size_t j=(i+1)%ring.size();
        if (ring[i]==ring[j]) return false;
        for (size_t k=i+1; k<ring.size(); ++k) {
            size_t l=(k+1)%ring.size();
            if (k==j || l==i) continue;
            const auto a=ring[i], b=ring[j], c=ring[k], d=ring[l];
            if (std::max(a[0],b[0]) < std::min(c[0],d[0]) || std::max(c[0],d[0]) < std::min(a[0],b[0]) ||
                std::max(a[1],b[1]) < std::min(c[1],d[1]) || std::max(c[1],d[1]) < std::min(a[1],b[1])) continue;
            if (orient(a,b,c)*orient(a,b,d)<=0 && orient(c,d,a)*orient(c,d,b)<=0) return false;
        }
    }
    return true;
}

struct Zone2D {
    int id = 0, rule = 0; // legacy, all, layer, radial, interface, closed surface
    bool allow_growth = true, allow_stop = true, allow_division = true, allow_death = true;
    double growth_rate = 5, target_size = 1.5, division_rate = 0.8, death_rate = 0.1;
    double lower = 0, upper = 0;
    bool positive = true;
    std::vector<PlanarPoint> points;
};

struct Training2D {
    int mode = 0; // 0 legacy, 1 tissue, 2 cell
    std::vector<Zone2D> zones;
    std::map<int, size_t> index;
    std::vector<int> agents;
    std::vector<double> born;
    int next_id = 1;
    void birth(int slot, double time) {
        agents[slot] = next_id++;
        born[slot] = time;
    }
} training2d;

struct BoundaryDistance {
    double squared = std::numeric_limits<double>::infinity(), side = 0;
    PlanarPoint closest{};
    bool inside = false;
    void edge(PlanarPoint p, PlanarPoint a, PlanarPoint b) {
        double x = b[0]-a[0], y = b[1]-a[1];
        double dx = p[0]-a[0], dy = p[1]-a[1];
        double length2 = x*x+y*y;
        if (length2==0) return;
        double t = std::clamp((dx*x+dy*y)/length2,0.,1.);
        PlanarPoint q{a[0]+t*x,a[1]+t*y};
        double distance = (p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]);
        if (distance < squared) {
            squared = distance;
            closest = q;
            side = x*dy-y*dx;
        }
        if ((a[1]>p[1]) != (b[1]>p[1]) && p[0]<a[0]+dy*x/y) inside = !inside;
    }
};

BoundaryDistance wall_distance(PlanarPoint point) {
    BoundaryDistance result;
    for (int i=1; i<=e_D[0][0]; ++i) {
        const auto a=v_D[e_D[i][0]], b=v_D[e_D[i][1]];
        result.edge(point,{a[1],a[2]},{b[1],b[2]});
    }
    return result;
}

bool zone_matches(const Zone2D& z, PlanarPoint p) {
    if (z.rule==1) return true;
    auto a = z.points.front();
    if (z.rule==2 || z.rule==3) {
        double dx=p[0]-a[0], dy=p[1]-a[1];
        double value=z.rule==2 ? dx*z.points[1][0]+dy*z.points[1][1] : std::hypot(dx,dy);
        return value>=z.lower && value<=z.upper;
    }
    BoundaryDistance boundary;
    for (size_t i=1; i<z.points.size(); ++i) boundary.edge(p,z.points[i-1],z.points[i]);
    if (z.rule==5) boundary.edge(p,z.points.back(),a);
    return boundary.squared==0 || ((z.rule==4 ? boundary.side>=0 : boundary.inside)==z.positive);
}

int trained_zone(int slot) {
    if (!training2d.zones.front().rule) return get_subdomain(slot);
    PlanarPoint center{c_center[slot][1],c_center[slot][2]};
    for (const auto& z : training2d.zones) if (zone_matches(z,center)) return z.id;
    throw std::runtime_error("cell center outside configured zones");
}
