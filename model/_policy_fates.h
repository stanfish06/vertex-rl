// Policy cuts are constructed in temporary polygons before the parent is retired.
bool closed_daughter(const Ring& ring) {
    return ring.size()<198 && simple_ring(ring) && polygon_area(ring)>0;
}

Ring half_ring(const Ring& parent, PlanarPoint center, PlanarPoint normal, int sign, double gap) {
    auto distance=[&](PlanarPoint p) {
        return sign*((p[0]-center[0])*normal[0]+(p[1]-center[1])*normal[1])-gap;
    };
    Ring result;
    int crossings=0;
    auto a=parent.back();
    for (const auto& b : parent) {
        double da=distance(a), db=distance(b);
        if ((da>0)!=(db>0)) {
            double t=da/(da-db);
            result.push_back({static_cast<float>(a[0]+t*(b[0]-a[0])),
                              static_cast<float>(a[1]+t*(b[1]-a[1]))});
            ++crossings;
        }
        if (db>0) result.push_back(b);
        a=b;
    }
    if (crossings!=2) return {};
    result.erase(std::unique(result.begin(),result.end()),result.end());
    if (result.size()>1 && result.front()==result.back()) result.pop_back();
    return result;
}

// The cut leaves each daughter one long straight edge facing its sibling.
// Contact acts between a vertex and the interior of an opposing edge, and the
// two cut edges share their end points, so without interior vertices the
// siblings slide through each other. Add points spaced at most lthUP apart and
// offset by half a spacing between siblings so each projects inside an edge.
void refine_cut(Ring& ring, PlanarPoint center, PlanarPoint normal, int sign, double gap, double stagger) {
    auto distance=[&](PlanarPoint p) {
        return sign*((p[0]-center[0])*normal[0]+(p[1]-center[1])*normal[1])-gap;
    };
    const double tolerance=1e-4;
    for (size_t i=0; i<ring.size(); ++i) {
        auto a=ring[i], b=ring[(i+1)%ring.size()];
        if (std::abs(distance(a))>tolerance || std::abs(distance(b))>tolerance) continue;
        double length=std::hypot(b[0]-a[0],b[1]-a[1]);
        int pieces=std::max(2,static_cast<int>(std::ceil(length/lthUP)));
        Ring inserted;
        for (int k=0; k<pieces; ++k) {
            double t=(k+stagger)/pieces;
            if (t<=1e-6 || t>=1-1e-6) continue;
            inserted.push_back({static_cast<float>(a[0]+t*(b[0]-a[0])),
                                static_cast<float>(a[1]+t*(b[1]-a[1]))});
        }
        if (ring.size()+inserted.size()<198) ring.insert(ring.begin()+i+1,inserted.begin(),inserted.end());
        return;
    }
}

int materialize_ring(const Ring& ring, int slot) {
    std::vector<int> vertices;
    for (const auto& p : ring) vertices.push_back(make_vertex(p[0],p[1]));
    std::vector<int> edges(1,0);
    for (size_t i=0; i<vertices.size(); ++i) edges.push_back(make_edge(vertices[i],vertices[(i+1)%vertices.size()],0));
    slot=cell_sides(ring.size(),edges.data(),slot);
    make_cell(slot);
    calc_cell_center(slot);
    return slot;
}

bool divide_at_angle(int slot, double angle, double time) {
    Ring parent;
    for (int j=1; j<=c_vertices[slot][0]; ++j) {
        int id=c_vertices[slot][j];
        parent.push_back({v[id][1],v[id][2]});
    }
    if (!closed_daughter(parent)) return false;
    const auto o=parent.front();
    double sum=0, x=0, y=0;
    for (size_t i=2; i<parent.size(); ++i) {
        auto a=parent[i-1], b=parent[i];
        double weight=orient(o,a,b);
        sum+=weight;
        x+=weight*((a[0]-o[0])+(b[0]-o[0]));
        y+=weight*((a[1]-o[1])+(b[1]-o[1]));
    }
    PlanarPoint center{o[0]+x/(3*sum),o[1]+y/(3*sum)};
    PlanarPoint normal{-std::sin(angle),std::cos(angle)};
    double gap=std::min(0.001,0.001*std::sqrt(sum/2));
    auto one=half_ring(parent,center,normal,1,gap);
    auto two=half_ring(parent,center,normal,-1,gap);
    if (!closed_daughter(one) || !closed_daughter(two)) return false;
    if (training2d.mode) {
        refine_cut(one,center,normal,1,gap,0);
        refine_cut(two,center,normal,-1,gap,0.5);
        if (!closed_daughter(one) || !closed_daughter(two)) return false;
    }
    int extra=static_cast<int>(one.size()+two.size())-static_cast<int>(parent.size());
    if (extra>array_max-Nv+ids[1][0] || extra>array_max-Ne+ids[2][0] || Nc-ids[3][0]>=array_max) return false;

    float total=c_A0[slot];
    float first=total*(polygon_area(one)/(polygon_area(one)+polygon_area(two)));
    if (!(first>0 && total-first>0)) return false;
    int growing=c_growth[slot];
    std::vector<int> edges(c[slot]+1,c[slot]+c[slot][0]+1);
    std::vector<int> vertices(c_vertices[slot]+1,c_vertices[slot]+c_vertices[slot][0]+1);
    dissolve_cell(slot);
    for (int edge : edges) dissolve_edge(std::abs(edge),0);
    for (int vertex : vertices) dissolve_vertex(vertex);
    materialize_ring(one,slot);
    int other=materialize_ring(two,0);
    c_dead[slot]=c_dead[other]=0;
    c_growth[slot]=c_growth[other]=growing;
    c_A0[slot]=first; c_A0[other]=total-first;
    training2d.birth(slot,time);
    training2d.birth(other,time);
    return true;
}
