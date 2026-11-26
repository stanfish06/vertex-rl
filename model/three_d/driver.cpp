#include "driver.h"
#include "config.h"
#include "export.h"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <iostream>
#include <optional>
#include <stdexcept>

namespace vertex3d {
namespace {
Json status(const char* name) {
    Json result(json_object_new_object(),json_object_put);
    json_object_object_add(result.get(),"status",json_object_new_string(name));
    return result;
}

Json error(const std::string& message) {
    Json result=status("error");
    json_object_object_add(result.get(),"message",json_object_new_string(message.c_str()));
    return result;
}

class Session {
    int threads_;
    Path output_;
    std::unique_ptr<Solver> solver_;
    Collection collection_;
    double limit_=0;
    bool failed_=false;

    const char* phase() const {
        if (failed_) return "solver_error";
        if (solver_->cells().empty()) return "terminated";
        if (limit_ == solver_->time()) return "truncated";
        return "ok";
    }
    Json report() const {
        auto result=state(*solver_);
        json_object_object_add(result.get(),"status",json_object_new_string(phase()));
        return result;
    }
    void active() const {
        if (failed_ || solver_->time() >= limit_ || solver_->cells().empty()) {
            throw std::invalid_argument("reset required after episode end or solver failure");
        }
    }
public:
    Session(int threads, const Path& output) : threads_(threads),output_(output),collection_(output) {}

    Json handle(json_object* request) {
        std::string op=string(field(request,"op"));
        if (op=="ping" || op=="close") {
            keys(request,{"op"});
            return status(op=="ping" ? "ready" : "closed");
        }
        if (op=="reset") {
            keys(request,{"op","seed","time_limit","experiment"});
            unsigned seed=seed_value(field(request,"seed"));
            double limit=number(field(request,"time_limit"));
            if (!(limit>0)) throw std::invalid_argument("time_limit must be positive");
            auto* spec=field(request,"experiment");
            Path base=".";
            Json file(nullptr,json_object_put);
            if (json_object_is_type(spec,json_type_string)) {
                Path path=string(spec);
                file=read_json(path); spec=file.get(); base=path.parent_path();
            }
            auto next=std::make_unique<Solver>(experiment(spec,base),seed,threads_);
            solver_=std::move(next);
            limit_=limit;
            failed_=false;
            collection_=Collection(output_);
            return report();
        }
        if (!solver_) throw std::invalid_argument("reset required");
        if (op=="measure") { keys(request,{"op"}); return report(); }
        if (op=="export") {
            keys(request,{"op","frame_id"});
            int frame=integer(field(request,"frame_id"),1);
            collection_.write(*solver_,frame);
            auto result=status("ok");
            json_object_object_add(result.get(),"frame_id",json_object_new_int(frame));
            json_object_object_add(result.get(),"time",json_object_new_double(solver_->time()));
            return result;
        }
        if (op=="divide") {
            keys(request,{"op","cell_id","normal"});
            active();
            bool accepted=solver_->divide(integer(field(request,"cell_id"),1),vector3(field(request,"normal")));
            auto result=report();
            json_object_object_add(result.get(),"division_accepted",json_object_new_boolean(accepted));
            return result;
        }
        if (op!="step") throw std::invalid_argument("unknown operation " + op);
        keys(request,{"op","deadline"});
        active();
        double deadline=number(field(request,"deadline"));
        if (deadline<=solver_->time()) throw std::invalid_argument("deadline must be later than time");
        try {
            solver_->advance(std::min(deadline,limit_));
        } catch (const std::exception& failure) {
            failed_=true;
            auto result=report();
            json_object_object_add(result.get(),"message",json_object_new_string(failure.what()));
            return result;
        }
        return report();
    }
};

void worker(int threads, const Path& output) {
    Session session(threads,output);
    std::string line;
    while (std::getline(std::cin,line)) {
        Json reply(nullptr,json_object_put);
        try { reply=session.handle(parse(line).get()); }
        catch (const std::exception& ex) { reply=error(ex.what()); }
        std::cout << json_object_to_json_string_ext(reply.get(),JSON_C_TO_STRING_PLAIN) << std::endl;
        if (string(field(reply.get(),"status"))=="closed") break;
    }
}

struct TimedCut { double time; int id; Vec3 normal; };

void run(const Path& path, const Path& output, int threads) {
    auto file=read_json(path);
    auto* obj=file.get();
    unsigned seed=field(obj,"seed",false) ? seed_value(field(obj,"seed")) : 646;
    double limit=number(field(obj,"time_limit"));
    double cadence=number(field(obj,"output_interval"));
    if (limit<=0 || cadence<=0) throw std::invalid_argument("time_limit and output_interval must be positive");
    std::vector<TimedCut> cuts;
    if (auto* events=field(obj,"divisions",false)) {
        if (!json_object_is_type(events,json_type_array)) throw std::invalid_argument("divisions must be an array");
        for (size_t i=0; i<json_object_array_length(events); ++i) {
            auto* d=json_object_array_get_idx(events,i);
            keys(d,{"time","cell_id","normal"});
            TimedCut cut{number(field(d,"time")),integer(field(d,"cell_id"),1),vector3(field(d,"normal"))};
            if (cut.time<0 || cut.time>limit) throw std::invalid_argument("division time outside run");
            cuts.push_back(cut);
        }
        std::stable_sort(cuts.begin(),cuts.end(),[](const auto& a,const auto& b){return a.time<b.time;});
    }
    Solver solver(experiment(obj,path.parent_path()),seed,threads);
    Collection collection(output);
    size_t event=0, sample=0;
    int frame=0;
    while (true) {
        while (event<cuts.size() && cuts[event].time==solver.time()) {
            auto& cut=cuts[event++];
            if (!solver.divide(cut.id,cut.normal)) throw std::runtime_error("scheduled division rejected");
        }
        collection.write(solver,++frame);
        if (solver.time()==limit || solver.cells().empty()) break;
        if (sample*cadence<=solver.time()) ++sample;
        double deadline=std::min(limit,sample*cadence);
        if (event<cuts.size()) deadline=std::min(deadline,cuts[event].time);
        solver.advance(deadline);
    }
    std::cout << "Wrote " << frame << " frames to " << (output/"tissue.pvd") << '\n';
}
}

int driver(int argc, char** argv) {
    try {
        int threads=1;
        std::vector<std::string> args;
        for (int i=1; i<argc; ++i) {
            std::string arg=argv[i];
            if (arg=="--threads" && i+1<argc) {
                std::string value=argv[++i];
                auto result=std::from_chars(value.data(),value.data()+value.size(),threads);
                if (result.ec!=std::errc() || result.ptr!=value.data()+value.size() || threads<=0) {
                    throw std::invalid_argument("--threads must be a positive integer");
                }
            } else args.push_back(arg);
        }
        if (args.size()>=1 && args[0]=="worker" && args.size()<=2) {
            worker(threads,args.size()==2 ? args[1] : "output");
        } else if (args.size()==3 && args[0]=="run") {
            run(args[1],args[2],threads);
        } else throw std::invalid_argument("usage: vertex3d run EXPERIMENT.json OUTPUT_DIR | worker [OUTPUT_DIR] [--threads N]");
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 1;
    }
}
}
