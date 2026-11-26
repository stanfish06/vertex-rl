#include <charconv>
#include <cmath>
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>

#include "solver_2d.h"
#ifdef VERTEX_WITH_3D
#include "three_d/driver.h"
#endif
#include <omp.h>
#include <iostream>
#include <string>
//*******************FUNCTION****************************
void sample_cli(float displacement) {
    if (Time >= TimeCount) {
        export_scheduled_frame();
        printf("Time=%g\tw=%.16g\t||\t%.16g\n", Time, wA+wS+wL, displacement);
        TimeCount += outFreq;
    }
}

void parseAndSet(const std::string &input) {
    std::istringstream iss(input);
    std::string command, parameter;
    float value;
    if (!(iss >> command >> parameter >> value) || !std::isfinite(value)) {
        std::cerr << "Expected set NAME NUMBER" << std::endl;
        return;
    }
    if (auto* option = find_parameter(parameter)) option->set(value);
    else if (parameter == "tMAX") tMAX = value;
    else if (parameter == "outFreq") outFreq = value;
    else if (parameter == "max_num_edges") max_num_edges = int(value);
    else if (parameter == "min_division_n_edges") min_division_n_edges = int(value);
    else if (parameter == "min_exist_n_edges") min_exist_n_edges = int(value);
    else if (parameter == "number_of_threads") {
        if (value < 1 || floor(value) != value || double(value) > std::numeric_limits<int>::max()) {
            std::cerr << "Thread count must be a positive integer" << std::endl;
            return;
        }
        set_solver_threads(static_cast<int>(value));
    }
    else std::cout << "Unknown parameter: " << parameter << std::endl;
    std::cout << std::endl;
}

void parameter_list() {        
printf(
    " Name                           Description                                Current value   \n"
    "  1. k0                         Cell-cell repulsion                        %.3g            \n"
    "  2. kS                         Cell-ECM repulsion                         %.3g            \n"
    "  3. lS                         ECM width                                  %.3g            \n"
    "  4. kA                         Area elasticity                            %.3g            \n"
    "  5. e_g0                       Line tension                               %.3g            \n"
    "  6. tMAX                       Final simulation time                      %.3g            \n"
    "  7. outFreq                    Output time interval                       %.3g            \n"
    "  8. maxDR                      Maximum displacement per step              %.3g            \n"
    "  9. lthUP                      Upper edge length bound                    %.3g            \n"
    " 10. lthDOWN                    Lower edge length bound                    %.3g            \n"
    " 11. lALERT                     Cell-cell repulsion width                  %.3g            \n"
    " 12. max_num_edges              Upper edge number bound                    %d              \n"
    " 13. min_division_n_edges       Lower edge number bound at division        %d              \n"
    " 14. min_exist_n_edges          Lower edge number bound normally           %d              \n"
    " 15. number_of_threads          Program parallelism                        %d		\n"
    " 16. resting_CV                 Resting zone cell volume                   %.3g 		\n"
    " 17. proliferative_CV           Proliferative zone cell volume             %.3g 		\n"
    " 18. hypertrophic_CV            Hypertrophic zone cell volume              %.3g 		\n"
    " 19. ossification_CV            Ossification zone cell volume              %.3g 		\n"
    " 20. resting_ECM                Resting zone ECM scaling                   %.3g 		\n"
    " 21. proliferative_ECM          Proliferative zone ECM scaling             %.3g 		\n"
    " 22. hypertrophic_ECM           Hypertrophic zone ECM scaling              %.3g 		\n"
    " 23. ossification_ECM           Ossification zone ECM scaling              %.3g 		\n"
    " 24. resting_T                  Resting zone line tension scaling          %.3g 		\n"
    " 25. proliferative_T            Proliferative zone line tension scaling    %.3g 		\n"
    " 26. hypertrophic_T             Hypertrophic zone line tension scaling     %.3g 		\n"
    " 27. ossification_T             Ossification zone line tension scaling     %.3g 		\n",
    k0, kS, lS, kA, e_g0, tMAX, outFreq, maxDR, lthUP, lthDOWN, lALERT, max_num_edges, min_division_n_edges, min_exist_n_edges, num_threads_force, CV[0], CV[1], CV[2], CV[3], ECM[0], ECM[1], ECM[2], ECM[3], EG[0], EG[1], EG[2], EG[3] 
);
}

void domain_setting() {
    if (v_D == nullptr) {
    	printf("Domain not set yet. Run specify domain first.\n");
	return;
    }
    std::cout << "Domain: " << std::endl;
    for (int i = 1; i <= v_D[0][0]; i++) {
	std::cout << "Domain edge " << i << " :" << std::endl;
	printf("Edge coords: (%f, %f)--(%f, %f)\n", v_D[e_D[i][0]][1], v_D[e_D[i][0]][2], v_D[e_D[i][1]][1], v_D[e_D[i][1]][2]);
	printf("Edge normal vectors: (%f, %f)\n", e_normal_unit_vector[e_D[i][0]][1], e_normal_unit_vector[e_D[i][0]][2]);
	printf("Edge tangent vectors: (%f, %f)\n", e_tangent_unit_vector[e_D[i][0]][1], e_tangent_unit_vector[e_D[i][0]][2]);
        std::cout << std::endl;
    }
    std::cout << std::endl;    
}

void interface_setting() {
    if (v_interface == nullptr) {
        printf("Interface not set yet. Run specify interface first.\n");
        return;
    }
    std::cout << "Interface: " << std::endl;
    for (int i = 1; i <= v_interface[0][0]; i++) {
        std::cout << "Interface " << i << " :" << std::endl;
        printf("Interface coords: (%f, %f)\n", v_interface[i][1], v_interface[i][2]);
        printf("Interface normal vectors: (%f, %f)\n", e_normal_unit_vector_interface[i][1], e_normal_unit_vector_interface[i][2]);
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


std::string checkUserInput() { 
    std::string userInput;
    //printf("waiting");
    std::cin >> userInput;
    std::cout << userInput << std::endl;
    return userInput;
}

//********************MAIN********************************
int main(int argc, char *argv[]){
#ifdef VERTEX_WITH_3D
    if (argc > 1 && std::string(argv[1]) == "--3d") return vertex3d::driver(argc - 1, argv + 1);
#endif


    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--reference") == 0) {
            reference_output = true;
        } else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            const char* value = argv[++i];
            const char* end = value + strlen(value);
            auto parsed = std::from_chars(value, end, placement_seed);
            if (parsed.ec != std::errc() || parsed.ptr != end) {
                fprintf(stderr, "Invalid unsigned placement seed: %s\n", value);
                return EXIT_FAILURE;
            }
        } else {
            fprintf(stderr, "Usage: %s [--seed UINT] [--reference]\n", argv[0]);
            return EXIT_FAILURE;
        }
    }

const char* helpManual = R"(
Usage: 2D vertex model for cartilage simulation

Commands:
  help               			Display this help message
  ls parameters      			List all parameters
  set [s]            			Set value for parameter s. Use <ls parameters> to see full list of parameters.
  specify domain     			Specify the domain of simulation. Remember to enter points in clockwise order. Don't enter same point twice. 
  ls domain				Print current domain setting.
  specify interface                     Specify the interface of simulation. Remember to enter points in clockwise order. Don't enter same point twice.
  ls interface                          Print current interface setting.
  run		     			Run simulation with current setting.    
  resume [path2snapshot]                Resume simulation with the setting provided in the file.
Examples:
)";

    printf("Press Ctrl+c or enter exit to quit. Enter help to get started.\n");
    Solver2D solver(false);
    char* buf;
    while ((buf = readline(">> ")) != nullptr) {
	std::string buf_str = buf;
	if (!buf_str.empty()) {
            add_history(buf);
        }
        free(buf);
        if (buf_str == "help") {
	    std::cout << helpManual << std::endl;
	} else if (buf_str == "exit") {
	    return 0;
	} else if (buf_str == "ls parameters") {
 	    //std::cout << parameterManual << std::endl;
	    parameter_list();
    	} else if (buf_str.substr(0, 3) == "set") {
	    parseAndSet(buf_str);
	} else if (buf_str == "specify domain") {
            if (read_and_set_domain() == 1) {
	    	printf("Domain setting fails.\n");
	    }
	} else if (buf_str == "ls domain") {
	    domain_setting();
	} else if (buf_str == "specify interface") {
	    if (read_and_set_interface() == 1) {
                printf("Interface setting fails.\n");
            }
        } else if (buf_str == "ls interface") {
            interface_setting();
        } else if (buf_str == "run" || buf_str.substr(0, 6) == "resume") {
            try {
                if (buf_str == "run") {
                    printf("Starting...\n");
                    solver.reset(placement_seed);
                } else {
                    if (reference_output) throw std::runtime_error("reference runs require fresh placement");
                    std::istringstream command(buf_str);
                    std::string verb, path;
                    command >> verb >> path;
                    auto suffix = path.find_last_of('_');
                    int frame = std::stoi(path.substr(suffix + 1));
                    solver.reset(placement_seed, path, frame);
                }
                solver.advance(tMAX, sample_cli);
            } catch (const std::exception& error) {
                std::cerr << error.what() << std::endl;
                return EXIT_FAILURE;
            }
            return 0;
        } else {
            printf("Unknown command\n");
        }
    }
    return 0;
}
