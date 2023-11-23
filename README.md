# What's this for?
This repo contains a modified version of Matej's vertex model, which is used for cartilage simulation.

# How to run?
-  Clone the repo
-  Go into `cartilage` folder
```
cd matej-vertex-model/cartilage
```
-  Compile the program
```
g++ main.c -O3 -o main.o -fopenmp -lreadline
```
-  Run the matlab script [simulation_initialization_and_visualization.mlx](matej-vertex-model/simulation_initialization_and_visualization.mlx) under `matej-vertex-model` folder
-  Follow the steps in the matlab script
