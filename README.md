# What's this for?
This repo contains a modified version of Matej Krajnc's vertex model, which is used for cartilage simulation.
[Matej's research home page](https://matejkrajnc.splet.arnes.si/research/)

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
-  Run `main.o`. Enter `help` for more details.

# Example
<p align="center">
  <img src="flow.png" alt="example" width="300"/>
  <img src="simulation.gif" alt="example animation" width="300"/>
</p>
