# ProjectHelium
Calculator and plotter for probability distribution for helium ion in 3D that numerically calculates data and plots them in either C++ or Python enviroment. Implemented different orbitals manipulation as well as slider that allows to slice the atom to see whats inside.

# Installation and setup

First clone the repository:

```bash
git clone https://github.com/Mikocoded/ProjectHelium.git
```

after cloning you will have to cmake the CMakeFile to build the executive:

```
cmake -B build && cmake --build build
```

after that you're ready to run the C++ code:

```
./build/scatter
```

# Overview

After running the scipt you will have to provide 3 quantum numbers. Then the code will save the csv file with plotting data. You will also be asked weather you want to open the plot in C++ Qt, if not just close the command.

You can use the Python PyVista plotter by simply running:

```
python3 Plotter.py
```

After that you can choose between a website plot or just a local plot that includes either scatter plot or volume just like the C++ plotter.
