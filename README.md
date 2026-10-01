# ProjectHelium
Calculator and plotter for probability distribution for helium ion in 3D

The different orbitals can be achieved by changing quantum numbers that you will be asked to provide in terminal after running the cpp code. The code can be compiled by "g++ Numerical_calculation.cpp" command. Make sure you have g++ compiler installed. After that the code will be ran after "./a.out" command in terminal. After running the code a CSV file should appear in your folder. 

In order to plot the results the output CSV file will be plugged into python code. You can change the name of the file just remember to also change it in read file code in Plotter.py if you want to plot the results. Any CSV file in similar structure can be read allowing the code to plot any data. After running the Plotter.py you will be asked to select either local plotter or website one. If you chose the local you can choose either scatter plot or volume plot (Volume looks better however needs more computing so for weaker systems I would advise scatter plot).
