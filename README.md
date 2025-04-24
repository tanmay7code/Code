\#Radme.md  
Description of the files and how to use the code  
**\#\# Contents**  
**\#Matlab**  
**\-**nsga2\_model1\_CT.m ( Carbon tax policy)  
\-nsga2\_model2\_CCO.m ( Carbon cap and offset policy)  
\-nsga2\_model3\_CCT.m( Carbon cap and trade policy)  
**\-**mopso\_model1\_CT.m (Carbon tax policy)  
\-mopso\_model2\_CCO.m  ( Carbon cap and offset policy)  
\-mopso\_model3\_CCT.m ( Carbon cap and trade policy)  
\-topsis.m (choose the best solution from pareto front)  
\-skill\_level.m (calculation for skill level for the next cycle)  
**\#Python**  
\-hyper\_volume (calculate the hypervolume from pareto front using Pymoo package)  
\-pareto\_front.csv  
\#\# **Requirements**  
\- MATLAB R2022b or later   
\- Optimization Toolbox  
\- PYTHON 3.12 or later   
\-Pymoo package

\#\# How to Run  
1\. Open the corresponding MATLAB script file based on the selected carbon emissions policy (for either NSGA-II or MOPSO).  
2\. Run the script to load the sample input data and execute the optimization algorithm.  
3\. For \*\*NSGA-II\*\*:  
   \- The output objective values will be saved in the workspace folders \`fp1\` (first objective) and \`fp2\` (second objective).  
4\. For \*\*MOPSO\*\*:  
   \- The output objective values will be saved in the workspace folders \`f1\` (first objective) and \`f2\` (second objective).  
5\. Plot the Pareto fronts over five consecutive cycles using the values from the first and second objectives.  
6\. Save the resulting Pareto front data into a file named \`pareto\_front.csv\`.  
7\. Run the Python script \`hyper\_volume.py\` to compute the hypervolume metric from the saved \`pareto\_front.csv\`.  
8\. To select the best solution from the Pareto front, open and run the \`topsis.m\` file in MATLAB, providing it with the objective values as input.

**\#\#\# Optional:**  
You may change parameters within to test different scenarios or datasets.

