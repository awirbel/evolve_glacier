## Glacier Tests

This folder contains the files to run a Glacier example with using the stl file for creating a new mesh from the manuscript. 
The mesh geometry is created with data from Farinotti et al. (2017) based on Wilson et al. (2013) and Maussion et al. (2018), see the manuscript for details.

There are three options for time discretization: Crank-Nicholson and Backward Euler (folder: Standard_Setup) and a 2nd order TVD Runge Kutta scheme (folder: Runge_Kutta).
To choose the Backward-Euler scheme, change in cfg_datadir.py in line 37 TS='CN' to TS='BE'.

### File Structure

All python files are located in your working directory, initial mesh data and results are then stored in a subfolder called data that has to be created beforehand. There are many more options for visualization of results, just check the code and commented save commands to find what datasets are also simulated and can be used for visualization. 



