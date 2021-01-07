# Semidiscrete-Central-Upwind-Schemes
Semidiscrete central upwind schemes for the compressible Euler Equations. Allows for solving 2D Reimann problems for Gas Dynamics without reimann problem solvers. 

# <sub>Mesh Reading Capability</sub>
The main python script is capable of reading in a mesh, defined by Gmsh's [3] .msh (v4) file format.
Mesh topology data required for the scheme is obtained from the file and stored in the Mesh class. 
<br /><br />
Following is output from the terminal to show this capability. Boundary cell indices have been 
chosen as the standard as they rely on all previous topology data to have been collected correctly.
<br /><br />
![](images/boundaryCellReading.png)

# <sub>Initialisation and Plotting</sub>
The Mesh is capable of reading in vectors of primitive variables and initialising the data accordingly.
Following image is the initial density distribution given the following initial conditions and as such, 
the initialisation process is verified.
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>1</sub> = [&rho;<sub>1</sub>, u<sub>1</sub>. v<sub>1</sub>, p<sub>1</sub>] = [1.0, 0.0, 0.0, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>2</sub> = [&rho;<sub>2</sub>, u<sub>2</sub>. v<sub>2</sub>, p<sub>2</sub>] = [0.5197, -0.7259, 0.0, 0.4]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>3</sub> = [&rho;<sub>3</sub>, u<sub>3</sub>. v<sub>3</sub>, p<sub>3</sub>] = [0.1072, -0.7259, -1.4045, 0.0439]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>4</sub> = [&rho;<sub>4</sub>, u<sub>4</sub>. v<sub>4</sub>, p<sub>4</sub>] = [0.2579, 0.0, -1.4045, 0.15]
<br />
<img src="images/initialisationVerification.png" width="480" height="568">

# <sub>Numerical Experiment Setups and Discussion</sub>
All the following numerical results have been obtained using a 100x100 mesh in a 1x1 square subected to an initial configuration listed below. 
All of the numerical experiments have been obtained from Ref. [2]. In all cases, the results lack precision due to the coarse 100x100 mesh. Simulation 
time has been a large limitation in going for more fine meshes. A refactoring of the scripts into a more suitable language like C, C++ or FORTRAN is required. 
<br /><br />


# <sub><sub>Configuration 3.</sub></sub>
Initial Conditions: 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>1</sub> = [&rho;<sub>1</sub>, u<sub>1</sub>. v<sub>1</sub>, p<sub>1</sub>] = [1.5, 0.0, 0.0, 1.5]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>2</sub> = [&rho;<sub>2</sub>, u<sub>2</sub>. v<sub>2</sub>, p<sub>2</sub>] = [0.5323, 1.206, 0.0, 0.3]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>3</sub> = [&rho;<sub>3</sub>, u<sub>3</sub>. v<sub>3</sub>, p<sub>3</sub>] = [0.138, 1.206, 1.206, 0.029]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>4</sub> = [&rho;<sub>4</sub>, u<sub>4</sub>. v<sub>4</sub>, p<sub>4</sub>] = [0.5323, 0.0, 1.206, 0.3]
<br /><br />
![](images/Config3.gif)

# <sub><sub>Configuration 4.</sub></sub>
Initial Conditions: 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>1</sub> = [&rho;<sub>1</sub>, u<sub>1</sub>. v<sub>1</sub>, p<sub>1</sub>] = [1.1, 0.0, 0.0, 1.1]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>2</sub> = [&rho;<sub>2</sub>, u<sub>2</sub>. v<sub>2</sub>, p<sub>2</sub>] = [0.5065, 0.8939, 0.0, 0.35]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>3</sub> = [&rho;<sub>3</sub>, u<sub>3</sub>. v<sub>3</sub>, p<sub>3</sub>] = [1.1, 0.8939, 0.8939, 1.1]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>4</sub> = [&rho;<sub>4</sub>, u<sub>4</sub>. v<sub>4</sub>, p<sub>4</sub>] = [0.5065, 0.0, 0.8939, 0.35]
<br /><br />
![](images/Config4.gif)

# <sub><sub>Configuration 5.</sub></sub>
Initial Conditions: 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>1</sub> = [&rho;<sub>1</sub>, u<sub>1</sub>. v<sub>1</sub>, p<sub>1</sub>] = [1.0, -0.75, -0.75, 1.0] 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>2</sub> = [&rho;<sub>2</sub>, u<sub>2</sub>. v<sub>2</sub>, p<sub>2</sub>] = [2.0, -0.75, 0.5, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>3</sub> = [&rho;<sub>3</sub>, u<sub>3</sub>. v<sub>3</sub>, p<sub>3</sub>] = [1.0, 0.75, 0.5, 1.0] 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>4</sub> = [&rho;<sub>4</sub>, u<sub>4</sub>. v<sub>4</sub>, p<sub>4</sub>] = [3.0, 0.75, -0.5, 1.0]
<br /><br />
![](images/Config5.gif)

# <sub><sub>Configuration 6.</sub></sub>
Initial Conditions: 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>1</sub> = [&rho;<sub>1</sub>, u<sub>1</sub>. v<sub>1</sub>, p<sub>1</sub>] = [1.0, 0.75, -0.5, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>2</sub> = [&rho;<sub>2</sub>, u<sub>2</sub>. v<sub>2</sub>, p<sub>2</sub>] = [2.0, 0.75, 0.5, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>3</sub> = [&rho;<sub>3</sub>, u<sub>3</sub>. v<sub>3</sub>, p<sub>3</sub>] = [1.0, -0.75, 0.5, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>4</sub> = [&rho;<sub>4</sub>, u<sub>4</sub>. v<sub>4</sub>, p<sub>4</sub>] = [3.0, -0.75, -0.5, 1.0]
<br /><br />
![](images/Config6.gif)

# <sub><sub>Configuration 8.</sub></sub>
Initial Conditions: 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>1</sub> = [&rho;<sub>1</sub>, u<sub>1</sub>. v<sub>1</sub>, p<sub>1</sub>] = [0.5197, 0.1, 0.1, 0.4]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>2</sub> = [&rho;<sub>2</sub>, u<sub>2</sub>. v<sub>2</sub>, p<sub>2</sub>] = [1.0, -0.6259, 0.1, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>3</sub> = [&rho;<sub>3</sub>, u<sub>3</sub>. v<sub>3</sub>, p<sub>3</sub>] = [0.8, 0.1, 0.1, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>4</sub> = [&rho;<sub>4</sub>, u<sub>4</sub>. v<sub>4</sub>, p<sub>4</sub>] = [1.0, 0.1, -0.6259, 1.0]
<br /><br />
![](images/Config8.gif)

# <sub><sub>Configuration 9.</sub></sub>
Initial Conditions: 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>1</sub> = [&rho;<sub>1</sub>, u<sub>1</sub>. v<sub>1</sub>, p<sub>1</sub>] = [1.0, 0.0, 0.0, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>2</sub> = [&rho;<sub>2</sub>, u<sub>2</sub>. v<sub>2</sub>, p<sub>2</sub>] = [0.5197, -0.7259, 0.0, 0.4]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>3</sub> = [&rho;<sub>3</sub>, u<sub>3</sub>. v<sub>3</sub>, p<sub>3</sub>] = [0.1072, -0.7259, -1.4045, 0.0439]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>4</sub> = [&rho;<sub>4</sub>, u<sub>4</sub>. v<sub>4</sub>, p<sub>4</sub>] = [0.2579, 0.0, -1.4045, 0.15]
<br /><br />
![](images/Config9.gif)

# <sub><sub>Configuration 11.</sub></sub>
Initial Conditions: 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>1</sub> = [&rho;<sub>1</sub>, u<sub>1</sub>. v<sub>1</sub>, p<sub>1</sub>] = [1.0, 0.1, 0.0, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>2</sub> = [&rho;<sub>2</sub>, u<sub>2</sub>. v<sub>2</sub>, p<sub>2</sub>] = [0.5313, 0.8276, 0.0, 0.4]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>3</sub> = [&rho;<sub>3</sub>, u<sub>3</sub>. v<sub>3</sub>, p<sub>3</sub>] = [0.8, 0.1, 0.0, 0.4]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>4</sub> = [&rho;<sub>4</sub>, u<sub>4</sub>. v<sub>4</sub>, p<sub>4</sub>] = [0.5313, 0.1, 0.7276, 0.4]
<br /><br />
![](images/Config11.gif)

# <sub><sub>Configuration 12.</sub></sub>
Initial Conditions: 
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>1</sub> = [&rho;<sub>1</sub>, u<sub>1</sub>. v<sub>1</sub>, p<sub>1</sub>] = [0.5313, 0.0, 0.0, 0.4]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>2</sub> = [&rho;<sub>2</sub>, u<sub>2</sub>. v<sub>2</sub>, p<sub>2</sub>] = [1.0, 0.7276, 0.0, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>3</sub> = [&rho;<sub>3</sub>, u<sub>3</sub>. v<sub>3</sub>, p<sub>3</sub>] = [0.8, 0.0, 0.0, 1.0]
<br />
&nbsp;&nbsp;&nbsp;&nbsp; W<sub>4</sub> = [&rho;<sub>4</sub>, u<sub>4</sub>. v<sub>4</sub>, p<sub>4</sub>] = [1.0, 0.0, 0.7276, 1.0]
<br /><br />
![](images/Config12.gif)

# <sub>References</sub>
[1] Alexander Kurganov; Eitan Tadmor (2002). <em>Solution of two-dimensional Riemann problems for gas dynamics without Riemann problem solvers.</em> , 18(5), 584–608. doi:10.1002/num.10025 <br />
[2] Kurganov, Alexander; Noelle, Sebastian; Petrova, Guergana  (2001). <em>Semidiscrete Central-Upwind Schemes for Hyperbolic Conservation Laws and Hamilton--Jacobi Equations.</em> SIAM Journal on Scientific Computing, 23(3), 707–740. doi:10.1137/s1064827500373413     
[3] Christophe Geuzaine; Jean-François Remacle (2009). <em>Gmsh: A 3-D finite element mesh generator with built-in pre- and post-processing facilities.</em> , 79(11), 1309–1331. doi:10.1002/nme.2579
