## 1-Dimensional Numerical Column Test

### General Structure
This repository contains an interactive, analytical [model]([https://jangei.github.io/1D_colum_test_analytical/](https://jangei.github.io/1D_column_test_numerical/)) that solves 1-dimensional solute transport through a cylinder, representative of a column experiment. 
The basic structure of the model is developed in Python, while the interactive nature of this model functions through JavaScript. 
This is realized through the Python library [bokeh](https://bokeh.org/) which provides customizable JavaScript callbacks in the Python language.
Contrary to the analytical model, which can be found in this [repository](https://github.com/JanGei/1D_column_test_analytical), the numerical model introduces non-linear sorption (Freundlich and Langmuir) sorption.
Therefore it is more complicated and contains an additional JS callback. 
The first callback contains the numerical model and computes the concentration profile on demand.
The second callback takes the results from the first callback and allows the user to navigate through the experiment time and space domain.

### Numerical Equations
![Numerical Equation](https://user-images.githubusercontent.com/99887101/195810360-752db23d-9f93-43ce-b3e8-3b0ce4e0ca3c.PNG)

To solve the reactive advection dispersion equation (ADE), an operator split is used in this model.
That means, that advection, reactive dispersion and sorption are accounted for in sequence.
The order within the sequence is identical as presented here:

#### Advection
The advective, dispersive, and reactive transport, including sorption, is computed through a fully implicit Crank-Nicholson scheme. 
To operate this scheme the domian is discretized in space into 100 identical cells. 
The time step size is set up in a way so that the concentration of one cell is passed to the next consecutive cell during one time step:  
<pre>
dt = dx / v  
with
- Distance between cell centers dx
- Seepage velocity v
</pre>

#### Reactive Dispersion
The set of equations shown above accounts for reactive dispersion.
All internal cells, *i.e.* all cells excluding the inlet (Top Cell) and the outlet (Bottom Cell), are governed by the same equation.
The notation is as follows:
<pre>
- Concentration c
- Cell index i (starting at 1 at the inlet cell)
- Time index k 
- Ammount of cells N
- Combined dispersive parameter s = D (dt/dx)<sup>2</sup>
  - Dispersion coefficient D
- Combined reactive parameter r = λ dt
  - First order rate constant λ
</pre>

#### Sorption

The model allows the user to choose between no sorption, linear sorption, and non-linear Langmuir and Freundlich sorption.
Therefore the model tracks the concentration in the aqueous phase, as well as in the solid phase.
Initializing the sorption step, the concentrations in each cell are summed up:
<pre>
Linear and Langmuir sorption:      C<sub>tot</sub> = s ρ<sub>d</sub> (1 - n) + n C<sub>aq</sub>
Freundlich sorption:               C<sub>tot</sub> = (ρ<sub>d</sub> (1 - n) K<sub>Fr</sub> C<sub>aq</sub><sup>n<sub>Fr</sub>-1</sup> + n) C<sub>aq</sub>
with
- Sorbed concentration s
- Total concentration C<sub>tot</sub>
- Linear Partitioning Coefficient K<sub>d</sub>
- Solid Density ρ<sub>d</sub>
- Porosity n
- Freundlich K K<sub>Fr</sub>
- Freundlich n n<sub>Fr</sub>
</pre>

Following that the solute is partitioned between the aqueous phase and the solid phase:
<pre>
Linear sorption:                  C<sub>aq</sub> = C<sub>tot</sub> / (K<sub>d</sub> (1-n) ρ<sub>d</sub> + n)
                                  s<sub>  </sub> = C<sub>aq</sub> K<sub>d</sub>
Langmuir sorption:                C<sub>aq</sub> = (-β + (β<sup>2</sup> - 4 n γ)<sup>0.5</sup>) / (2 n)
                                  s<sub>  </sub> = s<sub>max</sub> C<sub>aq</sub> / (K<sub>ads</sub> + C<sub>aq</sub>)
Freundlich sorption:              C<sub>aq</sub> = C<sub>tot</sub> / (ρ<sub>d</sub> (1-n) K<sub>Fr</sub> C<sub>aq</sub><sup>n<sub>Fr</sub>-1</sup> + n)
                                  s<sub>  </sub> = K<sub>Fr</sub> C<sub>aq</sub><sup>n<sub>Fr</sub></sup>
with
- Langmuir β = (1-n) ρ<sub>d</sub> s<sub>max</sub> + n K<sub>ads</sub> - C<sub>tot</sub>
- Langmuir γ = -C<sub>tot</sub> K<sub>ads</sub>
- Half saturation concentration K<sub>ads</sub>
- Specific sorption capacity s<sub>max</sub>
</pre>

#### User Interface (still in work)

![AnalyticalUI](https://user-images.githubusercontent.com/99887101/195622686-1e3190a3-8ecf-486a-8a09-4605bf15db6a.PNG)

The UI of the model is essentially divided into four parts.  
The plot in the upper left part shows the normed solute concentration as a function of space, *i.e.* length of the column, at a given point in time.  
The plot in the lower left part shows a breakthroughcurve, *i.e.* the normed solute concentration as  a function of time at a given location. 
This location is specified by the grey diamond in the upper plot, which can be dragged in order to change its position.  
The upper right part of the UI contains the widgets through which the user can change the models mode, geometry, and the solutes properties.
Initially the user can choose between a continuous injection and a pulse injection. The picture above depicts the model in its 'Continuous Injection' state. ![Numerical Equation](https://user-images.githubusercontent.com/99887101/195810319-2074bbae-f531-44a5-a996-1245fdc28b52.PNG)

If the mode of the model is switched to 'Pulse Injection', an additional slider appears, that controlls the duration of injection.
Furthermore, the user can choose between 'No Sorption' and 'Linear Sorption', as linear sorption of a solute can be modelled analytically. 
For this, the retardation factor R is computed:  
R  = 1 + (1-n)/n K<sub>d</sub> ρ<sub>d</sub>  
with:
- Porosity n
- Linear Partitioning Coefficient K<sub>d</sub>
- Solid Density ρ<sub>d</sub>

This retardataion factor is multiplied to the seepage velocity in order to model linear sorption. 
If the option 'Linear Sorption' is selected, three sliders, containing K<sub>d</sub> and ρ<sub>d</sub> appear.  
An extension of the problem to non-linear sorption can be found in this [repository]([https://bokeh.org/](https://github.com/JanGei/1D_column_test_numerical))

Below these options the 'Time' slider is found, through which the elapsed time within the simulation can be changed. 
Note that for practical reasons time is converted into pore volumes, but a conversion rate from pore volumes to hours is depicted in the slider title.
Following the 'Time' slider, all necessary and adjustable model parameters are displayed with a slider. 
Through the sliders the numerical values of the columns length and radius, the solutes first-order reaction coefficient, dispersion coefficient, and the setups flow rate and porosity can be altered.  

The lower right part of the UI contains a unit selection for the reaction, dispersion, and flow rate term. 
However, after selecting the unit of choice, the user needs to interact with the slider again, as this triggers the unit switch.  

Below the unit selection are two buttons that, when clicked, download a .csv file, containing x and y values of the upper and lower plot to the download folder of the computer.

#### Uncertainty

In experimental column setups the exact value of the first-order reaction coefficient and/or the dispersion coefficient is not known exactly.
To tackle this problem, the model allows the user to define a range with a lowest and highest possible value for both variables.
Within this range, a [latin hypercube](https://en.wikipedia.org/wiki/Latin_hypercube_sampling) is used to sample 100 parameter combinations containing both variables.
Thus, the model re-computes the aforeshown equation 100 times and computes the ensemble mean, as well as the upper and lower quartile and the minimum and maximum.
These values are also depicted in the figure in the upper left part of the UI.
Collapsing the range slider minimum and maximum into one point for both variables disables this option.

#### Future Work and Extensions

An extended (and numerical) version of this model can be found [here](https://github.com/JanGei/1D_column_test_numerical). 
It contains non-linear sorption (Freundlich and Langmuir). Further possible extensions and future work is discussed there.
