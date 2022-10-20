## 1-Dimensional Numerical Column Test

### General Structure
This repository contains an interactive, analytical [model](https://jangei.github.io/1D_column_test_numerical/) that solves 1-dimensional solute transport through a cylinder, representative of a column experiment. 
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

#### User Interface 

![Numerical UI](https://user-images.githubusercontent.com/99887101/196895783-cfc204af-c264-47c6-9900-2bffbf21ed7a.PNG)

The UI of the model is essentially divided into four parts.  
The plot in the upper left part shows the normed solute concentration as a function of space, *i.e.* length of the column, at a given point in time.  
The plot in the lower left part shows a breakthroughcurve, *i.e.* the normed solute concentration as  a function of time at a given location. 
This location is specified by the grey diamond in the upper plot, which can be dragged in order to change its position.  
The upper right part of the UI contains the widgets through which the user can change the models mode, geometry, and the solutes properties.
Initially the user can choose between a continuous injection and a pulse injection. The picture above depicts the model in its 'Continuous Injection' state.

If the mode of the model is switched to 'Pulse Injection', an additional slider appears, that controlls the duration of injection.
Furthermore, the user can choose between 'Linear Sorption', 'Langmuir Sorption', and 'Freundlich Sorption' (see above)

Below these options the 'Time' slider is found, through which the elapsed time within the simulation can be changed. 
Note that for practical reasons time is converted into pore volumes, but a conversion rate from pore volumes to hours is depicted in the slider title.
Following the 'Time' slider, all necessary and adjustable model parameters are displayed with a slider. 
Through the sliders, the numerical values of the columns length and radius, the solutes first-order reaction coefficient, dispersion coefficient, and the setup's flow rate and porosity can be altered.  

The lower right part of the UI contains sliders for adjusting the sorption relevant parameters. Depending on the type of sorption selected,
these sliders will change. Below that, a unit selection for the reaction, dispersion, and flow rate slider is found.
However, after selecting the unit of choice, the user needs to interact with the slider again, as this triggers the unit switch.  

Lastly, two buttons are found that, when clicked, download a .csv file, containing x and y values of the upper and lower plot to the download folder of the computer.

#### Future Work and Extensions

As mentioned above, an operator split is used in this model to account for the different processes individually. This yields the modular character of the model, meaning that it takes approximately 5 minutes to change the code to omit sorption completely from the computation.
Contrary to that a modular extension of the model is feasable and practical. In this section four possible extensions are presented:

<pre>
n-Molecular and Microbial Degradation

The model in its current state allows the user to specify a first-order reaction rate. This is equivalent to photolytic degradation 
of first order or radioactive decay. However, there are more, and more complicated ways of compound degradation. 
One way of accounting for a more flexible and more representative setup is the addition of several compounds that react in a 
bi-, tri-, or higher molecular reaction. Adding complexity to the system, the user could specify the interactions between differnt 
compounds. 
Furthermore, especially organic compunds are subject to microbial degradataion. Different models of microbial growth, thus 
compund degradation, are described, e.g. Monod-kinetics. This could also be implemented into the model, which would
simultaneously allow to track the ammount of micro-organisms, thus modelling microbial growth.
</pre>

![redox zone](https://user-images.githubusercontent.com/99887101/196898647-2dd437dd-12ac-4eff-b7fb-b7e2a680ccdf.png) <sup>[1]</sup>

<pre>
Redox Zone

Adding to the idea of microbial growth, the redox zonation (as depicted above) could be investigated. If, initially a soil column 
is filled with oxygen saturated water and a microbial degradable compound is introduced into to the system, oxygenic microbes
will use the oxygen until it is depleted to degrade and feed of the introduced compund. As soon as the oxygen is depleted in 
one cell, the next best energy metabolism takes over until its resources are depleted and so on. This leads to a zonation along
the axis of transport. By keeping track of different constituents of a typical soil, this redox zonation could be implemented 
into the model.
</pre>

<pre>
Sampling Time

An important part of laboratory work is smapling. In the case of column experiments the order of magnitude of the time needed to
get a breaktrhough curve can differ greatly. Therefore it is important to know, when to sample the columns. If a setup is thought
to last 10 pore volumes it would not be reccommended to sample ever 0.5 pore volumes, as one would miss cruicial information around
the actual breaking through of the curve. Therefore un-evenly distributed sampling times, with a maximum around the inflection point
are advised. This could be implemented by taking the data of the lower plot and with the help of statistics the model could propose
a set of objectively good sampling times to the user. However, it would have to be implemented in a way to account for variability
in the setup, as well as for uncertainty in the users input.
</pre>

<pre>
Inverse Modelling

While the aforementioned suggestions can be implemented modular, the idea of inverse modellling requires a complete overhaul of 
the model. The idea is that, after obtaining measurement data, the user is able to feed the model with tabular data, which are
then analyzed. The user would still remain in controll of the model modes, i.e. requiring bi-molecular reaction, monod kinetics
and Langmuir sorption and the model would in turn, propose a set of parameters that objectively solve the problem best. However,
the structure of the code is not conditioned to be easily invertable, but an inverse application of the model at hand would ease
the evaluation of such experiments.
</pre>

#### Bugs

As no code is perfect, this one also contains bugs. Here is a list of known ones:

<pre>
Column Length Bug

After changing the length of the column, an additional interaction with the time slide will change the x-axis of the upper plot,
which yields false information, as the data displayed fits to a column of the former length.
</pre>

<pre>
BTC location bug

Within the JS callback, the postioion of the grey diamong in the upper plot that in turn determines the location of the BTC in the
lower plot, is forced on to the x-axis during each callback event, i.e. changing the diamonds location. Most of the times, this
works as intended but in other times the correction of the y-coordinate of diamond seems to be skipped, resulting in an off-(y=0)-
position of the diamond. Interacting with the diamond again fixes this problem.
</pre>

#### References

<sup>[1]</sup> Vodyanitskii, Yu. (2015). The influence of Fe(III) on oil biodegradation in excessively moistened soils and sediments. Eurasian Soil Science. 48. 10.1134/S1064229315070121. 
