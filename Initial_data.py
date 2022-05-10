from cmath import nan
from logging import PlaceHolder
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource,FuncTickFormatter, CustomJS, Slider, Panel, Range1d, Tabs, Button, RangeSlider, RadioButtonGroup, PointDrawTool
from bokeh.plotting import Figure, output_file, show
from bokeh.events import Tap, Pan
import numpy as np
import math
from math import exp, erfc, sqrt, ceil, floor
from bokeh.embed import components
import os
import sys


# Setting working sirectory to current folder
os.chdir(os.path.dirname(sys.argv[0]))

# Some functions
def get_gamma(reac,Dis,sep_vel):
  res = sqrt(1 + 4 * reac * Dis / sep_vel**2)
  return res

def transport(c,c_in):
  for j in range(len(c)):
    c[len(c)-j-1] = c[len(c)-j-2]
  c[0] = c_in
  
  return c

def reactive_dispersion(c_arr,Disp, dx, dt, nX, A_cn, b_cn, k, rf):
  dt_sub = dt / rf
  s = Disp*dt_sub/dx**2
  r = -k * dt_sub 

  for j in range(rf):
    for i in range(nX):
      # Left hand side matrix A_CN
      if (i > 0 and i < nX-1):
        A_cn[i][i-1]  = - s/2
        A_cn[i][i]    = 1 + s - r
        A_cn[i][i+1]  = - s/2
      elif (i == 0):
        A_cn[i][i]    = 1 + s/2 - r
        A_cn[i][i+1]  = - s/2
      elif (i == nX-1):
        A_cn[i][i-1]  = - s/2 
        A_cn[i][i]    = 1 + s/2 - r
      # Right hand side vector 
      if (i > 0 and i < nX-1):
        b_cn[i]   = s/2*c_arr[i-1] + (1-s)* c_arr[i] + s/2*c_arr[i+1]
      elif (i == 0):
        b_cn[i]    = (1-s/2)*c_arr[i] + s/2*c_arr[i+1]
      elif (i == nX-1):
        b_cn[i]    = (1-s/2)*c_arr[i] + s/2*c_arr[i-1]

  res = np.linalg.solve(A_cn,b_cn)
  return res

def reactive_dispersion_implicit(c_arr,Disp, dx, dt, nX, A_cn, k, rf):
  dt_sub = dt / rf
  s = Disp*dt_sub/dx**2
  r = -k * dt_sub 

  for j in range(rf):
    for i in range(nX):
      # Left hand side matrix A_CN
      if (i > 0 and i < nX-1):
        A_cn[i][i-1]  = - s
        A_cn[i][i]    = 1 + 2*s - r
        A_cn[i][i+1]  = - s
      elif (i == 0):
        A_cn[i][i]    = 1 + s - r
        A_cn[i][i+1]  = - s
      elif (i == nX-1):
        A_cn[i][i-1]  = - s
        A_cn[i][i]    = 1 + s - r

  res = np.linalg.solve(A_cn,c_arr)
  return res

# Initial slider parameters (min, max, step, value)
# Pore Volume [-]
pore_vol  = [np.log(0.001), np.log(7), (np.log(7)-np.log(0.001)) / 1000, np.log(0.5)]
# Column radius [m]
col_rad   = [0.005, 0.2, 0.0001, 0.05]
# Flow rate [ml/h]
flow      = [1, 50, 0.1, 10]
# Porosity [-]
poros     = [0.01, 1, 0.01, 0.5]
# Duration of pulse injection [min]
puls_inj  = [30, 43200, 30, 18000]
# Column length [m]
col_len   = [0.01, 0.5, 0.001, 0.2]
# Dispersion coefficient [ln(m2/h)]
disp      = [np.log(1e-6), np.log(1e-1), (np.log(1e-1)-np.log(1e-6))/300, np.log(1e-5)]
# First order reaction coefficient [ln(1/h)]
reac      = [np.log(1e-4), np.log(1), (np.log(1)-np.log(1e-4))/300, np.log(1e-2)]

# Initial slider parameter (min, max, step, value) for sorption
# Solid densitiy [kg/m3]
rho_s     = [2000, 3000, 1, 2650]
# Linear partitioning coefficient [m3/kg]
Kd        = [5e-5, 5e-3, 5e-5, 2e-3]
# Half saturation concentration [mol/m3]
K_ads     = [0.01, 10, 0.01, 1]
# Specific sorption capacity [mol/kg]
s_max     = [0.0001, 1, 0.0001, 0.1]
# Freundlich Sorption Parameter [mmol^(1-n)*L^(n) / kg]
K_Fr      = [0.01, 10, 0.01, 1]
# Frfeundlich n [-]
Fr_n      = [1, 2, 0.01, 1.3]

# Parameters for plot initialization + adjusting units
poros_ini   = poros[3]                                      # [-]
col_len_ini = col_len[3]                                    # [m]
col_rad_ini = col_rad[3]                                    # [m]
flow_ini    = flow[3]/1000/1000/3600                        # [m3/s]     
disp_ini    = exp(disp[3])                                  # [m2/s]
reac_ini    = exp(reac[3])                                  # [1/s]

# Pore space in the column [m3]
porespace   = col_len_ini * math.pi * col_rad_ini**2 * poros_ini
# Seepage velocity [m/s]
velocity_ini    = flow_ini/(col_rad_ini**2*math.pi*poros_ini)
# Time needed to fully flush the column
porevolume  = col_len_ini / velocity_ini
# Normed inlet concentration
c0 = 1

# number of spatial nodes in the domain
nX = 100
# spatial discretization
dx = col_len_ini/nX
# temporal discretization
dt = dx/velocity_ini
# Neumann-Number
Ne = 4 * disp_ini / velocity_ini / dx
# Refinement factor
rf = ceil(Ne)
# End time
tN = exp(pore_vol[1]) * porevolume
# number of temporal nodes
nT = floor(tN/dt)

# Initializing lists
c_tot_array   = np.zeros((nX,nT))
c_intermed    = np.zeros(nX)
A_cn          = np.zeros((nX,nX))
b_cn          = np.zeros(nX)
# Subdividing the column into 1000 equally long parts
x           = np.linspace(0,col_len[3],nX)
# Subdividing the duration of the experiment into nT equally long parts
PVspan      = np.linspace(exp(pore_vol[0]),exp(pore_vol[1]),nT)

# Numerical Crank-Nicholson-Scheme to determine concentrations
for i in range(nT):
  c_intermed = transport(c_intermed,c0)
  #c_intermed = reactive_dispersion(c_intermed,disp_ini,dx,dt,nX,A_cn,b_cn,reac_ini,rf)
  c_intermed = reactive_dispersion_implicit(c_intermed,disp_ini,dx,dt,nX,A_cn,reac_ini,rf)
  # no sorption here
  # store results
  for j in range(nX):
    c_tot_array[j][i] = c_intermed[j]
  print("We are "+ str(i+1) + " of " + str(nT) + " into the computation")


np.savetxt('Initial_Data_reaction_implicit.csv', c_tot_array, delimiter=';')

# first run took 55 mins
# fully implicit (700) 25 min
# fully implicit (500) 