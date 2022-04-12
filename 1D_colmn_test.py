from cmath import nan
from logging import PlaceHolder
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource,FuncTickFormatter, CustomJS, Slider, Panel, Range1d, Tabs, Button, RangeSlider, RadioButtonGroup, PointDrawTool
from bokeh.plotting import Figure, output_file, show
from bokeh.events import Tap, Pan
import numpy as np
import math
from math import exp, erfc, sqrt
from bokeh.embed import components
from scipy.stats import qmc
import os
import sys

from zmq import PUSH

# Setting working sirectory to current folder
os.chdir(os.path.dirname(sys.argv[0]))

# Some functions
def get_gamma(reac,Dis,sep_vel):
  res = sqrt(1 + 4 * reac * Dis / sep_vel**2)
  return res

def getc_cont(x,c,vel,t,L1,L2,reac_l,reac_h,disp_l,disp_h):
  for j in range(len(L1)):
    r_intermed = reac_l + (reac_h-reac_l)*L1[j]
    D_intermed = disp_l + (disp_h-disp_l)*L2[j]
    gam_intermed = get_gamma(r_intermed,D_intermed,vel)

    for i in range(len(x)):

      if x[i] <= 0:
        c[j,i] = 1
      else: 
        c[j,i] = 1/2 * exp(x[i]*vel/(2*D_intermed))*(exp(-x[i]*vel*gam_intermed/(2*D_intermed))*erfc((x[i]-vel*t*gam_intermed)/sqrt(4*D_intermed*t))+exp(x[i]*vel*gam_intermed/(2*D_intermed))*erfc((x[i]+vel*t*gam_intermed)/sqrt(4*D_intermed*t)))
  
  return c

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

# Initial slider parameter (min, max, step, value) for numerical model
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

# number of nodes in the domain
num_n = 200

# Initial slider parameters (min, max, step, minvalue, maxvalue)
# Dispersion coefficient [ln(m2/h)]
disp      = [np.log(1e-6), np.log(1e-1), (np.log(1e-1)-np.log(1e-6))/300, np.log(1e-5), np.log(5e-5)]
# First order reaction coefficient [ln(1/h)]
reac      = [np.log(1e-4), np.log(1), (np.log(1)-np.log(1e-4))/300, np.log(1e-3), np.log(5e-3)]

# Subdividing the column into 1000 equally long parts
x           = np.linspace(-col_len[3]*0.02,col_len[3],num_n)
# Subdividing the duration of the experiment into 1000 equally long parts
PVspan      = np.linspace(exp(pore_vol[0]),exp(pore_vol[1]),num_n)

# Parameters for plot initialization + adjusting units
poros_ini   = poros[3]                                      # [-]
col_len_ini = col_len[3]                                    # [m]
col_rad_ini = col_rad[3]                                    # [m]
flow_ini    = flow[3]/1000/1000/3600                        # [m3/s]     
disp_ini    = np.mean([exp(disp[3]), exp(disp[4])])/3600    # [m2/s]
reac_ini    = np.mean([exp(reac[3]), exp(reac[4])])/3600    # [1/s]

# Pore space in the column [m3]
porespace   = col_len_ini * math.pi * col_rad_ini**2 * poros_ini
# Seepage velocity [m/s]
velocity_ini    = flow_ini/(col_rad_ini**2*math.pi*poros_ini)
# Time needed to fully flush the column
porevolume  = col_len_ini / velocity_ini
# Initial time point
time_ini    = porevolume * exp(pore_vol[3]) 
# Intial point for breakthrough curve
xBTC_ini    = col_len_ini/2
# Normed inlet concentration
c0 = 1
# Area [m2]
Area =  math.pi * col_rad_ini**2 

# Latin Hypercube Sampling 
sampler = qmc.LatinHypercube(d=1)
# 2 Samples for uncertain variables with 100 values each
Lcube1 = sampler.random(100)
Lcube2 = sampler.random(100)

# Concentration list
c      = np.zeros((len(Lcube1),len(x)))
# Concentration lists for the highest and smallest concentration at each point in the column
c_mean = np.empty((1,len(x)))
c_min  = np.empty((len(x)))
c_max  = np.empty((len(x)))
# Concnetration list for breakthrough curve
c_t = np.empty((len(PVspan)))

# Solving 1D transport equation in space
c1 = getc_cont(x,c,velocity_ini,time_ini,Lcube1,Lcube2,exp(reac[3])/3600,exp(reac[4])/3600,exp(disp[3])/3600,exp(disp[4])/3600)
c_mean = getc_cont(x,c_mean,velocity_ini,time_ini,[0.5],[0.5],exp(reac[3])/3600,exp(reac[4])/3600,exp(disp[3])/3600,exp(disp[4])/3600)

for j in range(len(x)):
  c_min[j] = np.min(c1[:,j])
  c_max[j] = np.max(c1[:,j])

# Gamma coefficient needed for BTC
gam_m = get_gamma(reac_ini,disp_ini,velocity_ini)

# Solving 1D transport equation in time
for j in range(len(PVspan)):
  c_t[j] = c0/2 * exp(xBTC_ini*velocity_ini/(2*disp_ini))*(exp(-xBTC_ini*velocity_ini*gam_m/(2*disp_ini))*erfc((xBTC_ini-velocity_ini*PVspan[j]*porevolume*gam_m)/sqrt(4*disp_ini*PVspan[j]*porevolume))+exp(xBTC_ini*velocity_ini*gam_m/(2*disp_ini))*erfc((xBTC_ini+velocity_ini*PVspan[j]*porevolume*gam_m)/sqrt(4*disp_ini*PVspan[j]*porevolume)))

# Defining data sources with dictionary
source1 = ColumnDataSource(data = dict(x=x, y=c_mean[0], ymin = c_min, ymax = c_max))
source2 = ColumnDataSource(data = dict(x2=PVspan, y2=c_t))
source3 = ColumnDataSource(data = dict(xBTC = [col_len[3]/2], yBTC = [0]))

# Concentrtation Plot
COLp = Figure(min_height = 400, y_axis_label='c(t)/c0',
            x_axis_label='x [m]',sizing_mode="stretch_both")
COLp.line('x', 'y', source = source1, line_width = 3, line_alpha = 0.6, line_color = 'red')
COLp.line('x', 'ymin', source = source1, line_width = 3, line_alpha = 0.6, line_color = 'black', line_dash = 'dashed')
COLp.line('x', 'ymax', source = source1, line_width = 3, line_alpha = 0.6, line_color = 'black', line_dash = 'dashed')
COLp.y_range = Range1d(-0.03, 1.05)
COLp.xaxis.axis_label_text_font_size = "17pt"
COLp.yaxis.axis_label_text_font_size = "17pt"
COLp.xaxis.major_label_text_font_size = "12pt"
COLp.yaxis.major_label_text_font_size = "12pt" 

# Initializing PointDrawTool
BTCcircle = COLp.diamond(x='xBTC',y = 'yBTC', source=source3 , size=18, color = 'black', fill_alpha=0.6 )
COLp.add_tools(PointDrawTool(renderers=[BTCcircle], num_objects = 1))
COLp.toolbar.active_multi = COLp.select_one(PointDrawTool)

# BTC plot -- Now consistent with JS side
BTCp = Figure(min_height = 400, y_axis_label='c(t)/c0',
            x_axis_label='Pore Volume',sizing_mode="stretch_both")
BTCp.line('x2', 'y2', source = source2, line_width = 3, line_alpha = 0.6, line_color = 'red')
BTCp.y_range = Range1d(0, 1.05)
BTCp.x_range = Range1d(0, 7)
BTCp.title = "Breakthrough Curve at x = 0.100 m (Drag diamond in upper plot to change)"
BTCp.xaxis.axis_label_text_font_size = "17pt"
BTCp.yaxis.axis_label_text_font_size = "17pt"
BTCp.xaxis.major_label_text_font_size = "12pt"
BTCp.yaxis.major_label_text_font_size = "12pt" 
BTCp.title.text_font_size = "13pt"

# sliders 
pore_vol_sl   = Slider(start=pore_vol[0], end=pore_vol[1], value=pore_vol[3], step=pore_vol[2], title="Pore Volume (1PV = " + str("%.2f" %(porevolume/3600)) + " h)",
                    format=FuncTickFormatter(code="""return (Math.exp(tick)).toFixed(4)+' [PV]'"""),sizing_mode="stretch_width")
pulse_inj_sl  = Slider(title = "Duration of Injection", start = puls_inj[0], end = puls_inj[1], step = puls_inj[2], value = puls_inj[3],
                    format=FuncTickFormatter(code="""return (tick/60).toFixed(1)+' [min]'"""),sizing_mode="stretch_width")
col_len_sl    = Slider(title = "Column length", start = col_len[0], end = col_len[1], step = col_len[2], value = col_len[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(3)+' [m]'"""),sizing_mode="stretch_width")
col_rad_sl    = Slider(title = "Column radius", start = col_rad[0], end = col_rad[1], step = col_rad[2], value = col_rad[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(3)+' [m]'"""),sizing_mode="stretch_width")
disp_sl       = RangeSlider(title = "Dispersion coefficient ", start = disp[0], end = disp[1], step = disp[2], value =(disp[3], disp[4]),
                    format=FuncTickFormatter(code="""return Math.exp(tick).toExponential(1).toString()+' [m2/h]'"""),sizing_mode="stretch_width")
reac_sl       = RangeSlider(title = "Reaction coefficient ", start = reac[0], end = reac[1], step = reac[2], value = (reac[3], reac[4]),
                    format=FuncTickFormatter(code="""return Math.exp(tick).toExponential(1).toString()+' [1/h]'"""),sizing_mode="stretch_width")
flow_sl       = Slider(title = "Flow Rate", start = flow[0], end = flow[1], step = flow[2], value = flow[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(1)+' [mL/h]'"""),sizing_mode="stretch_width")
poros_sl      = Slider(title = "Porosity", start = poros[0], end = poros[1], step = poros[2], value = poros[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(2)+' [-]'"""),sizing_mode="stretch_width")
# sliders for numerical model
rho_s_sl      = Slider(title = "Solid Density", start = rho_s[0], end = rho_s[1], step = rho_s[2], value = rho_s[3],
                    format=FuncTickFormatter(code="""return (tick/1000).toFixed(2)+' [kg/L]'"""),sizing_mode="stretch_width")
Kd_sl         = Slider(title = "Linear Partinioning Coefficient", start = Kd[0], end = Kd[1], step = Kd[2], value = Kd[3],
                    format=FuncTickFormatter(code="""return (tick*1000).toFixed(2)+' [L/kg]'"""),sizing_mode="stretch_width")
Kads_sl       = Slider(title = "Half Saturation Concentration", start = K_ads[0], end = K_ads[1], step = K_ads[2], value = K_ads[3],
                    format=FuncTickFormatter(code="""return (tick).toFixed(2)+' [mmol/L]'"""),sizing_mode="stretch_width")
s_max_sl      = Slider(title = "Specific Sorption Capacity", start = s_max[0], end = s_max[1], step = s_max[2], value = s_max[3],
                    format=FuncTickFormatter(code="""return (tick*1000).toFixed(1)+' [mmol/L]'"""),sizing_mode="stretch_width")
K_Fr_sl       = Slider(title = "Freundlich Sorption Parameter", start = K_Fr[0], end = K_Fr[1], step = K_Fr[2], value = K_Fr[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(2)+' [mmol^(1-n)*L^(n) / kg]'"""),sizing_mode="stretch_width")
Fr_n_sl       = Slider(title = "Freundlich n", start = Fr_n[0], end = Fr_n[1], step = Fr_n[2], value = Fr_n[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(2)+' [L/kg]'"""),sizing_mode="stretch_width")

Labels1 = ["Analytical Model", "Numerical Model"]
Labels2 = ["Continuous Injection", "Pulse Injection"]
Labels3 = ["Linear Sorption", "Langmuir Sorption", "Freundlich Sorption"]

rg_AN = RadioButtonGroup(labels = Labels1, active = 0)
rg_CP = RadioButtonGroup(labels = Labels2, active = 0)
rg_ST = RadioButtonGroup(labels = Labels3, active = 0)

computebutton = Button(label="Compute Numerical Model", button_type="success",sizing_mode="stretch_width")

with open ('callback_compute_numerical.js', 'r') as file2:
  cbCode_numerical = file2.read()
callback_compute_numerical = CustomJS(args=dict(
                            source1=source1,
                            source2 = source2,
                            source3 = source3,
                            pore_vol_sl = pore_vol_sl,
                            col_len_sl = col_len_sl,
                            reac_sl = reac_sl,
                            disp_sl = disp_sl,
                            col_rad_sl = col_rad_sl,
                            flow_sl = flow_sl,
                            poros_sl = poros_sl,
                            rho_s_sl = rho_s_sl,
                            Kd_sl = Kd_sl,
                            Kads_sl = Kads_sl,
                            s_max_sl = s_max_sl,
                            K_Fr_sl = K_Fr_sl,
                            Fr_n_sl = Fr_n_sl,
                            rg_CP = rg_CP,
                            rg_ST = rg_ST,
                            pulse_inj_sl = pulse_inj_sl,
                            BTCp = BTCp,
                            computebutton = computebutton
                            ),
    code=cbCode_numerical)

computebutton.js_on_click(callback_compute_numerical)

with open ('callback.js', 'r') as file1:
  cbCode = file1.read()
callback = CustomJS(args=dict(
                            source1=source1,
                            source2 = source2,
                            source3 = source3,
                            Lcube1 = Lcube1,
                            Lcube2 = Lcube2,
                            pore_vol_sl = pore_vol_sl,
                            col_len_sl = col_len_sl,
                            reac_sl = reac_sl,
                            disp_sl = disp_sl,
                            col_rad_sl = col_rad_sl,
                            flow_sl = flow_sl,
                            poros_sl = poros_sl,
                            rho_s_sl = rho_s_sl,
                            Kd_sl = Kd_sl,
                            Kads_sl = Kads_sl,
                            s_max_sl = s_max_sl,
                            K_Fr_sl = K_Fr_sl,
                            Fr_n_sl = Fr_n_sl,
                            rg_CP = rg_CP,
                            rg_AN = rg_AN,
                            rg_ST = rg_ST,
                            pulse_inj_sl = pulse_inj_sl,
                            BTCp = BTCp,
                            computebutton = computebutton
                            ),
    code=cbCode)

savebutton1 = Button(label="Save (Upper Plot)", button_type="success",sizing_mode="stretch_width")
savebutton1.js_on_click(CustomJS(args=dict(source=source1),code=open(os.path.join(os.path.dirname(__file__),"download.js")).read()))
savebutton2 = Button(label="Save (Lower Plot)", button_type="success",sizing_mode="stretch_width")
savebutton2.js_on_click(CustomJS(args=dict(source=source2),code=open(os.path.join(os.path.dirname(__file__),"download.js")).read()))
#credit: https://stackoverflow.com/questions/31824124/is-there-a-way-to-save-bokeh-data-table-content

# callbacks for widgets
pore_vol_sl.js_on_change('value', callback)
col_len_sl.js_on_change('value', callback)
col_rad_sl.js_on_change('value', callback)
reac_sl.js_on_change('value', callback)
disp_sl.js_on_change('value', callback)
flow_sl.js_on_change('value', callback)
poros_sl.js_on_change('value', callback)
pulse_inj_sl.js_on_change('value', callback)
rg_CP.js_on_change('active',callback)
rg_AN.js_on_change('active',callback)
rg_ST.js_on_change('active',callback)
# Make a button that needs to be pressed in order to compute numerical model?
COLp.js_on_event(Tap, callback)
COLp.js_on_event(Pan, callback)

layout1 = column(rg_AN,rg_CP,pore_vol_sl,col_len_sl,col_rad_sl,reac_sl,disp_sl,flow_sl,poros_sl,pulse_inj_sl,sizing_mode="stretch_width")
layout2 = column(rg_ST,rho_s_sl,Kd_sl,Kads_sl,s_max_sl,K_Fr_sl,Fr_n_sl,computebutton, savebutton1, savebutton2, sizing_mode="stretch_width")

pulse_inj_sl.visible = False
rho_s_sl.visible = False
Kd_sl.visible = False
rg_ST.visible = False
computebutton.visible = False
Kads_sl.visible = False
s_max_sl.visible = False
K_Fr_sl.visible = False
Fr_n_sl.visible = False

tab1 = Panel(child=COLp, title="ADRE")
plots = Tabs(tabs=[tab1])


# Work with template in order to modify html code
script, (div1, div2, div3, div4) = components((COLp,layout1,BTCp,layout2))

# Add hmtl lines
f = open("./themodel.js", "w")
script = "\n".join(script.split("\n")[2:-1])
f.write(script)
f.close()

# read in the template file
with open('template', 'r') as file :
  filedata = file.read()

# replace the target strings (object in html is "placeholder")
filedata = filedata.replace('+placeholder1+', div1)
filedata = filedata.replace('+placeholder2+', div2)
filedata = filedata.replace('+placeholder3+', div3)
filedata = filedata.replace('+placeholder4+', div4)

# write to html file
with open('index.html', 'w') as file:
  file.write(filedata)
