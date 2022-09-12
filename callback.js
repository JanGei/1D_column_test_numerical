// This callback organizes the investigation of the data source and visibility of widgets

// Extracting data sources
var x           = source1.data['x'] 
var y           = source1.data['y']
var x2          = source2.data['x2']
var y2          = source2.data['y2']
var x3          = source3.data['xBTC']
var y3          = source3.data['yBTC']
var c_tot_array = sourcetot.data.c_tot_array

// Array to be displayed, initially filled with 0's
var c_tot     = Array(nX).fill(0)
// Extending Array in the time dimension
for (let i = 0; i < nT; i++) {
  c_tot[i] = new Array(nX).fill(0)
}

// Replacing zeroes with associated data from the source
for (let i = 0; i < nX; i++) {
  for (let j = 0; j < nT; j++) {
    c_tot[i][j] = c_tot_array[(j + i*nT)]
  }
}

// Column length [m]
const col_len   = col_len_sl.value; 

// Extracting status of RadioButtonGroup
var rg_SType = rg_ST.active
var rg_CPval = rg_CP.active

// Finding corresponding time and location for both plots
var xBTC  = x3[0];   
var tPV   = timestep_sl.value

// Fix point draw tool to x-axis and limit its range on x-axis
y3[0] = 0
if (x3[0]<=0.0001) {
  x3[0] = 0.001
} else if (x3[0] > col_len) {
  x3[0] = col_len
}

// Finding the closest value in spatial discretization to xBTC (x3)
var closestx = x.reduce(function(prev, curr) {
  return (Math.abs(curr - x3) < Math.abs(prev - x3) ? curr : prev);
});
// Finding corresponding index in x
var idx = x.indexOf(closestx)

// Assign new data for upper plot
for (let i = 0; i < x.length; i++){
  y[i] = c_tot[i][Math.round(tPV*100)]
}
// Assign new data for lower plot
for (let i = 0; i < x2.length; i++){
  y2[i] = c_tot[idx][i]
}

// Update Sliders
BTCp.title.text   = 'Breakthrough Curve at x = ' + xBTC.toFixed(3) + ' m (Drag diamond in upper plot to change)'

// Change UI for numerical model -- Computation is preformed in callback_compute_numerical.js

// Visibility conditions
if (rg_SType == 0) {
  rho_s_sl.visible = true
  Kd_sl.visible = true
  Kads_sl.visible = false
  s_max_sl.visible = false
  K_Fr_sl.visible = false
  Fr_n_sl.visible = false
} else if (rg_SType == 1) {
  rho_s_sl.visible = false
  Kd_sl.visible = false
  Kads_sl.visible = true
  s_max_sl.visible = true
  K_Fr_sl.visible = false
  Fr_n_sl.visible = false
} else if (rg_SType == 2) {
  rho_s_sl.visible = false
  Kd_sl.visible = false
  Kads_sl.visible = false
  s_max_sl.visible = false
  K_Fr_sl.visible = true
  Fr_n_sl.visible = true
}

if (rg_CPval==0) {
  pulse_inj_sl.visible = false
} else {
  pulse_inj_sl.visible = true
}

// Displaying units
var r_format  = r_dict[r_us.value]
var D_format  = D_dict[D_us.value]
var fl_format = fl_dict[fl_us.value]
reac_sl.format = r_format
disp_sl.format = D_format
flow_sl.format = fl_format

source1.change.emit();
source2.change.emit();