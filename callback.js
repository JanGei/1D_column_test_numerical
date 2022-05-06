// Extracting data sources
var x   = source1.data['x'] 
var y   = source1.data['y']
var x2  = source2.data['x2']
var y2  = source2.data['y2']
var x3  = source3.data['xBTC']
var y3  = source3.data['yBTC']

const col_len   = col_len_sl.value;                 // [m]

var rg_SType = rg_ST.active
// array with the entire data
//const c_tot = c_tot
// finding corresponding time and location for both plots
var xBTC  = x3[0];   
var tPV   = timestep_sl.value

// Fix point draw tool to x-axis and limit its range on x-axis
y3[0] = 0
if (x3[0]<=0.0001) {
  x3[0] = 0.001
} else if (x3[0] > col_len) {
  x3[0] = col_len
}

// Finding the closest value in x to xBTC (x3)
var closestx = x.reduce(function(prev, curr) {
  return (Math.abs(curr - x3) < Math.abs(prev - x3) ? curr : prev);
});
// Finding corresponding index in x
var idx = x.indexOf(closestx)

// assign new data for upper plot
for (let i = 0; i < x.length; i++){
  y[i] = c_tot[i][Math.round(tPV*100)]
}
// assign new data for lower plot
for (let i = 0; i < x2.length; i++){
  y2[i] = c_tot[idx][i]
}

// Update Sliders
BTCp.title.text   = 'Breakthrough Curve at x = ' + xBTC.toFixed(3) + ' m (Drag diamond in upper plot to change)'


// Change UI for numerical model -- Computation is preformed in different fiel
// Make sorption type options visible

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


if (rg_CP==0) {
  pulse_inj_sl.visible = false
} else {
  pulse_inj_sl.visible = true
}

source1.change.emit();
source2.change.emit();