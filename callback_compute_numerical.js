function transport_num_CN(c_arr,Disp, sep_vel, dx_CN, dt_CN, nX, c_in, A_cn, b_cn, k){
  // Coefficients
  var p1 = Disp*dt_CN/dx_CN**2
  var p2 = sep_vel*dt_CN/(4*dx_CN)
  var p3 = sep_vel*dx_CN/Disp
  // First order reaction coefficient
  var p4 = k * dt_CN / 2

  for (let i = 0; i < nX; i++) {
    // Left hand side matrix A_CN
    if (i > 0 && i < nX-1) { // internal cells
      A_cn[i][i-1]  = -p1/2-p2
      A_cn[i][i]    = 1+p1+p4
      A_cn[i][i+1]  = -p1/2+p2
    } else if (i == 0) { // inflow cell 
      A_cn[i][i]    = p1*p3 + 2*p2*p3 + 1 + p1 +p4
      A_cn[i][i+1]  = -p1
    } else if (i == nX-1) { // outflow cell
      A_cn[i][i-1]  = -p1 - p1*p3 - 2*p2*p3
      A_cn[i][i]    = 1 + p1 + p1*p3 - 2*p2*p3 + p4
    }
    // Right hand side vector --> produces NaNs atm after 4 entries
    if (i > 0 && i < nX-1 && 1==1) { // internal cells
      b_cn[i]   = (p1/2+p2)*c_arr[i-1] + (1-p1-p4)*c_arr[i] + (p1/2-p2)*c_arr[i+1]
    } else if (i == 0) { // inflow cell
      b_cn[i]    = (-p1*p3 - 2*p2*p3 + 1 - p1 - p4)*c_arr[i] + p1*c_arr[i+1] + 2*(p1*p3 + 2*p2*p3)* c_in
    } else if (i == nX-1) { // outflow cell
      b_cn[i]    = (p1 + p1*p3 - 2*p2*p3)*c_arr[i-1] + (1 - p1 - p1*p3 + 2*p2*p3 - p4) *c_arr[i]
    }
  }
  const res = math.lusolve(A_cn,b_cn)
  return res
}

function total_conc(c,rg_SType,sorbed,rho_s,poros,K_Fr,Fr_n) {
  // This function sums up the concentration in the aqueous and solid phase per node
  var c_tot = Array(c.length).fill(0)
  if (rg_SType == 2) {
    for (let i = 0; i < c.length; i++) {
      c_tot[i] = (rho_s*(1-poros)*K_Fr*c[i]**(Fr_n-1)+ poros)*c[i]
    }
  } else {
    for (let i = 0; i < c.length; i++) {
      c_tot[i] = sorbed[i]*rho_s*(1-poros) + poros*c[i]
    }
  }
  return c_tot
}

// Extracting data sources
var x   = source1.data['x'] 
var y   = source1.data['y']
var ymin= source1.data['ymin']
var ymax= source1.data['ymax']
var x2  = source2.data['x2']
var y2  = source2.data['y2']
var x3  = source3.data['xBTC']
var y3  = source3.data['yBTC']

var rg_CP     = rg_CP.active                      // [0]
var rg_SType  = rg_ST.active                      // [0]

// Values needed for all models
const col_len   = col_len_sl.value;                 // [m]
const rad       = col_rad_sl.value;                 // [m]
const reac_l    = Math.exp(reac_sl.value[0])/3600;  // [1/s]
const reac_h    = Math.exp(reac_sl.value[1])/3600;  // [1/s]
const disp_l    = Math.exp(disp_sl.value[0])/3600;  // [m2/s]
const disp_h    = Math.exp(disp_sl.value[1])/3600;  // [m2/s]
const Q         = flow_sl.value/1000/1000/3600;     // [m3/s]
const poros     = poros_sl.value;                   // [-]
const t_inj     = pulse_inj_sl.value                // [s]
var xBTC        = x3[0];                            // [m]
const rho_s     = rho_s_sl.value                    // [kg/m3]
const Kd        = Kd_sl.value                       // [m3/kg]
const K_ads     = Kads_sl.value                     // [mol/m3]
const s_max     = s_max_sl.value                    // [mol/kg]
const K_Fr      = K_Fr_sl.value                     // [mmol^(1-n)*L^(n) / kg]
const Fr_n      = Fr_n_sl.value                     // [-]

// Derived entities
const A       = math.PI * rad**2;             // [m2]
const vel     = Q/A;                          // [m/s]
const sep_vel = vel / poros                   // [m/s]
const reac    = (reac_l + reac_h)/2           // [1/s] 
const Dis     = (disp_l + disp_h)/2           // [m2/s]  
const PS      = col_len * A * poros           // [m3]
const PV      = col_len/sep_vel               // [s] VEL oder SEP_VEL?
const c0      = 1;                            // [-] 
const nX      = x.length                      // [-]
const dx_CN   = col_len/nX                    // [m]
const dt_CN   = dx_CN / sep_vel               // [s]
const t_end   = PV * (x2[x2.length-1])        // [s]
const nT      = math.floor(t_end/dt_CN)       // [-]

// Initialize lists
var c_array     = Array(nX).fill(0)
var s_array     = Array(nX).fill(0)
var c_tot_array = Array(nT).fill(0)
var s_tot_array = Array(nT).fill(0)
// Creating multi-dimensional arrays
for (let i = 0; i < nT; i++) {
  c_tot_array[i] = new Array(nX).fill(0)
  s_tot_array[i] = new Array(nX).fill(0)
}
var A_CN        = Array(nX).fill(0)
var B_CN        = Array(nX).fill(0)
for (let i = 0; i < nX; i++) {
  A_CN[i]     = new Array(nX).fill(0)
}


console.log("Crank - Nicholson: \n"
              +"Spatial Discretization is " + dx_CN + " m\n"
              +"Temporal Discretization is " + dt_CN + " s\n" 
              +"Seepage velocity is "+ sep_vel + " m/s\n"
              +"Dispersion Coefficient is "+ Dis + " m2/s\n"
              +"There are " + nX + " spatial nodes\n"
              +"There are " + t_end/dt_CN + " temporal nodes\n"
              +"The Courant Number equals " + dt_CN * sep_vel /dx_CN
  )


for (let i = 0; i < 50; i++) {  //needs to be nT
  // Set inlet concentration 
  if (rg_CP == 0) {
      var c_in = c0
  } else if (rg_CP == 1 && i*dt_CN < t_inj) {
      var c_in = c0
  } else {
      var c_in = 0
  }

  // Transport
  // Transport with Crank Nicholson scheme
  var c_array = transport_num_CN(c_array,Dis,sep_vel,dx_CN,dt_CN,nX,c_in,A_CN,B_CN,reac)

  // Sorption
  if (1 == 1) {
  if (rg_SType == 0) { // Linear Sorption
      var c_tot_lin = total_conc(c_array,rg_SType,s_array,rho_s,poros)
      for (let j = 0; j < c_array.length; j++) {
        c_array[j][0] = c_tot_lin[j]/(Kd*(1-poros)*rho_s+poros)
        s_array[j]    = c_array[j][0]*Kd
      }

  } else if (rg_SType == 1) { // Langmuir Sorption
      var c_tot_lang = total_conc(c_array,rg_SType,s_array,rho_s,poros)

      for (let j = 0; j < c_array.length; j++) {
        var beta_lang  = (1-poros)*rho_s*s_max + poros*K_ads - c_tot_lang[j]
        var gamma_lang = -c_tot_lang[j]*K_ads

        c_array[j][0] = (-beta_lang + math.sqrt(beta_lang**2 - 4*poros*gamma_lang))/(2*poros)
        s_array[j]    = s_max*c_array[j]/(K_ads+c_array[j][0]) //[mol/kg]
      }

  } else if (rg_SType == 2) { // Freundlich Sorption
      var c_tot_Fr = total_conc(c_array,rg_SType,s_array,rho_s,poros,K_Fr,Fr_n)
      
      for (let j = 0; j < c_array.length; j++) {
        // Picard Iteration: Guess that everything is in the aqueous phase
        var c_old = c_tot_Fr[j]/poros
        // Loop until convergence criterion is met
        while (math.abs(c_old-c_array[j])>1e-9) {
          c_old = c_array[j][0]
          c_array[j][0] = c_tot_Fr[j]/(rho_s*(1-poros)*K_Fr*c_array[j][0]**(Fr_n-1)+poros)
        }
        s_array[j] = K_Fr * c_array[j]**Fr_n //[mmol/kg] --> different from other sorption types
      }
  }}

  // Store results
  for (let j = 0; j < c_array.length; j++) {
    c_tot_array[i][j] = c_array[j][0] 
    s_tot_array[i][j] = s_array[j]
  }
}

console.log(c_tot_array)

