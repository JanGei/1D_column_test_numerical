function transport(c_arr,c_in){
  // Move concentrations by one cell (Cr == 1)
  for (let i = 0; i < c_arr.length; i++){
    c_arr[c_arr.length-i-1] = c_arr[c_arr.length-i-2]
  }
  // Assign inlet concentration to first cell
  c_arr[0] = [c_in]

  return c_arr
}

function reactive_dispersion(c_arr,Disp, dx_CN, dt_CN, nX, A_cn, b_cn, k, rf){
  const dt_sub = dt_CN/rf
  // Coefficient for dispersion
  var p1 = Disp*dt_sub/dx_CN**2
  // Coefficient for first order reaction
  var p2 = k * dt_sub 

  for (let j = 0; j < rf; j++) {
    for (let i = 0; i < nX; i++) {
      // Left hand side matrix A_cn
      if (i > 0 && i < nX-1) { // internal cells
        A_cn[i][i-1]  = - p1/2
        A_cn[i][i]    = 1 + p1
        A_cn[i][i+1]  = - p1/2
      } else if (i == 0) { // inflow cell 
        A_cn[i][i]    = 1 + p1/2 
        A_cn[i][i+1]  = - p1/2
      } else if (i == nX-1) { // outflow cell
        A_cn[i][i-1]  = - p1/2 
        A_cn[i][i]    = 1 + p1/2 
      }
      // Right hand side vector 
      if (i > 0 && i < nX-1) { // internal cells
        b_cn[i]   = p1/2*c_arr[i-1] + (1-p1)* c_arr[i] + p1/2*c_arr[i+1]
      } else if (i == 0) { // inflow cell
        b_cn[i]    = (1-p1/2)*c_arr[i] + p1/2*c_arr[i+1]
      } else if (i == nX-1) { // outflow cell
        b_cn[i]    = (1-p1/2)*c_arr[i] + p1/2*c_arr[i-1]
      }
    }
    c_arr = math.lusolve(A_cn,b_cn) //this gives us an array with an array for each entry
  }
  return c_arr
}

function reactive_dispersion_fully_implicit(c_arr,Disp, dx_CN, dt_CN, nX,  A_cn, k){
  // Coefficient for dispersion
  var p1 = Disp*dt_CN/dx_CN**2
  // Coefficient for first order reaction
  var p2 = k * dt_CN 
  
    for (let i = 0; i < nX; i++) {
      // Left hand side matrix A_cn
      if (i > 0 && i < nX-1) { // internal cells
        A_cn[i][i-1]  = - p1
        A_cn[i][i]    = 1 + 2*p1 +p2
        A_cn[i][i+1]  = - p1
      } else if (i == 0) { // inflow cell 
        A_cn[i][i]    = 1 + p1 + p2
        A_cn[i][i+1]  = - p1 
      } else if (i == nX-1) { // outflow cell
        A_cn[i][i-1]  = - p1 
        A_cn[i][i]    = 1 + p1 + p2
      }
    }
    c_arr = math.lusolve(A_cn,c_arr)
  
  return c_arr
}

function total_conc(c,rg_SType,sorbed,rho_s,poros,K_Fr,Fr_n) {
  // This function sums up the concentration in the aqueous and solid phase per node
  var c_tot_intermed = Array(c.length).fill(0)
  if (rg_SType == 2) {
    for (let i = 0; i < c.length; i++) {
      c_tot_intermed[i] = (rho_s*(1-poros)*K_Fr*c[i]**(Fr_n-1)+ poros)*c[i]
    }
  } else {
    for (let i = 0; i < c.length; i++) {
      c_tot_intermed[i] = sorbed[i]*rho_s*(1-poros) + poros*c[i]
    }
  }
  return c_tot_intermed
}

var rg_CP     = rg_CP.active                      // [0]
var rg_SType  = rg_ST.active                      // [0]

// Values needed for all models
const col_len   = col_len_sl.value;                 // [m]
const rad       = col_rad_sl.value;                 // [m]
const reac      = Math.exp(reac_sl.value)/3600;     // [1/s]
const Dis       = Math.exp(disp_sl.value)/3600;     // [m2/s]
const Q         = flow_sl.value/1000/1000/3600;     // [m3/s]
const poros     = poros_sl.value;                   // [-]
const t_inj     = pulse_inj_sl.value                // [s]
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
const PS      = col_len * A * poros           // [m3]
const PV      = col_len/sep_vel               // [s] VEL oder SEP_VEL?
const c0      = 1;                            // [-] 
const dx_CN   = col_len/nX                    // [m]
const dt_CN   = dx_CN / sep_vel               // [s]
const Ne      = 4* Dis / sep_vel / dx_CN 	    // [-]
const rf      = math.ceil(Ne)
const t_end   = PV * PVspan        // [s]
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
var A_cn        = Array(nX).fill(0)
var B_CN        = Array(nX).fill(0)
for (let i = 0; i < nX; i++) {
  A_cn[i]     = new Array(nX).fill(0)
}


console.log("Crank - Nicholson: \n"
              +"Spatial Discretization is " + dx_CN + " m\n"
              +"Temporal Discretization is " + dt_CN + " s\n" 
              +"Seepage velocity is "+ sep_vel + " m/s\n"
              +"Dispersion Coefficient is "+ Dis + " m2/s\n"
              +"There are " + nX + " spatial nodes\n"
              +"There are " + t_end/dt_CN + " temporal nodes\n"
              +"The Neumann Number is equal to " + Ne/4 +"\n"
              +"The refinement factor is " + rf +"\n"
              +"The Courant Number equals " + dt_CN * sep_vel /dx_CN 
  )


for (let i = 0; i < nT; i++) {  //needs to be nT
  // Set inlet concentration 
  if (rg_CP == 0) {
      var c_in = c0
  } else if (rg_CP == 1 && i*dt_CN < t_inj) {
      var c_in = c0
  } else {
      var c_in = 0
  }

  // Transport
  var c_array = transport(c_array,c_in)

  // Reactive Dispersion with Crank Nicholson scheme rf times
  if (1==1) {
  //var c_array = reactive_dispersion(c_array,Dis,dx_CN,dt_CN,nX,A_cn,B_CN,reac,rf)
  // Fully implicit calculation does not need sub-timestepping
  var c_array = reactive_dispersion_fully_implicit(c_array,Dis,dx_CN,dt_CN,nX,A_cn,reac)
  }
  console.log((i+1) + " out of " + nT + " computational steps have been perfromed")

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
    c_tot_array[j][i] = c_array[j][0] 
    s_tot_array[j][i] = s_array[j]
  }
}

console.log(c_tot_array)

// entering results in big source
// extracting big data source (next 12 lines are not really needed)
var sourcetot_array = sourcetot.data.c_tot_array

for (let i = 0; i < nX; i++) {
  for (let j = 0; j < nT; j++) {
    sourcetot_array[(j + i*nT)] = c_tot_array[i][j] 
  }
}

timestep_sl.title = 'Pore Volume (1PV =' + (PV/3600).toFixed(2) +'h)';


sourcetot.change.emit();