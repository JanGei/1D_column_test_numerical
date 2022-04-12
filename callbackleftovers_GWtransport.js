function transport_num(c,c_in,disp,dx,dt){
    // This function approximates advection and dispersion numerically by a FVM
    
    // Move concentration by 1 cell
    for (let i = 1; i < c.length; i++) {
      c[c.length-i] = c[c.length-i-1]
    }
    // First cell gets inlet concentration
    c[0] = c_in
    var Jd = []
    // Dispersive fluxes between cells
    for (let i = 0; i < c.length+1; i++) {
      if (i == 0) {
        Jd[i] = 0
      } else if (i == c.length) {
        Jd[i] = Jd[i-1]
      } else {
        Jd[i] = (c[i-1] - c[i])/dx*disp
      }
    }
  
    for (let i = 0; i < Jd.length-1; i++) {
      c[i] = c[i] + dt/dx * (Jd[i] - Jd[i+1])
    }
  
    return c
}

// This discretization follows the scheme of the GW transport code
// Requirement 1: Neumann-Number = 1/3
// Requirement 2: Courant-Number = 1
const dx    = 3*Dis/sep_vel
const dt    = 3*Dis/sep_vel**2

// Lists for spatial and temporal span
var x_num = []
var t_num = []

for (let i = 0; ((i+1)*dx) < col_len; i++) {
    x_num[i] = (i+1)*dx
}
for (let i = 0; ((i+1)*dt) < t_end; i++) {
    t_num[i] = (i+1)*dt
}
console.log("GW_Transport: \n"
              +"Spatial Discretization is " + dx + " m\n"
              +"Temporal Discretization is " + dt + " s\n" 
              +"Seepage velocity is "+ sep_vel + "m/s\n"
              +"Dispersion Coefficient is "+ Dis + "m2/s\n"
              +"There are " + x_num.length + " spatial nodes\n"
              +"There are " + t_num.length + " temporal nodes\n"
              +"The Courant Number equals " + dt * sep_vel /dx
)