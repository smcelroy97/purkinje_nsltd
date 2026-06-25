NEURON {
POINT_PROCESS Ampa_PC
RANGE tau0, tau1, gmax, e, i, g
NONSPECIFIC_CURRENT i
} 

UNITS {
(nA) = (nanoamp)
(mV) = (millivolt)
(nS) = (nanosiemens) 
} 

PARAMETER {
: onset = 0 (ms)
tau0 = 0.5 (ms)
tau1 = 1.2 (ms)
gmax = 0.7 (nS) <0.0002, 0.05>
e = 0 (mV)
}

ASSIGNED {
v (mV)
i (nA)
g (nS)
factor
: total
}

STATE {
   A :(nanosiemens) 
   B :(nanosiemens) 
}

INITIAL {LOCAL tpeak
if(tau0/tau1>0.999) {
tau0 = 0.999*tau1
}
A=0 :(nanosiemens) 
B=0 : (nanosiemens) 
UNITSOFF
tpeak=0.8 : (tau0*tau1)/(tau1-tau0)*log(tau1/tau0)
factor=-exp(-tpeak/tau0)+exp(-tpeak/tau1)
factor = 1/factor
}
UNITSON
BREAKPOINT {
   SOLVE state METHOD cnexp
   g = gmax*(B - A)
   i = (0.001)*g*(v - e) 
 }
UNITSOFF
DERIVATIVE state {
   A' = -A/tau0
   B' = -B/tau1
}


NET_RECEIVE (weight (nanosiemens) ) {
A = A + weight*factor
B = B + weight*factor
}
UNITSON



