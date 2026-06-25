 NEURON {
	POINT_PROCESS GABA
	NONSPECIFIC_CURRENT i
	RANGE g, i, e, tauRise, tauFall, gmax
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
    
	tauRise = 0.9 (ms)
	tauFall = 26.5 (ms)
	e = -80 (mV)
        gmax (microsiemens)
	}

ASSIGNED {
    
	v (mV)
	i (nA)
	g (microsiemens)
	factor
        
}

STATE {
	A (microsiemens)
	B (microsiemens)
}

INITIAL {
	LOCAL tp
	if (tauRise/tauFall > .9999) {
		tauRise = .9999*tauFall
	}
	A = 0
	B = 0
	tp = 3.2 :(tauRise*tauFall)/(tauFall - tauRise) * log(tauFall/tauRise)
	factor = -exp(-tp/tauRise) + exp(-tp/tauFall)
	factor = 1/factor

    
}

BREAKPOINT {
    : Here the conductance is updated each time step, while the NET_RECEIVE block
    : is only invoked by being contacted by a NetCon object.
  	SOLVE state METHOD cnexp
   	g = gmax*(B - A)
   	i = g*(v - e)   	
}

DERIVATIVE state {
	A' = -A/tauRise
	B' = -B/tauFall
}

NET_RECEIVE(weight (microsiemens)) {
      A = A + weight*factor
      B = B + weight*factor
}
