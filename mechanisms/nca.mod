TITLE nca.mod  
 
COMMENT
konduktivitas valtozas hatasa- somaban 
ENDCOMMENT
 
UNITS {
    (mA) =(milliamp)
    (mV) =(millivolt)
}
 
? interface 
NEURON { 
	SUFFIX nca
	USEION nca READ enca WRITE inca
	RANGE  gnca
	RANGE gncabar
	RANGE cinf, ctau, dinf, dtau, inca
}
 
PARAMETER {
	v (mV)
	celsius = 6.3 (degC)
	gncabar (mho/cm2)
}
 
STATE {
	c d
}
 
ASSIGNED {
	gnca (mho/cm2)
	inca (mA/cm2)
	enca (mV)

	cinf dinf
	ctau (ms) dtau (ms)     
} 

? currents
BREAKPOINT {
	SOLVE states METHOD cnexp
    gnca = gncabar*c*c*d
	inca = gnca*(v-enca)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	c = cinf
	d = dinf
}

? states
DERIVATIVE states {	:Computes state variables m, h, and n 
    rates(v)	:      at the current v and dt.
	c' = (cinf-c)/ctau
	d' = (dinf-d)/dtau
}

? rates
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
    LOCAL  alpha, beta, sum
    :"c" NCa activation system
	
    alpha = -0.19*vtrap(v-19.88,-10)
	beta = 0.046*exp(-v/20.73)
	sum = alpha+beta        
	ctau = 1/sum      cinf = alpha/sum
                :"d" NCa inactivation system
	alpha = 0.00016/exp(-v/48.4)
	beta = 1/(exp((-v+39)/10)+1)
	sum = alpha+beta        
	dtau = 1/sum      dinf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

