COMMENT

Sodium current for the axon

Conductances taken from Traub & Miles 1995 paper where ratio of Na for 
soma:axon is 1:5

References:

1.	Martina, M., Vida, I., and Jonas, P.  Distal initiation and active
	propagation of action potentials in interneuron dendrites,
	Science, 287:295-300, 2000.

			soma	axon-lacking dend	axon-bearing dend
Na+	gmax	    107 ps/um2	   117 ps/um2		   107 ps/um2
	slope 	    10.9 mV/e	   11.2 mV/e		   11.2 mV/e
	V1/2        -37.8 mV       -45.6 mV                -45.6 mV



2.	Marina, M. and Jonas, P.  Functional differences in Na+ channel
	gating between fast-spiking interneurons and principal neurones of rat
	hippocampus, J. Physiol., 505.3:593-603, 1997.

*Note* The interneurons here are basket cells from the dentate gyrus.

Na+	Activation V1/2				-25.1 mV
	slope			 		11.5
	Activation t (-20 mV)	 		0.16 ms
	Deactivation t (-40 mV)	 		0.13 ms
 	Inactivation V1/2			-58.3 mV
	slope			 		6.7
	onset of inactivation t (-20 mV)	1.34 ms
	onset of inactivation t (-55 mV)	18.6 ms
	recovery from inactivation t		2.0 ms
	(30 ms conditioning pulse)
	recovery from inactivation t		2.7 ms
	(300 ms conditioning pulse)

ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
    SUFFIX Naaxon
    USEION na READ ena WRITE ina
    NONSPECIFIC_CURRENT il
    RANGE gnaaxon, gl, el, ina
	GLOBAL minf, hinf, mtau, htau
}
 
PARAMETER {
	v (mV)
	celsius = 24 (degC)
    gnaaxon = .0107  (mho/cm2)
    ena     = 90     (mV)
    gl      = .00005 (mho/cm2)
    el      = -70    (mV)
}
 
STATE { m h }
 
ASSIGNED {
    ina (mA/cm2)
    il (mA/cm2)
    minf 
	hinf
	mtau (ms)
	htau (ms)
}
 
INITIAL {
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnaaxon*pow(minf, 3)*h*(v - ena)    
	il = gl*(v - el)
}

DERIVATIVE states {	:exact when v held constant
	rates(v)
	h' = (hinf - h) / htau 
}
UNITSOFF

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at 
	:current v.
	:Call once from HOC to initialize inf at resting v.
	LOCAL alpha, beta
	
	alpha = 0.1*vtrap(-(v+38),10)
	beta  = 4*exp(-(v+63)/18)
	mtau  = 1/(alpha + beta)
	minf  = alpha*mtau
	alpha = 0.07*Exp(-(v+63)/20)
	beta  = 1/(1+Exp(-(v+33)/10))
	htau  = 1/(alpha + beta)
	hinf  = alpha*htau
}

FUNCTION vtrap(x,y) {
	:Traps for 0 in denominator of rate eqns.
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(Exp(x/y) - 1)
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}
UNITSON
