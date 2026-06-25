
TITLE Voxel of Nitric Oxide Diffusion

COMMENT Equations from
        Pablo Fernandez Lopez, Patricio Garcia Baez, and Carmen Paz Suarez Araujo. Nitric Oxide Diﬀusion and Multi-compartmental
        Systems: Modeling and Implications DOI: 10.1007/978-3-319-26555-1 59

        Instituto Universitario de Ciencias y Tecnologias Cibernticas,
        Universidad de Las Palmas de Gran Canaria

        Adapted to be used in NEURON and NetPyNE by Scott McElroy
        SUNY Downstate scott.mcelroy@downstate.edu
ENDCOMMENT

NEURON {
    THREADSAFE
    POINT_PROCESS no_voxel
    RANGE conc, conc0, dx_pos, dx_neg, dy_pos, dy_neg, dz_pos, dz_neg, lam, F
    POINTER conc_xp, conc_xn, conc_yp, conc_yn, conc_zp, conc_zn
}

UNITS {
    (molar) = (1/liter)
    (nM) = (nanomolar)
}

PARAMETER {
    conc0  = 0 (nM)
    dx_pos = 0 (1/ms)
    dx_neg = 0 (1/ms)
    dy_pos = 0 (1/ms)
    dy_neg = 0 (1/ms)
    dz_pos = 0 (1/ms)
    dz_neg = 0 (1/ms)
    lam    = 0 (1/ms)
    F      = 0 (nM/ms)
}

ASSIGNED {
    conc_xp (nM)
    conc_xn (nM)
    conc_yp (nM)
    conc_yn (nM)
    conc_zp (nM)
    conc_zn (nM)
}

STATE {
    conc (nM)
}

INITIAL {
    conc = conc0
}

BREAKPOINT {
    SOLVE diffusion METHOD cnexp
}

DERIVATIVE diffusion {
    LOCAL dCdiff
    dCdiff = 0
    dCdiff = dCdiff + dx_pos*(conc_xp - conc)
    dCdiff = dCdiff + dx_neg*(conc_xn - conc)
    dCdiff = dCdiff + dy_pos*(conc_yp - conc)
    dCdiff = dCdiff + dy_neg*(conc_yn - conc)
    dCdiff = dCdiff + dz_pos*(conc_zp - conc)
    dCdiff = dCdiff + dz_neg*(conc_zn - conc)

    conc' = dCdiff - lam*conc + F
}
