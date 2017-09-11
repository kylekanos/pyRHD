# simply RHD solver

from math import sqrt

# full vars
GAMMA = 5.0/3.0
GAMM1 = GAMMA - 1.0
LENGTH = 1.0
NCELLS = 400
DX = 1.0 / NCELLS

class StateObj(object):
    def __init__(self, nvars):
        if nvars == 5:
            # primitive vars
            self.p   = 0.0
            self.v   = 0.0
            self.rho = 0.0
            self.gam = 0.0
            self.eps = 0.0
        elif nvars == 3:
            # conservative vars
            self.mass = 0.0
            self.mom  = 0.0
            self.erg  = 0.0
        else:
            raise ValueError("cannot create a start vector with %d states!" % nvars)

#!> get the position
def getX(i):
    return LENGTH * (i + 0.5) / NCELLS

#!> get the wavespeed
def getWavespeed(P):
    return sqrt(GAMM1 * GAMMA * P.eps / (1.0 + GAMMA * P.eps))

#!> get the lorentz factor
def getLorentz(vel):
    return 1.0 / sqrt(1.0 - vel**2)

#!> get the conservative vars from the primitive
def getConservative(P):
    u = StateObj(3)
    h = 1.0 + P.eps + P.p / P.rho
    u.mass = P.rho * P.gam
    u.mom  = P.rho * h * P.gam**2 * P.v
    u.erg  = P.rho * h * P.gam**2 - P.p - u.mass
    return u

#!> get the primitive vars from the conservative
def getPrimitive(U, P):
    for i in range(NCELLS):
        delta = 10.0
        P[i].p = 10.0
        # do newton raphson iteration to get pressure
        while abs(delta) > 1e-13:
            if P[i].p < (abs(U[i].mom) - U[i].erg - U[i].mass):
                P[i].p = abs(U[i].mom) - U[i].erg - U[i].mass
            P[i].v   = U[i].mom / (U[i].erg + U[i].mass + P[i].p)
            P[i].gam = getLorentz(P[i].v)
            P[i].rho = U[i].mass / P[i].gam
            P[i].eps = (U[i].erg + U[i].mass * (1.0 - P[i].gam) + P[i].p * (1.0 - pow(P[i].gam,2))) / (U[i].mass * P[i].gam)

            # trying to match P_{eos} with P_{guess}
            f = GAMM1 * P[i].rho * P[i].eps - P[i].p
            cs = getWavespeed(P[i])
            dfdx = pow(P[i].v * cs, 2) - 1.0
            delta = -f / dfdx
            P[i].p += delta
        # update rest of primitives
        P[i].v   = U[i].mom / (U[i].erg + U[i].mass + P[i].p)
        P[i].gam = getLorentz(P[i].v)
        P[i].rho = U[i].mass / P[i].gam
        P[i].eps = (U[i].erg + U[i].mass * (1.0 - P[i].gam) + P[i].p*(1.0 + pow(P[i].gam,2)))/(U[i].mass * P[i].gam)


#!> initialize the system
def initialize(P, U):
    for i in range(NCELLS):
        oneByZ = 1.0 / (i + 1)
        P[i].rho = 1.00 * oneByZ
        P[i].v   = 0.00
        P[i].p   = 0.17 * oneByZ
        P[i].gam = getLorentz(P[i].v)
        P[i].eps = P[i].p/(GAMM1 * P[i].rho)

        U[i] = getConservative(P[i])

    with open("init.dat", "w") as f:
        for i in range(NCELLS):
            f.write("%0.5f\t%0.15f\t%0.15f\t%0.15f\n" % (getX(i), P[i].rho, P[i].v, P[i].p))

def updateFlux(U, P, F):
    for i in range(NCELLS):
        F[i].mass = U[i].mass * P[i].v
        F[i].mom  = U[i].mom * P[i].v + P[i].p
        F[i].erg  = U[i].mom - U[i].mass * P[i].v

def updateBCs(U, P):
    # fixed bc's
    P[0].rho = 1.00
    P[0].p   = 0.17
    P[0].v   = 0.0
    P[0].gam = getLorentz(P[0].v)
    P[0].eps = P[0].p / (GAMM1 * P[0].rho)
    U[0] = getConservative(P[0])
    # outlow
    P[-1].rho = 2.0 * P[-2].rho - P[-3].rho
    P[-1].v   = 2.0 * P[-2].v   - P[-3].v
    P[-1].p   = 2.0 * P[-2].p   - P[-3].p
    P[-1].gam = getLorentz(P[-1].v)
    P[-1].eps = P[-1].p / (GAMM1 * P[-1].rho)
    #P[-1].gam = 2.0 * P[-2].gam - P[-3].gam
    #P[-1].eps = 2.0 * P[-2].eps - P[-3].eps
    U[-1] = getConservative(P[-1])

def getFluxSpeed(Pl, Pu, up):
    vel = 0.5 * (Pl.v + Pu.v)
    cs  = 0.5 * (getWavespeed(Pl) + getWavespeed(Pu))
    if up == 1:
        return max(0.0, ((vel + cs) / (1.0 + vel * cs)))
    else:
        return max(0.0, ((vel - cs) / (1.0 - vel * cs)))

def getHLL(i, U, P, F):
    ap = getFluxSpeed(P[i], P[i+1], 1)
    an = getFluxSpeed(P[i], P[i+1], 0)

    den = 1.0 / (ap + an)

    Fhll = StateObj(3)
    Fhll.mass = (ap * F[i].mass - an * F[i+1].mass + ap * an * (U[i+1].mass - U[i].mass)) * den
    Fhll.mom  = (ap * F[i].mom  - an * F[i+1].mom  + ap * an * (U[i+1].mom  - U[i].mom )) * den
    Fhll.erg  = (ap * F[i].erg  - an * F[i+1].erg  + ap * an * (U[i+1].erg  - U[i].erg )) * den

    return Fhll

#!> step in time-step
def timeStep(dt, U, F, P):
    dtdx = dt / DX
    H = [StateObj(3) for i in range(NCELLS)]
    # get HLL flux
    for i in range(NCELLS-1):
        H[i] = getHLL(i, U, P, F)
    # update conservative
    for i in range(NCELLS-1):
        U[i].mass -= dtdx * H[i].mass
        U[i].mom  -= dtdx * H[i].mom
        U[i].erg  -= dtdx * H[i].erg
        U[i+1].mass += dtdx * H[i].mass
        U[i+1].mom  += dtdx * H[i].mom
        U[i+1].erg  += dtdx * H[i].erg
    # update primitives
    getPrimitive(U, P)
    updateBCs(U, P)
    updateFlux(U, P, F)

# state variables
prim = [StateObj(5) for i in range(NCELLS)]
cons = [StateObj(3) for i in range(NCELLS)]
flux  = [StateObj(3) for i in range(NCELLS)]
#!> run main program
def run():
    tend = 0.04
    dt = tend / 250
    initialize(prim, cons)
    updateFlux(cons, prim, flux)
    t = 0.0
    for i in range(100):
        timeStep(dt, cons, flux, prim)
        t += dt

    with open('soln.txt', 'w') as f:
        for i in range(NCELLS):
            f.write("%0.5f\t%0.15f\t%0.15f\t%0.15f\n" % (getX(i), prim[i].rho, prim[i].v, prim[i].p))
    print('done')

if __name__ == '__main__':
    run()
