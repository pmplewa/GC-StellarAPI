from astropy import units, constants
import numpy as np
from scipy.optimize import newton


R0 = 8.33e3 # distance to the Galactic Center
R0_unit = units.pc

M0 = 4.31e6 # mass of the central black hole
M0_unit = units.solMass

m_unit = M0 * M0_unit # mass unit
a_unit = units.arcsec # angle unit
t_unit = units.yr # time unit

l_unit = a_unit.to(units.rad) * (R0 * R0_unit) # length unit

G = float(constants.G.cgs * m_unit * t_unit**2 / l_unit**3)


class Particle():
    """A particle with a (3D) position and velocity.

    :param float x: x-position
    :param float y: y-position
    :param float z: z-position
    :param float vx: x-velocity
    :param float vy: y-velocity
    :param float vz: z-velocity ("radial velocity")
    :param float m=0: mass of the particle
    """
    def __init__(self, x, y, z, vx, vy, vz, mass=0):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.mass = mass

class Orbit():
    """A Kepler orbit represented by its orbital elements.

    :param float a: semi-major axis
    :param float e: eccentricity
    :param float inc: inclination
    :param float Omega: longitude of ascending node
    :param float omega: argument of pericenter
    :param float tp: time of pericenter
    :param float m=0: mass of the orbiting body
    """
    def __init__(self, a, e, inc, Omega, omega, tp, mass=0):
        self.a = a
        self.e = e
        self.inc = inc
        self.Omega = Omega
        self.omega = omega
        self.tp = tp
        self.mass = mass

class ProperMotion():
    """A proper motion up to second order.

    :param float t0: reference epoch
    :param float x0: x-position at t0
    :param float y0: x-position at t0
    :param float vx0: x-velocity at t0
    :param float vy0: y-velocity at t0
    :param float ax0=None: x-acceleration at t0
    :param float ay0=None: y-acceleration at t0
    """
    def __init__(self, t0, x0, y0, vx0, vy0, ax0=None, ay0=None):
        self.t0 = t0
        self.x0 = x0
        self.y0 = y0
        self.vx0 = vx0
        self.vy0 = vy0
        self.ax0 = ax0
        self.ay0 = ay0

class Body():
    """A body, typically a star, with an orbit or proper motion measurement.

    :param string name: unique id of the body
    :param Particle primary: reference particle
    :param float m=0: mass of the body
    :param float magnitude=None: brightness magnitude
    :param string spec_type=None: stellar spectral type

    :var Orbit orbit: orbit (if at least a, e, inc, Omega, omega and tp in
        kwargs)
    :var ProperMotion promo: proper motion (if at least x0, y0, vx0, vy0 and
        t0 in kwargs)
    """
    def __init__(self, name, primary, mass=0, magnitude=None, spec_type=None,
                 names=[], orbit=None, promo=None, **kwargs):
        self.id = name
        self.primary = primary
        self.mass = mass
        self.magnitude = magnitude
        self.spec_type = spec_type
        self.names = names
        if orbit is not None:
            self.orbit = Orbit(mass=mass, **orbit)
        elif promo is not None:
            self.promo = ProperMotion(**promo)

    def __repr__(self):
      return f"<{self.__class__.__name__} {self.id}>"

    def mean_motion(self):
        """Return the mean motion of the body's orbit.

        :rtype: float
        :raises Exception: if the star has no measured orbit
        """
        if self.orbit is not None:
            mu = G * (self.mass + self.primary.mass)
            return mean_motion(mu, self.orbit.a)
        else:
            raise Exception(f"Body {self.id} has no measured orbit.")

    def period(self):
        """Return the orbital period of the body's orbit.

        :rtype: float
        :raises Exception: if the star has no measured orbit
        """
        if self.orbit is not None:
            mu = G * (self.mass + self.primary.mass)
            return period(mu, self.orbit.a)
        else:
            raise Exception(f"Body {self.id} has no measured orbit.")

    def locate(self, t):
        """Return the position and velocity of the body at a specific time,
        based on a a measurement of its orbit or its proper motion.

        :param float t: time of observation

        :rtype: Particle
        :raises Exception: if the body can not be located
        """
        if self.orbit is not None:
            return orbit_to_particle(self.orbit, self.primary, t)
        elif self.promo is not None:
            return promo_to_particle(self.promo, self.primary, t)
        else:
            raise Exception(f"Body {self.id} can not be located.")

def mod2pi(x):
    """Wrap an angle to the range [-pi, pi].

    :param float x: input angle

    :rtype: float
    """
    return (x + np.pi) % (2 * np.pi) - np.pi

def arccos2(x, sign):
    """Returns arccos(x) choosing the correct quadrant.

    :param float x: input value

    :rtype: float
    """
    return np.arccos(x) if (sign >= 0) else 2 * np.pi - np.arccos(x)

def mean_motion(mu, a):
    """Return the mean motion of a Kepler orbit.

    :param float mu: gravitational parameter G*(M+m)
    :param float a: semi-major axis

    :rtype: float
    """
    return np.sign(a) * np.sqrt(np.fabs(mu / a**3))

def period(mu, a):
    """Return the period of a Kepler orbit.

    :param float mu: gravitational parameter G*(M+m)
    :param float a: semi-major axis

    :rtype: float
    """
    n = mean_motion(mu, a)
    return 2*np.pi/n

def eccentric_anomaly(e, M, *args, **kwargs):
    """Convert a mean anomaly into an eccentric anomaly.

    :param float e: eccentricity
    :param float M: mean anomaly

    :rtype: float
    """
    if e < 1:
        f = lambda E: E - e * np.sin(E) - M
        fp = lambda E: 1 - e * np.cos(E)
        E0 = M if e < 0.8 else np.sign(M) * np.pi
        E = mod2pi(newton(f, E0, fp, *args, **kwargs))
    else:
        f = lambda E: E - e * np.sinh(E)-M
        fp = lambda E: 1 - e * np.cosh(E)
        E0 = np.sign(M) * np.log(2 * np.fabs(M) / e + 1.8)
        E = newton(f, E0, fp, *args, **kwargs)
    return E

def true_anomaly(e, M):
    """Convert a mean anomaly into a true anomaly.

    :param float e: eccentricity
    :param float M: mean anomaly

    :rtype: float
    """
    E = eccentric_anomaly(e, M)
    if e > 1:
        return 2 * np.arctan(np.sqrt((1 + e) / (e - 1)) * np.tanh(E/2))
    else:
        return 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))

def particle_to_orbit(particle, primary, t, tol=1e-8):
    """Find the orbital elements corresponding to a set of positions
    and velocities.

    :param Particle particle: particle instance to convert
    :param Particle primary: reference particle
    :param float t: time of observation

    :rtype: Orbit
    """
    assert primary.mass > 0, "Primary particle has no mass."

    mu = G * (particle.mass + primary.mass)

    dx = particle.x - primary.x
    dy = particle.y - primary.y
    dz = particle.z - primary.z
    dvx = particle.vx - primary.vx
    dvy = particle.vy - primary.vy
    dvz = particle.vz - primary.vz
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    assert r > 0, "Particle is on top of primary particle."

    v2 = dvx**2 + dvy**2 + dvz**2
    vc2 = mu / r
    a = -mu / (v2 - 2 * vc2)

    hx = (dy * dvz - dz * dvy)
    hy = (dz * dvx - dx * dvz)
    hz = (dx * dvy - dy * dvx)
    h = np.sqrt(hx**2 + hy**2 + hz**2)

    dv2 = v2 - vc2
    vr = (dx * dvx + dy * dvy + dz * dvz) / r
    rvr = r*vr

    ex = (dv2 * dx - rvr * dvx) / mu
    ey = (dv2 * dy - rvr * dvy) / mu
    ez = (dv2 * dz - rvr * dvz) / mu
    e = np.sqrt(ex**2 + ey**2 + ez**2)

    inc = np.arccos(hz / h)

    nx = -hy
    ny = hx
    n = np.sqrt(nx**2 + ny**2)

    Omega = arccos2(nx / n, ny)

    if (e < 1):
        E = arccos2((1 - r / a) / e, vr)
        M = E - e * np.sin(E)
    else:
        E = np.sign(vr) * np.arccosh((1 - r / a) / e)
        M = e * np.sinh(E) - E

    tp = t - M / np.fabs(mean_motion(mu, a))

    if (inc < tol) or (inc > np.pi-tol):
        pomega = arccos2(ex/e, ey)
        if (inc < np.pi / 2):
            omega = pomega-Omega
        else:
            omega = Omega - pomega
    else:
        omega = arccos2((nx * ex + ny * ey) / (n * e), ez)

    return Orbit(a, e, inc, Omega, omega, tp, particle.mass)

def orbit_to_particle(orbit, primary, t):
    """Predict the position an velocity of a body from its orbital elements.

    :param Orbit orbit: orbit instance to convert
    :param Particle primary: reference particle
    :param float t: time of observation

    :rtype: Particle
    """
    assert orbit.e != 1, "Can not initialize a radial orbit (e = 1)."
    assert orbit.e >= 0, "A valid orbit must have e >= 0."
    if orbit.e > 1:
        assert orbit.a < 0, "A bound orbit (a > 0) must have e < 1."
    else:
        assert orbit.a > 0, "An unbound orbit (a < 0) must have e > 1."

    mu = G * (orbit.mass + primary.mass)

    n = mean_motion(mu, orbit.a)
    M = mod2pi(n * (t - orbit.tp))
    f = true_anomaly(orbit.e, M)
    assert orbit.e * np.cos(f) > -1, "An unbound orbit can not have f set beyond \
        the range allowed by the parabolic asymptotes."

    r = orbit.a * (1 - orbit.e**2) / (1 + orbit.e * np.cos(f))
    v = np.sqrt(mu / orbit.a / (1 - orbit.e**2))

    cO = np.cos(orbit.Omega)
    sO = np.sin(orbit.Omega)
    co = np.cos(orbit.omega)
    so = np.sin(orbit.omega)
    cf = np.cos(f)
    sf = np.sin(f)
    ci = np.cos(orbit.inc)
    si = np.sin(orbit.inc)

    x = primary.x + r * (cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    y = primary.y + r * (sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    z = primary.z + r * (so*cf + co*sf)*si

    vx = primary.vx + v * ((orbit.e+cf)*(-ci*co*sO-cO*so) - sf*(co*cO-ci*so*sO))
    vy = primary.vy + v * ((orbit.e+cf)*(ci*co*cO-sO*so) - sf*(co*sO+ci*so*cO))
    vz = primary.vz + v * ((orbit.e+cf)*co*si - sf*si*so)

    return Particle(x, y, z, vx, vy, vz, orbit.mass)

def promo_to_particle(promo, primary, t):
    """Predict the position and velocity of a body from its proper motion.

    :param ProperMotion promo: proper motion instance to convert
    :param Particle primary: reference particle
    :param float t: time of observation

    :rtype: Particle
    """
    vx = primary.vx + promo.vx0
    vy = primary.vy + promo.vy0

    x = primary.x + promo.x0 + vx * (t - promo.t0)
    y = primary.y + promo.y0 + vy * (t - promo.t0)

    if promo.ax0:
        x += promo.ax0 / 2 * (t - promo.t0)**2
        vx += promo.ax0 * (t - promo.t0)
    if promo.ay0:
        y += promo.ay0 / 2 * (t - promo.t0)**2
        vy += promo.ax0 * (t - promo.t0)

    return Particle(x, y, None, vx, vy, None)
