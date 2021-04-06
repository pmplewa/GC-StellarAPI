import numpy as np
import os

from astropy import units, constants
from astropy.units import Unit
from flask import jsonify, request, render_template, send_from_directory
from werkzeug.exceptions import HTTPException
from werkzeug.routing import NumberConverter

from . import app
from .data import stars, black_hole
from .orbit_tools import (R0, M0, R0_unit, M0_unit, G, Orbit, Particle,
    orbit_to_particle, particle_to_orbit)


class InvalidUsage(HTTPException):
    code = 400 # Bad Request

    def __init__(self, description):
        super(InvalidUsage, self).__init__(description)

class FloatConverter(NumberConverter):
    regex = r"-?\d+(\.\d+)?" # match negative floats as well as integers
    num_convert = float

    def __init__(self, map, min=None, max=None):
        NumberConverter.__init__(self, map, 0, min, max)

app.url_map.converters["float"] = FloatConverter


@app.route("/")
@app.route("/<path:path>")
def index(path="index.html"):
    return send_from_directory(os.path.join(app.root_path, "static"), path)

@app.route("/api/parameters/")
def get_parameters():
    """Return the values of main parameters.

    :>json float R0: distance to the Galactic Center
    :>json string R0_unit: unit of R0
    :>json float M0: mass of the central black hole
    :>json string M0_unit: unit of M0

    .. sourcecode:: json

        {
            "R0": 8340.0,
            "R0_unit": "pc",
            "M0": 4400000.0,
            "M0_unit": "solMass"
        }
    """
    return jsonify({"R0": R0, "R0_unit": R0_unit.name, "M0": M0, "M0_unit": M0_unit.name})

@app.route("/api/stars/")
def get_stars():
    """Return a list of stars for which information is available. If a star
    has multiple names, it will appear multiple times in this list.

    :>json list names: all star names that are recognized

    .. sourcecode:: json

        {"names": ["S9", "S8", "S0-18", "S0-19", "S0-14", "S2", "..." ]}
    """
    return jsonify({"names": list(stars.keys())})

@app.route("/api/stars/<string:name>/")
def get_star(name):
    """Return the orbital elements or the proper motion parameters of a
    particular star.

    :param string name: name of the star

    **Case 1**

    :>json string id: unique id of the star
    :>json string data_type: "orbit"
    :>json float a: semi-major axis
    :>json float e: eccentricity
    :>json float inc: inclination
    :>json float Omega: longitude of ascending node
    :>json float omega: argument of pericenter
    :>json float tp: time of pericenter

    **Example**

    .. sourcecode:: html

        /api/stars/S2

    .. sourcecode:: json

        {
            "id": "S2",
            "data_type": "orbit",
            "a": 0.123,
            "e": 0.88,
            "inc": 2.36056,
            "Omega": 3.9338,
            "omega": 1.10933,
            "tp": 2002.32
        }

    **Case 2**

    :>json string id: unique id of the star
    :>json string data_type: "proper motion"
    :>json float t0: reference epoch
    :>json float x0: x-position at t0
    :>json float y0: y-position at t0
    :>json float vx0: x-velocity at t0
    :>json float vy0: y-velocity at t0
    :>json float ax0: x-acceleration at t0
    :>json float ay0: y-acceleration at t0

    **Example**

    .. sourcecode:: html

        /api/stars/S7

    .. sourcecode:: json

        {
            "id": "S7",
            "data_type": "proper_motion",
            "t0": 2000.41,
            "x0": -0.0298,
            "y0": 0.5332,
            "vx0": -0.00308,
            "vy0": -0.00403,
            "ax0": null,
            "ay0": null
        }

    :raises InvalidUsage: if the star is not found
    """
    try:
        star = stars[name]
    except:
        raise InvalidUsage("A star with this name was not found.")
    if star.orbit is not None:
        orbit = star.orbit
        return jsonify(
        {
            "id": star.id,
            "data_type": "orbit",
            "a": orbit.a,
            "e": orbit.e,
            "inc": orbit.inc,
            "Omega": orbit.Omega,
            "omega": orbit.omega,
            "tp": orbit.tp
        })
    elif star.promo is not None:
        promo = star.promo
        return jsonify(
        {
            "id": star.id,
            "data_type": "proper_motion",
            "x0": promo.x0,
            "y0": promo.y0,
            "vx0": promo.vx0,
            "vy0": promo.vy0,
            "ax0": promo.ax0,
            "ay0": promo.ay0,
            "t0": promo.t0
        })
    else:
        return jsonify({"id": star.id})

@app.route("/api/stars/<string:name>/<float:t>")
def get_location(name, t):
    """Return the location of a particular star at a certain time.

    :param string name: name of the star
    :param float t: time of observation

    :>json string id: unique id of the star
    :>json float t: time of observation
    :>json float x: x-position at t
    :>json float y: y-position at t
    :>json float z: z-position at t
    :>json float vx0: x-velocity at t
    :>json float vy0: y-velocity at t
    :>json float vz0: z-velocity at t ("radial velocity")

    **Example**

    .. sourcecode:: html

        /api/stars/S2/2016.0

    .. sourcecode:: json

        {
            "id": "S2",
            "t": 2016.0,
            "x": 0.0722,
            "y": -0.0659,
            "z": -0.0969,
            "vx": -0.0343,
            "vy": 0.0032,
            "vz": 0.0265
        }

    :raises InvalidUsage: if the star is not found
    :raises InvalidUsage: if the star can not be located
    """
    try:
        star = stars[name]
    except:
        raise InvalidUsage("A star with this name was not found.")
    try:
        particle = star.locate(t)
    except Exception as exc:
        raise InvalidUsage(str(exc))
    else:
        return jsonify(
        {
            "id": name,
            "t": t,
            "x": particle.x,
            "y": particle.y,
            "z": particle.z,
            "vx": particle.vx,
            "vy": particle.vy,
            "vz": particle.vz
        })

@app.route("/api/names/<string:name>")
def get_name(name):
    """Return alternative names for a particular star.

    :param string name: name of the star

    :>json string id: unique id of the star
    :>json list names: alternative names for the star

    **Example**

    .. sourcecode:: html

        /api/names/S0-2

    .. sourcecode:: json

        {
            "id": "S2",
            "names": ["S0-2", "S2"]
        }

    :raises InvalidUsage: if the star is not found
    """
    try:
        star = stars[name]
    except:
        raise InvalidUsage("A star with this name was not found.")
    else:
        return jsonify({"id": star.id, "names": star.names})

@app.route("/api/period/<string:name>")
def get_period(name):
    """Return the orbital period of a particular star.

    :param string name: name of the star

    :>json string id: unique id of the star
    :>json float period: orbital period of the star

    **Example**

    .. sourcecode:: html

        /api/period/S2

    .. sourcecode:: json

        {
            "id": "S2",
            "period": 15.78
        }

    :raises InvalidUsage: if the star is not found
    :raises InvalidUsage: if the star has no measured orbit
    """
    try:
        star = stars[name]
    except:
        raise InvalidUsage("A star with this name was not found.")
    try:
        period = star.period()
    except Exception as exc:
        raise InvalidUsage(str(exc))
    else:
        return jsonify({"id": star.id, "period": period})

@app.route("/api/distance/<string:name1>/<string:name2>/<float:t>")
def get_distance(name1, name2, t):
    """Return the on-sky separation of a pair of stars at a certain time.

    :param string name1: name of the first star
    :param string name2: name of the second star
    :param float t: time of observation

    :>json float distance_onsky: on-sky separation of the two stars at t
    :>json list vector_onsky: on-sky difference vector of the stellar
        positions at t

    **Example**

    .. sourcecode:: html

        /api/distance/IRS16C/S2/2016.0

    .. sourcecode:: json

        {
            "distance_onsky": 1.1967,
            "vector_onsky": [0.5005, 1.0870]
        }

    :raises InvalidUsage: if either star is not found
    :raises InvalidUsage: if either star can not be located
    """
    try:
        star1 = stars[name1]
        star2 = stars[name2]
    except:
        raise InvalidUsage("A star with this name was not found.")
    try:
        particle1 = star1.locate(t)
        particle2 = star2.locate(t)
    except Exception as exc:
        raise InvalidUsage(str(exc))
    else:
        dx = particle1.x-particle2.x
        dy = particle1.y-particle2.y
        return jsonify({"vector_onsky": [dx, dy], "distance_onsky": np.sqrt(dx**2 + dy**2)})

@app.route("/api/orbit_convert/", methods=["POST"])
def orbit_convert():
    """Convert between orbital elements and positions/velocities.

    **Case 1**

    :<json float t: time of observation
    :<json float a: semi-major axis
    :<json float e: eccentricity
    :<json float inc: inclination
    :<json float Omega: longitude of ascending node
    :<json float omega: argument of pericenter
    :<json float tp: time of pericenter

    :>json float t: time of observation
    :>json float x: x-position at t
    :>json float y: y-position at t
    :>json float z: z-position at t
    :>json float vx0: x-velocity at t
    :>json float vy0: y-velocity at t
    :>json float vz0: z-velocity at t ("radial velocity")

    **Case 2**

    :<json float t: time of observation
    :<json float x: x-position at t
    :<json float y: y-position at t
    :<json float z: z-position at t
    :<json float vx0: x-velocity at t
    :<json float vy0: y-velocity at t
    :<json float vz0: z-velocity at t ("radial velocity")

    :>json float a: semi-major axis
    :>json float e: eccentricity
    :>json float inc: inclination
    :>json float Omega: longitude of ascending node
    :>json float omega: argument of pericenter
    :>json float tp: time of pericenter
    """
    req = request.get_json()
    if all([key in req for key in
            ("t", "a", "e", "inc", "Omega", "omega", "tp")]):
        particle = orbit_to_particle(Orbit(**req), black_hole, req["t"])
        return jsonify(
        {
            "t": req["t"],
            "x": particle.x,
            "y": particle.y,
            "z": particle.z,
            "vx": particle.vx,
            "vy": particle.vy,
            "vz": particle.vz
        })
    elif all([key in req for key in
            ("t", "x", "y", "z", "vx", "vy", "vz")]):
        orbit = particle_to_orbit(Particle(**req), black_hole, req["t"])
        return jsonify(
        {
            "a": orbit.a,
            "e": orbit.e,
            "inc": orbit.inc,
            "Omega": orbit.Omega,
            "omega": orbit.omega,
            "tp": orbit.tp
        })
    else:
        raise InvalidUsage("Unable to perform the conversion.")

@app.route("/api/unit_convert/<float:value>")
def unit_convert(value):
    """Convert between different quantities.

    :param float value: value to convert

    :query string from: unit of value
    :query string to: target unit

    :>json float result: value in target unit
    :>json string result_unit: unit of result

    *Supported source units:*

        * arcsec - arc second
        * arcsec/yr - arc second per year

    *Supported target units:*

        * rad - radian
        * deg - degree
        * au - astronomical unit
        * lyr - light year
        * rg - Schwarzschild radius of the central black hole
        * km/s - kilometer per second

    **Example**

    .. sourcecode:: html

        /api/unit_convert/1.0?from=arcsec&to=pc

    .. sourcecode:: json

        {
          "result": 0.04028801690020244,
          "result_unit": "pc"
        }

    :raises InvalidUsage: if the conversion fails due to unsupported or
        incompatible units
    """
    arg_src = request.args.get("from")
    arg_dst = request.args.get("to")

    # Only these units are currently supported.
    src_units = ("arcsec", "arcsec/yr")
    dst_units = ("rad", "deg", "au", "pc", "lyr", "rg", "km/s")

    if arg_src not in src_units:
        raise InvalidUsage(f"Unit {arg_src} is unspported.")
    if arg_dst not in dst_units:
        raise InvalidUsage(f"Unit {arg_dst} is unspported.")

    # If necessary, convert angular to physical units using R0.
    if arg_src == "arcsec" and arg_dst in ("au", "pc", "lyr", "rg"):
        src_unit = Unit(arg_src).to(units.rad) * (R0 * R0_unit)
    elif arg_src == "arcsec/yr" and arg_dst in ("km/s"):
        src_unit = Unit(arg_src).to(units.rad / units.yr) * (R0 * R0_unit / units.yr)
    else:
        src_unit = Unit(arg_src)

    if arg_dst == "rg":
        rg = 2 * (M0 * M0_unit) * constants.G / constants.c**2
        try:
            result = float((value * src_unit) / rg)
            return jsonify({"result": result, "result_unit": "rg"})
        except Exception as exc:
            raise InvalidUsage(str(exc))

    try:
        dst_unit = Unit(arg_dst)
        result = (value*src_unit).to(dst_unit)
    except Exception as exc:
        raise InvalidUsage(str(exc))
    else:
        return jsonify({"result": result.value, "result_unit": f"{result.unit}"})
