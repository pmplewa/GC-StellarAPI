import json
from pkg_resources import resource_filename

from .orbit_tools import Particle, Body


black_hole = Particle(x=0, y=0, z=0, vx=0, vy=0, vz=0, mass=1)

stars = {}
with open(resource_filename(__package__, "data.json"), "r") as data_file:
    data = json.load(data_file)
    for star in data:
        for name in star["names"]:
            stars[name] = Body(star["id"], primary=black_hole, **star)
