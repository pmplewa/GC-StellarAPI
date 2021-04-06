import numpy as np
import pandas as pd
import pytest
import requests


url = "http://localhost:5000/api"

def test_parameters():
    r = requests.get(f"{url}/parameters")
    assert r.status_code == 200
    assert r.headers["content-type"] == "application/json"
    json = r.json()
    assert "M0" in json and "M0_unit" in json
    assert "R0" in json and "R0_unit" in json

def test_names():
    r = requests.get(f"{url}/names/S2")
    assert r.status_code == 200
    json = r.json()
    assert "names" in json and "S0-2" in json["names"]

@pytest.mark.parametrize("star", ["S2", "S4"])
def test_orbit(star):
    r = requests.get(f"{url}/stars/{star}")
    assert r.status_code == 200
    orbit = r.json()
    assert "data_type" in orbit and orbit["data_type"] == "orbit"

    r = requests.get(f"{url}/period/{star}")
    assert r.status_code == 200
    period = r.json()["period"]
    assert period > 0

    t_min = orbit["tp"] - period / 2
    t_max = orbit["tp"] + period / 2
    t_val = np.linspace(t_min, t_max, 42)
    data = pd.DataFrame([requests.get(f"{url}/stars/{star}/{t}").json() for t in t_val])
    assert len(data) == len(t_val)
