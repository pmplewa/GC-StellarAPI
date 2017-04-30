# Galactic Center Stellar Motion API

This web service can predict the exact position and velocities of many stars in the Galactic Center
at any particular time, based on the measurements by [Gillessen et al. (2009, ApJ 692/1075)](http://dx.doi.org/10.1088/0004-637X/692/2/1075).

It is currently hosted at the following URL, including some documentation and examples:
```
https://s-stars.herokuapp.com
```

### Local Install

Please make sure you are running Python version 3:
```
python --version
```
First, download the code:
```
git clone https://github.com/pmplewa/GC-StellarAPI
cd GC-StellarAPI
```
Second, install any missing required packages:
```
pip install -r requirements.txt
```
Third, start the server:
```
python server.py
```
Finally, point your web browser to this URL to access the documentation:
```
http://localhost:5000/
```
