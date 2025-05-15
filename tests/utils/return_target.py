from astropy import coordinates as coord
from astropy import time 
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u

def make_target(name, ra, dec):
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame="icrs")
    return FixedTarget(coord=coord, name=name)