import pytest
import warnings
from astropy import coordinates as coord
from astropy import time 
from astroplan import FixedTarget
from pyscope.telrun import tested_hardcodedblocks
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time


from pyscope.telrun import sch
def testHD_225095():
    # To test, we need following things
    # One: Start time 
    # Two: End time
    # Three: Target 
    # And, that is it :)))))
    start_time = Time('2025-05-09 7:00:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 8:15:00', scale='utc')  # this is end time for this test
    ra_str  = "00:03:27.15"          # hours, minutes, seconds
    dec_str = "+55:33:03.23"         # degrees, arcmin, arcsec
    our_coord = SkyCoord(ra=ra_str,
                        dec=dec_str,
                        unit=(u.hourangle, u.deg),   # RA in hours, Dec in degrees
                        frame="icrs")
    our_test_star = FixedTarget(coord=our_coord, name="HD 225095") # we made a new target here
    # calling the function with warnings suppressed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = tested_hardcodedblocks.main([our_test_star], start_time, end_time) # this is running while compressing warnings
    
    # Check the table exists and count rows
    assert table is not None
    num_rows = len(table)
    print(f"Table has {num_rows} rows")
    assert num_rows > 0
    assert num_rows == 3
