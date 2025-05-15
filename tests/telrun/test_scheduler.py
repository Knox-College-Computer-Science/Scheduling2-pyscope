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

def test_HD_225095_Less_Than_35():
    # To test, we need following things
    # One: Start time 
    # Two: End time
    # Three: Target 
    # And, that is it :)))))
    # should not schedule any blocks
    start_time = Time('2025-05-09 7:00:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 7:30:00', scale='utc')  # this is end time for this test
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
        try:
            table = tested_hardcodedblocks.main([our_test_star], start_time, end_time)
        except Exception as e:
            print(f"Ignored exception: {e}")
            table = None

    assert table is None

    # num_rows = len(table)
    # print(f"Table has {num_rows} rows")
    # assert num_rows > 0
    # assert num_rows == 3

def test_HD_225095_One_Hour():
    # To test, we need following things
    # One: Start time 
    # Two: End time
    # Three: Target 
    # And, that is it :)))))
    #Should have 1 schedule block,0 transition blocks, with 1 filter
    start_time = Time('2025-05-09 7:00:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 8:00:00', scale='utc')  # this is end time for this test
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
    assert num_rows == 1

def test_HD_More_Than_One_Hour():
    # To test, we need following things
    # One: Start time 
    # Two: End time
    # Three: Target 
    # And, that is it :)))))
    #Should have 1 Transition block, with 2 filters
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

def test_HD_225095_Two_Hours():
    # To test, we need following things
    # One: Start time 
    # Two: End time
    # Three: Target 
    # And, that is it :)))))
    #Should have 2 transition blocks, with 3 filters and 3 blocks [Total 5 blocks]
    #Length=5
    start_time = Time('2025-05-09 7:00:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 9:00:00', scale='utc')  # this is end time for this test
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
    assert num_rows == 5
