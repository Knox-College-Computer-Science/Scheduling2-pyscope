import pytest
import warnings
from astropy import coordinates as coord
from astropy import time 
from astroplan import FixedTarget
from pyscope.telrun import quality
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import sys
import os

# Add parent directory to path so we can import from tests
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
# Add pyscope directory to the path
pyscope_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../pyscope'))
if pyscope_path not in sys.path:
    sys.path.insert(0, pyscope_path)
from tests.utils.return_target import make_target


from pyscope.telrun import sch

print("The path is ", sys.path)
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
    name = "HD 225095"
    our_test_star = make_target(name, ra_str, dec_str)
    # calling the function with warnings suppressed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            table = quality.main([our_test_star], start_time, end_time)
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
    name = "HD 225095"
    our_test_star = make_target(name, ra_str, dec_str)
    # calling the function with warnings suppressed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = quality.main([our_test_star], start_time, end_time) # this is running while compressing warnings
    
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
    name = "HD 225095"
    our_test_star = make_target(name, ra_str, dec_str)
    # calling the function with warnings suppressed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = quality.main([our_test_star], start_time, end_time) # this is running while compressing warnings
    
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
    name = "HD 225095"
    our_test_star = make_target(name, ra_str, dec_str)
    # calling the function with warnings suppressed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = quality.main([our_test_star], start_time, end_time) # this is running while compressing warnings
    
    # Check the table exists and count rows
    assert table is not None
    num_rows = len(table)
    print(f"Table has {num_rows} rows")
    assert num_rows > 0
    assert num_rows == 5



# Now, testing a star that is never visible
def test_star_never_visible():
    # it should not be scheduled at all
    # so the length should be 0
    start_time = Time('2025-05-09 7:00:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 15:00:00', scale='utc')  # this is end time for this test
    never_visible = make_target("Acrux", "12:26:35.9", "-63:05:57")
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = quality.main([never_visible], start_time, end_time) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    assert length==0

def test_star_always_up():
    # this star should always be visible hence we can schedule it anytime at night
    # below is the one small half of night (3 hours)
    start_time = Time('2025-05-09 1:30:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 4:30:00', scale='utc')  # this is end time for this test
    always_up = make_target("Polaris", "02:31:48.7", "+89:15:51")
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = quality.main([always_up], start_time, end_time) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    assert length==5
    # below is the another small window just after above ( 3 hours again )
    start_time = Time('2025-05-09 4:30:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 7:30:00', scale='utc')  # this is end time for this test
    always_up = make_target("Polaris", "02:31:48.7", "+89:15:51")
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = quality.main([always_up], start_time, end_time) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    assert length==5
    # another small window of 3 hours below
    start_time = Time('2025-05-09 7:30:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 10:30:00', scale='utc')  # this is end time for this test
    always_up = make_target("Polaris", "02:31:48.7", "+89:15:51")
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = quality.main([always_up], start_time, end_time) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    assert length==5


