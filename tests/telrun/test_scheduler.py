import warnings
from astropy import coordinates as coord
from astroplan import FixedTarget
from pyscope.telrun import mercuryschedule
from tests.utils.return_target import make_target
from astropy.time import Time
import sys
import os

# Add parent directory to path so we can import from tests
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
# Add pyscope directory to the path
pyscope_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../pyscope'))
if pyscope_path not in sys.path:
    sys.path.insert(0, pyscope_path)

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
            table = mercuryschedule.main([our_test_star], start_time, end_time, ['B', 'G', 'R'])
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
        table = mercuryschedule.main([our_test_star], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
    
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
        table = mercuryschedule.main([our_test_star], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
    
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
    start_time = Time('2025-05-09 1:30:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 10:30:00', scale='utc')  # this is end time for this test
    ra_str  = "00:03:27.15"          # hours, minutes, seconds
    dec_str = "+55:33:03.23"         # degrees, arcmin, arcsec
    name = "HD 225095"
    our_test_star = make_target(name, ra_str, dec_str)
    # calling the function with warnings suppressed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([our_test_star], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
    
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
    start_time = Time('2025-05-09 1:30:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 10:30:00', scale='utc')  # this is end time for this test
    never_visible = make_target("Acrux", "12:26:35.9", "-63:05:57")
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = mercuryschedule.main([never_visible], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
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
            table = mercuryschedule.main([always_up], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    assert length==5
    # below is the another small window just after above ( 3 hours again )
    start_time = Time('2025-05-09 4:30:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 7:30:00', scale='utc')  # this is end time for this test
    always_up = make_target("Polaris", "02:31:48.7", "+89:15:51")
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = mercuryschedule.main([always_up], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    assert length==5
    # another small window of 3 hours below
    start_time = Time('2025-05-09 7:30:00', scale='utc') #this is start time for this test
    end_time = Time('2025-05-09 10:30:00', scale='utc')  # this is end time for this test
    always_up = make_target("Polaris", "02:31:48.7", "+89:15:51")
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = mercuryschedule.main([always_up], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    assert length==5

def test_star_near_moon():
    start_time = Time('2025-05-09 1:30:00', scale='utc') #this is start time for this test around 8:30 cst
    end_time = Time('2025-05-09 10:30:00', scale='utc')  # this is end time for this test around 5:30 cst
    near_moon_star = FixedTarget.from_name("Spica")
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = mercuryschedule.main([near_moon_star], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    assert length==0 #it should not get scheduled at all



def test_two_overlapping_stars():
    start_time = Time('2025-05-09 7:00:00', scale='utc') 
    end_time = Time('2025-05-09 7:40:00', scale='utc')  
    star_one = FixedTarget.from_name("Altair")
    star_two = FixedTarget.from_name("Vega")
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = mercuryschedule.main([star_one, star_two], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    # now, according to FIFO, only Altair Star should get scheduled in that time frame
    assert length==1 #only one should get scheduled
    # And, let's check the name of the target
    assert "Altair"==table[0]['target']
    # now, let's run the same code but putting star_two as first star
    with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            table = mercuryschedule.main([star_two, star_one], start_time, end_time, ['B', 'G', 'R']) # this is running while compressing warnings
    assert table is not None
    length = len(table)
    # now, according to FIFO, only Altair Star should get scheduled in that time frame
    assert length==1 #only one should get scheduled
    # And, let's check the name of the target 
    #Take exposure time into account, don't hardcode this in as it is different for each image
    assert "Vega"==table[0]['target']
    


def test_transitioner_blocks():
    start_time = Time('2025-05-09 7:00:00', scale='utc') 
    end_time = Time('2025-05-09 10:40:00', scale='utc')
    my_star = FixedTarget.from_name("Vega")
    table = mercuryschedule.main([my_star], start_time, end_time, ['B', 'B', 'B'])
    assert table is not None
    length = len(table)
    assert length == 3 # no transitioner blocks should be there
    # currently I am scheduling Blue filter 3 times, so ideally I should have no transitioner blocks



def test_star_visible_during_window():
    """
    Test a star that is clearly visible throughout the night from UIUC.
    Vega should be visible and schedulable.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')

    star = FixedTarget.from_name("Vega")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, ['B', 'G', 'R'])

    assert table is not None
    assert len(table) > 0
    assert 5 == len(table) 

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, ['B', 'G', 'R', 'B', 'G', 'R', 'B', 'G', 'R'])

    assert table is not None
    assert len(table) > 0
    assert 17 == len(table) 


def test_star_not_visible_at_all():
    """
    Test a star that never rises at UIUC (southern hemisphere star).
    Acrux is never visible from UIUC.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')

    star = make_target("Acrux", "12:26:35.9", "-63:05:57")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, ['B', 'G', 'R'])

    assert table is not None
    assert len(table) == 0



def test_transitioner_blocks():
    """
    Test scheduling with multiple filters for the same star, expecting transition blocks.
    Vega is used and scheduled with 3 filters, so:
    - 3 observing blocks
    - 2 transition blocks
    - Total: 5 blocks
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')

    star = FixedTarget.from_name("Vega")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, ['B', 'G', 'R'])

    assert table is not None
    assert len(table) == 5


def test_star_visible_only_on_other_dates():
    """
    Star visible on some dates but not this one — Canopus (visible in winter from UIUC).
    Should not be scheduled for May 29.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')
    star = make_target("Canopus", "06:23:57.1", "-52:41:44")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, ['B', 'G', 'R'])
    assert table is not None
    assert len(table) == 0


def test_star_visible_during_daytime():
    """
    Star is visible, but only during the day (simulate this by giving a daytime UTC window).
    Vega should be invisible from 14:00–16:00 UTC (9 AM–11 AM CST).
    """
    start_time = Time('2025-05-29 14:00:00', scale='utc')
    end_time = Time('2025-05-29 16:00:00', scale='utc')
    star = FixedTarget.from_name("Vega")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, ['B'])
    assert table is not None
    assert len(table) == 0



def test_two_star_transition_blocks():
    """
    Two stars with different filters should create transition block.
    Expect 3 blocks: star1, transition, star2
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')
    s1 = FixedTarget.from_name("Vega")
    s2 = FixedTarget.from_name("Altair")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([s1, s2], start_time, end_time, ['B', 'R'])
    assert table is not None
    assert len(table) >= 3


def test_quality_score_high_pass():
    """
    Use a star expected to have high quality score (like Vega), should be scheduled.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')
    star = FixedTarget.from_name("Vega")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, ['R'])
    assert table is not None
    assert len(table) > 0


def test_quality_score_low_fail():
    """
    Target near horizon at UIUC — Sirius in late May is setting.
    Likely to have QualityScore < 3.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')
    star = make_target("Sirius", "06:45:08.9", "-16:42:58")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, ['G'])
    assert table is not None
    assert len(table) == 0


def test_conflict_same_priority_fifo():
    """
    Conflict test with equal priority stars — FIFO order matters.
    """
    start_time = Time('2025-05-29 07:00:00', scale='utc')
    end_time = Time('2025-05-29 07:30:00', scale='utc')
    s1 = FixedTarget.from_name("Altair")
    s2 = FixedTarget.from_name("Vega")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table1 = mercuryschedule.main([s1, s2], start_time, end_time, ['B'])
        table2 = mercuryschedule.main([s2, s1], start_time, end_time, ['B'])
    assert len(table1) == 1
    assert table1[0]['target'] == "Altair"
    assert len(table2) == 1
    assert table2[0]['target'] == "Vega"


def test_star_multiple_observations_fit():
    """
    Test a single star requested multiple times.
    Should fit all 6 filters.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')
    star = FixedTarget.from_name("Vega")
    filters = ['B', 'G', 'R', 'B', 'G', 'R']
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, filters)
    assert table is not None
    assert len(table) == 11  # 6 obs + 5 transitions


def test_star_mixed_with_others():
    """
    Request multiple observations of one star and singular observations of others.
    Check packing.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')
    vega = FixedTarget.from_name("Vega")
    altair = FixedTarget.from_name("Altair")
    deneb = FixedTarget.from_name("Deneb")
    filters = ['B', 'G', 'R', 'B', 'G', 'R', 'B']
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([vega, altair, deneb], start_time, end_time, filters)
    assert table is not None
    assert len(table) >= 9



START_TIME = Time('2025-05-29 01:30:00', scale='utc')
END_TIME = Time('2025-05-29 10:30:00', scale='utc')


def test_multiple_stars_different_priority():
    """Test more stars than time allows. Confirm some get skipped."""
    vega = FixedTarget.from_name("Vega")
    altair = FixedTarget.from_name("Altair")
    deneb = FixedTarget.from_name("Deneb")
    m13 = make_target("M13", "16:41:41.24", "+36:27:35.5")
    targets = [vega, altair, deneb, m13]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main(targets, START_TIME, END_TIME, ['B'] * 12)

    assert table is not None
    assert any(row['target'] == 'Vega' for row in table)


def test_star_near_setting():
    """Sirius is setting in May. Should be unschedulable."""
    sirius = make_target("Sirius", "06:45:08.9", "-16:42:58")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([sirius], START_TIME, END_TIME, ['G'])

    assert table is not None
    assert len(table) == 0


def test_multiple_observations_of_single_star():
    """Request same star with different filters = multiple observations."""
    star = FixedTarget.from_name("Altair")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], START_TIME, END_TIME, ['B', 'G', 'R', 'B'])

    assert table is not None
    assert len(table) >= 4
    assert len(table) == 7


def test_observation_within_short_window():
    """Deneb visible briefly. Should be scheduled within short night window."""
    start = Time('2025-05-29 04:00:00', scale='utc')
    end = Time('2025-05-29 05:00:00', scale='utc')
    deneb = FixedTarget.from_name("Deneb")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([deneb], start, end, ['B'])

    assert table is not None
    assert len(table) > 0


def test_repeated_same_filter():
    """No transitions should occur for repeated use of one filter."""
    polaris = make_target("Polaris", "02:31:48.7", "+89:15:51")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([polaris], START_TIME, END_TIME, ['B', 'B', 'B'])

    assert table is not None
    assert len(table) == 3


def test_multiple_targets_overlap_time():
    """Two targets at same time: only one should be scheduled."""
    start = Time('2025-05-29 07:00:00', scale='utc')
    end = Time('2025-05-29 07:40:00', scale='utc')
    vega = FixedTarget.from_name("Vega")
    altair = FixedTarget.from_name("Altair")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([vega, altair], start, end, ['B'])

    assert table is not None
    assert len(table) == 1


def test_quality_score_high_success():
    """Deneb is typically high above the horizon — should be scheduled."""
    deneb = FixedTarget.from_name("Deneb")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([deneb], START_TIME, END_TIME, ['G'])

    assert table is not None
    assert len(table) > 0


def test_observation_overlap_and_transition():
    """Two stars with different filters = expect transition blocks."""
    vega = FixedTarget.from_name("Vega")
    altair = FixedTarget.from_name("Altair")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([vega, altair], START_TIME, END_TIME, ['B', 'G'])

    assert table is not None
    assert len(table) >= 3


def test_star_set_before_observation_window():
    """Canopus sets long before observing — should not be scheduled."""
    canopus = make_target("Canopus", "06:23:57.1", "-52:41:44")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([canopus], START_TIME, END_TIME, ['R'])

    assert table is not None
    assert len(table) == 0


def test_target_always_up_multiple_filters():
    """Polaris is circumpolar — should schedule cleanly with multiple filters."""
    polaris = make_target("Polaris", "02:31:48.7", "+89:15:51")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([polaris], START_TIME, END_TIME, ['B', 'G', 'R'])

    assert table is not None
    assert len(table) >= 3

def test_default_filters_used_when_none():
    """
    If filters=None is passed, the scheduler should use the default filter list.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')
    star = FixedTarget.from_name("Vega")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, None)  # filters = None

    assert table is not None
    assert len(table) > 0

def test_invalid_star_name():
    """
    If an invalid star is passed (e.g., not resolvable), scheduler should handle gracefully.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 10:30:00', scale='utc')
    invalid_star = make_target("FakeStarXYZ", "00:00:00", "+00:00:00")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([invalid_star], start_time, end_time, ['B'])

    assert table is None or len(table) == 0


def test_star_with_skipped_filters_due_to_time():
    """
    If time allows only one or two filters out of many, schedule only those.
    """
    start_time = Time('2025-05-29 01:30:00', scale='utc')
    end_time = Time('2025-05-29 02:10:00', scale='utc')  # 30 minutes
    star = FixedTarget.from_name("Vega")
    filters = ['B', 'G', 'R', 'B', 'G']
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = mercuryschedule.main([star], start_time, end_time, filters)

    assert table is not None
    assert 1 <= len(table) < len(filters)