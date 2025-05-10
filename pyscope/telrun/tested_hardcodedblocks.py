from astroplan import Observer, FixedTarget, ObservingBlock
from astroplan.constraints import (AirmassConstraint, AtNightConstraint,
                                   TimeConstraint, is_observable, is_always_observable)
from astroplan.scheduling import (Transitioner, PriorityScheduler, SequentialScheduler, Schedule)
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from zoneinfo import ZoneInfo
import numpy as np

# Use Central Time (Chicago timezone) with automatic DST handling
central = ZoneInfo("America/Chicago")

# ===================== CONFIGURATION =====================

SCHEDULE_START = Time('2025-05-09 7:00:00', scale='utc')
SCHEDULE_END = Time('2025-05-09 8:15:00', scale='utc')

print(SCHEDULE_START.to_datetime(timezone=central))
print(SCHEDULE_END.to_datetime(timezone=central))
READ_OUT = 20 * u.second
DENEB_EXP = 60 * u.second
M13_EXP = 100 * u.second
EXPOSURE_COUNT = 16
# FILTERS = ['B', 'G', 'R']
FILTERS = ['B','G','R']

# ===================== SETUP FUNCTIONS =====================

def setup_observer():
    # Coordinates for University of Illinois Astronomical Observatory
    illinois_location = EarthLocation(lat=40.1106*u.deg, lon=-88.2260*u.deg, height=222*u.m)
    return Observer(location=illinois_location, name="UIUC Observatory", timezone="America/Chicago")

def setup_targets():
    # our_star = SkyCoord(ra=00:03:27.15*u.deg, dec=55:33:03.23*u.deg)
    ra_str  = "00:03:27.15"          # hours, minutes, seconds
    dec_str = "+55:33:03.23"         # degrees, arcmin, arcsec

    our_coord = SkyCoord(ra=ra_str,
                        dec=dec_str,
                        unit=(u.hourangle, u.deg),   # RA in hours, Dec in degrees
                        frame="icrs")
    our_new_star = FixedTarget(coord=our_coord, name="HD 225095")
    # return [
    #     FixedTarget.from_name('Altair'),
    #     FixedTarget.from_name('Vega'),
    #     FixedTarget.from_name('Deneb'),
    #     FixedTarget.from_name('M13')
    # ]
    return [our_new_star]

def setup_constraints():
    # return [
        # AirmassConstraint(max=3, boolean_constraint=False),
    #     AtNightConstraint.twilight_civil()
    # ]
    return [
        AirmassConstraint(max=3, boolean_constraint=False),
        AtNightConstraint.twilight_civil()
    ]

def setup_transitioner():
    slew_rate = 0.8 * u.deg / u.second
    return Transitioner(slew_rate, {
        'filter': {
            ('B', 'G'): 10 * u.second,
            ('G', 'R'): 10 * u.second,
            'default': 30 * u.second
        }
    })

# ===================== HELPER FUNCTIONS =====================

def print_visibility_info(observer, targets, start_time, end_time):
    print("=== Sun & Moon Info ===")
    anchor_time = Time('2025-05-09 00:00:00', scale='utc')  # Anchor to local solar noon
    sunrise = observer.sun_rise_time(anchor_time, which='next')
    sunset = observer.sun_set_time(anchor_time, which='next')
    moonrise = observer.moon_rise_time(anchor_time, which='next')
    moonset = observer.moon_set_time(anchor_time, which='next')

    for label, time_obj in zip(
        ["Sunrise", "Sunset", "Moonrise", "Moonset"],
        [sunrise, sunset, moonrise, moonset]
    ):
        print(f"{label}: {time_obj.to_datetime(timezone=central)}")

    print("=== Target Visibility Summary ===")
    time_grid = start_time + np.arange(0, (end_time - start_time).to(u.minute).value, 10) * u.minute
    constraints = setup_constraints()

    for target in targets:
        print(f"\nTarget: {target.name}")
        up_mask = observer.target_is_up(time_grid, target)
        obs_mask = is_observable(constraints, observer, [target], times=time_grid)
        always_obs = is_always_observable(constraints, observer, [target], times=time_grid)

        print(f"  - Is up at start? {observer.target_is_up(start_time, target)}")
        print(f"  - Observable at some point? {obs_mask[0]}")
        print(f"  - Always observable? {always_obs[0]}")
        print(f"  - Rise: {observer.target_rise_time(start_time, target).to_datetime(timezone=central)}")
        print(f"  - Set: {observer.target_set_time(start_time, target).to_datetime(timezone=central)}")
        print(f"  - Transit: {observer.target_meridian_transit_time(start_time, target).to_datetime(timezone=central)}")

# ===================== SCHEDULING =====================

def create_blocks(targets, constraint):
    blocks = []
    exposures = {
        'Altair': (60 * u.second, 16),
        'Vega': (60 * u.second, 16),
        'Deneb': (60 * u.second, 16),
        'M13': (100 * u.second, 16),
        'HD 225095': (100 * u.second, 16),

    }

    for priority, band in enumerate(FILTERS):
        for target in targets:
            exp_time, n = exposures[target.name]
            blocks.append(ObservingBlock.from_exposures(
                target, priority, exp_time, n, READ_OUT,
                configuration={'filter': band}, constraints=[constraint]))
    return blocks

def run_scheduler(scheduler_class, observer, blocks, start_time, end_time, constraints, transitioner):
    schedule = Schedule(start_time, end_time)
    scheduler = scheduler_class(observer=observer, constraints=constraints, transitioner=transitioner)
    scheduler(blocks, schedule)
    return schedule

# ===================== MAIN =====================

def main():
    observer = setup_observer()
    targets = setup_targets()
    constraints = setup_constraints()
    transitioner = setup_transitioner()

    print_visibility_info(observer, targets, SCHEDULE_START, SCHEDULE_END)

    time_constraint = TimeConstraint(SCHEDULE_START, SCHEDULE_END)
    blocks = create_blocks(targets, time_constraint)

    print("\n=== Running Priority Scheduler ===")
    priority_schedule = run_scheduler(PriorityScheduler, observer, blocks, SCHEDULE_START, SCHEDULE_END, constraints, transitioner)
    print(priority_schedule.to_table())

if __name__ == '__main__':
    main()
