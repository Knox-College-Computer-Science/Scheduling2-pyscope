from mercuryschedule import create_blocks, main
from astropy.time import Time
from astroplan import FixedTarget
import tqdm
from astropy import units as u
from astroplan.scheduling import (Transitioner)
from astroplan import  AirmassConstraint
from astropy.table import Table
import astropy.units as u
from astroplan import TimeConstraint, Transitioner
from astropy.time import Time
import astropy.units as u

catalog=None,
queue=None,
add_only=False,
ignore_order=False,
date=None,
length=1,
observatory=None,
max_altitude=-12,
elevation=30,
airmass=3,
moon_separation=30,
scheduler=("", ""),
gap_time=60,
resolution=5,
name_format="{code}_{target}_{filter}_{exposure}s_{start_time}",
filename=None,
telrun=False,
plot=None,
yes=False,
quiet=False,
verbose=0,
reconfig_file=None,



def basic_schedule_table_with_transitions(blocks, start_time, end_time, transitioner):
    """
    Naïve, sequential scheduler for astroplan.ObservingBlock that
    explicitly inserts transition "blocks" (slew + filter change) between observations,
    and returns a Table with columns:
      target, start time (UTC), end time (UTC),
      duration (minutes), ra, dec, configuration
    """
    entries = []
    current = start_time.copy()
    prev = None

    for blk in blocks:
        # --- transition block ---
        if prev is not None:
            try:
                dt = transitioner.estimate_transition_time(prev, blk)
            except Exception:
                dt = 5 * u.second

            t0 = current
            t1 = t0 + dt
            if t1 > end_time:
                break

            # configuration for transition
            prev_f = prev.configuration.get('filter')
            next_f = blk.configuration.get('filter')
            cfg = [f"filter:{prev_f} to {next_f}"] if prev_f and next_f else []

            entries.append({
                'target':      'TransitionBlock',
                'start time (UTC)': t0.iso,
                'end time (UTC)':   t1.iso,
                'duration (minutes)': dt.to(u.minute).value,
                'ra':           '',
                'dec':          '',
                'configuration': cfg
            })
            current = t1

        # --- science block ---
        t0 = current
        t1 = t0 + blk.duration
        if t1 > end_time:
            break

        entries.append({
            'target':      blk.target.name,
            'start time (UTC)': t0.iso,
            'end time (UTC)':   t1.iso,
            'duration (minutes)': blk.duration.to(u.minute).value,
            'ra':           blk.target.coord.ra.deg,
            'dec':          blk.target.coord.dec.deg,
            'configuration': blk.configuration
        })

        current = t1
        prev = blk

    # build table with desired column order
    return Table(rows=entries, names=[
        'target',
        'start time (UTC)',
        'end time (UTC)',
        'duration (minutes)',
        'ra',
        'dec',
        'configuration'
    ])


target = [
    FixedTarget.from_name("Altair"),  # Will be high in the sky
    FixedTarget.from_name("Deneb"),   # Will be high in the sky
    FixedTarget.from_name("Polaris"), # Will be high in the sky
    FixedTarget.from_name("Acrux"),   # Will have high airmass
    FixedTarget.from_name("Canopus"), # Will have high airmass
    FixedTarget.from_name("Sirius"),  # Will have high airmass
]

# 1) Define your observing window:
SCHEDULE_START = Time('2025-05-09 01:30:00', scale='utc')
SCHEDULE_END   = Time('2025-05-09 10:30:00', scale='utc')

# 2) Specify filters & exposures:
filters = ['B', 'G', 'R']
exposures = {
    'Altair': (60 * u.second, 16),
    'Deneb': (60 * u.second, 16),
    'Polaris': (60 * u.second, 16),
    'Acrux': (100 * u.second, 16),
    'Canopus': (100 * u.second, 16),
    'Sirius': (100 * u.second, 16)
}
read_out = 20*u.second

# 3) Build constraints for block creation:
time_constraint = TimeConstraint(SCHEDULE_START, SCHEDULE_END)
airmass_constraint = AirmassConstraint(max=3.0)  # Only allow airmass <= 3.0

# 4) Generate your ObservingBlock list with constraints:
blocks = create_blocks(
    targets=target,
    constraint=[time_constraint, airmass_constraint],  # Add airmass constraint
    filters=filters,
    read_out=read_out,
    exposures=exposures
)

# 5) Build your Transitioner:
slew_rate = 0.8 * u.deg / u.second
transition_cfg = {
    'filter': {
        ('B','G'): 10*u.second,
        ('G','R'): 10*u.second,
    },
    'default': 30*u.second
}
transitioner = Transitioner(slew_rate, transition_cfg)

# 6) Now call the basic scheduler with the blocks:
print("\nBasic Scheduler Results (no airmass constraint):")
basic_table = basic_schedule_table_with_transitions(
    blocks,
    SCHEDULE_START,
    SCHEDULE_END,
    transitioner
)
print(basic_table)

# 7) Call mercury scheduler (which respects airmass constraint):
print("\nMercury Scheduler Results (with airmass constraint):")
output = main(target,
             Time('2025-05-09 01:30:00', scale='utc'), 
             Time('2025-05-09 10:30:00', scale='utc'))
print(output)

# Print comparison summary
print("\nComparison Summary:")
print("1. Basic Scheduler:")
print("   - Schedules all observations regardless of airmass")
print("   - May include observations with airmass > 3.0")
print("   - Total observations:", len(basic_table))

print("\n2. Mercury Scheduler:")
print("   - Respects airmass constraint (max 3.0)")
print("   - Skips observations with high airmass")
print("   - Total observations:", len(output) if output is not None else 0)

