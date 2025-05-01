import sch
import schedtab
import astroplan
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import time as astrotime
from astropy.time import Time
from astroplan import FixedTarget, Observer, Transitioner, ObservingBlock

noon_before = Time('2025-04-30 19:00')
noon_after = Time('2025-05-01 19:00')

from astroplan.constraints import AtNightConstraint, AirmassConstraint

# create the list of constraints that all targets must satisfy
global_constraints = [AirmassConstraint(max = 3, boolean_constraint = False),
                      AtNightConstraint.twilight_civil()]

# Set the schedule time
sched_time = astrotime.Time.now()

# Read the .sch file
block_groups = []
# Define priorities for each university
university_priorities = {
    "uiowa.edu": 1,  # University of Iowa - Highest priority
    "augustana.edu": 2,  # Augustana College
    "macalester.edu": 3,  # Macalester College
    "knox.edu": 4,  # Knox College
    "Other": 5  # Default priority for unknown universities
}

# Read and assign priorities to each schedule
for schedule_path in [
    r"C:\Users\dibya\Documents\GitHub\Scheduling2-pyscope\pyscope\telrun\schedules\Be.sch",
    r"C:\Users\dibya\Documents\GitHub\Scheduling2-pyscope\pyscope\telrun\schedules\Sch Files (Valid)\xpgQSOTest.sch",
    r"C:\Users\dibya\Documents\GitHub\Scheduling2-pyscope\pyscope\telrun\schedules\Sch Files (Valid)\alc_m57-2025-03-12 (1).sch",
    r"C:\Users\dibya\Documents\GitHub\Scheduling2-pyscope\pyscope\telrun\schedules\Sch Files (Valid)\xpgQSOs.sch"
]:
    blocks = sch.read(schedule_path)
    # Extract university from the first block's observer email
    if blocks and 'observer' in blocks[0]:
        observer_email = blocks[0]['observer']
        if isinstance(observer_email, list):
            observer_email = observer_email[0]  # Take first email if multiple observers
        university_domain = observer_email.split('@')[-1] if '@' in observer_email else 'Other'
    else:
        university_domain = 'Other'
    
    # Assign priority based on university domain
    priority = university_priorities.get(university_domain, university_priorities['Other'])
    
    # Assign priority to all blocks in the group
    for block in blocks:
        block['priority'] = priority
        block['university'] = university_domain
    
    block_groups.append(blocks)

print("sch parser did not gave any errors")
# Add IDs and sched_time to all blocks
for group in block_groups:
    for block in group:
        block.setdefault("ID", astrotime.Time.now().mjd)
        block["sched_time"] = sched_time

# Convert dictionary blocks to ObservingBlock objects
observing_blocks = []
for group in block_groups:
    for block in group:
        # Create FixedTarget from coordinates
        target = FixedTarget(coord=SkyCoord(block['target_ra']*u.deg, block['target_dec']*u.deg), 
                           name=block['name'])
        
        # Create ObservingBlock with the assigned priority
        obs_block = ObservingBlock(
            target=target,
            duration=block['duration'],
            priority=block['priority'],  # Use the university-based priority
            name=block['name'],
            configuration=block,
            constraints=block['constraints']
        )
        observing_blocks.append(obs_block)

# Sort observing blocks by priority (lower number = higher priority)
observing_blocks.sort(key=lambda x: x.priority)

# Initialize the observer and targets, and create observing blocks
apo = Observer.at_site('apo')
transitioner = Transitioner(slew_rate=2*u.deg/u.second)

seq_scheduler = astroplan.SequentialScheduler(constraints = [],
                                    observer = apo,
                                    transitioner = transitioner)

sequential_schedule = astroplan.Schedule(noon_before, noon_after)

# Print initial state
print("\nInitial Schedule:")
print(f"Start time: {sequential_schedule.start_time.iso}")
print(f"End time: {sequential_schedule.end_time.iso}")
print(f"Number of slots: {len(sequential_schedule.slots)}")

# Schedule the blocks
seq_scheduler(observing_blocks, sequential_schedule)

# Print final state
print("\nFinal Schedule:")
print(f"Start time: {sequential_schedule.start_time.iso}")
print(f"End time: {sequential_schedule.end_time.iso}")
print(f"Number of slots: {len(sequential_schedule.slots)}")

# Count different types of blocks
total_blocks = 0
observing_blocks = 0
transition_blocks = 0
university_counts = {univ: 0 for univ in university_priorities.keys()}
priority_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}  # Track blocks by priority

print("\nScheduled Blocks:")
for i, slot in enumerate(sequential_schedule.slots):
    if slot.block is not None:
        total_blocks += 1
        print(f"\nSlot {i}:")
        print(f"  Start: {slot.start.iso}")
        print(f"  End: {slot.end.iso}")
        if isinstance(slot.block, astroplan.ObservingBlock):
            observing_blocks += 1
            priority = slot.block.priority
            priority_counts[priority] += 1
            university = slot.block.configuration.get('university', 'Other')
            university_counts[university] += 1
            print(f"  Type: Observing Block")
            print(f"  Name: {slot.block.name}")
            print(f"  Target: {slot.block.target.name}")
            print(f"  Duration: {slot.block.duration.to(u.minute)}")
            print(f"  Priority: {priority}")
            print(f"  University: {university}")
        else:
            transition_blocks += 1
            print(f"  Type: Transition Block")
            print(f"  Duration: {(slot.end - slot.start).to(u.minute)}")

print(f"\nSchedule Summary:")
print(f"Total blocks: {total_blocks}")
print(f"Observing blocks: {observing_blocks}")
print(f"Transition blocks: {transition_blocks}")
print("\nBlocks by University:")
for university, count in university_counts.items():
    print(f"{university}: {count} blocks")
print("\nBlocks by Priority:")
for priority, count in priority_counts.items():
    print(f"Priority {priority}: {count} blocks")
print(f"Total observing time: {sum(slot.block.duration.to(u.minute) for slot in sequential_schedule.slots if isinstance(slot.block, astroplan.ObservingBlock))} minutes")
print(f"Total transition time: {sum((slot.end - slot.start).to(u.minute) for slot in sequential_schedule.slots if not isinstance(slot.block, astroplan.ObservingBlock))} minutes")

# Get the scheduled blocks from the schedule
scheduled_blocks = []
for slot in sequential_schedule.slots:
    if slot.block is not None:
        if isinstance(slot.block, astroplan.ObservingBlock):
            # Regular observing block
            config = slot.block.configuration
            config['start_time'] = slot.start
            config['end_time'] = slot.end
            scheduled_blocks.append(config)
        else:
            # Transition block
            transition_block = {
                'ID': astrotime.Time.now().mjd,
                'name': 'TransitionBlock',
                'start_time': slot.start,
                'end_time': slot.end,
                'target': None,
                'priority': 0,
                'observer': [],
                'code': '',
                'title': 'Transition',
                'filename': '',
                'type': 'transition',
                'backend': 0,
                'filter': '',
                'exposure': 0,
                'nexp': 0,
                'repositioning': [0, 0],
                'shutter_state': False,
                'readout': 0,
                'binning': [1, 1],
                'frame_position': [0, 0],
                'frame_size': [0, 0],
                'pm_ra_cosdec': 0,
                'pm_dec': 0,
                'comment': 'Transition between targets',
                'sch': '',
                'status': 'T',
                'message': 'Transition',
                'sched_time': sched_time
            }
            scheduled_blocks.append(transition_block)

# Generate table
table = schedtab.blocks_to_table(scheduled_blocks)

# 📅 Auto-name file based on current time
now = astrotime.Time.now()
formatted_time = now.strftime("%Y-%m-%dT%H-%M-%S")
output_path = f"telrun_{formatted_time}.ecsv"

# Save the table
table.write(output_path, format="ascii.ecsv", overwrite=True)

print(f"\nSaved table to {output_path} ✅")


