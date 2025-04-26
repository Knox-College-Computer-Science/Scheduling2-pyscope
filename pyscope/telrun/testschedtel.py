import sch
import schedtab
import astroplan
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import time as astrotime

# Set the schedule time
sched_time = astrotime.Time.now()

# Read the .sch file
block_groups = []
block_groups.append(
    sch.read(r"C:\Users\Santosh\Desktop\0-100\cs322\Sprint 2\v2\Scheduling2-pyscope\pyscope\telrun\schedules\1070 Project Images SP25.sch")
)

# Add IDs and sched_time to all blocks
for group in block_groups:
    for block in group:
        block.setdefault("ID", astrotime.Time.now().mjd)
        block["sched_time"] = sched_time

# Flatten the block_groups
all_blocks = [block for group in block_groups for block in group]

# Generate table
table = schedtab.blocks_to_table(all_blocks)

# 📅 Auto-name file based on current time
now = astrotime.Time.now()
formatted_time = now.strftime("%Y-%m-%dT%H-%M-%S")
output_path = f"telrun_{formatted_time}.ecsv"

# Save the table
table.write(output_path, format="ascii.ecsv", overwrite=True)

print(f"Saved table to {output_path} ✅")


