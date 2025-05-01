import configparser
import itertools
import logging
import os

import astroplan
import numpy as np
from astropy import coordinates as coord
from astropy.table import Table
from astropy import time as astrotime
from astropy import units as u

logger = logging.getLogger(__name__)


def blocks_to_table(blocks):
    """Convert a list of observing blocks to an astropy table."""
    # Initialize lists for each column
    ids = []
    names = []
    start_times = []
    end_times = []
    targets = []
    priorities = []
    observers = []
    codes = []
    titles = []
    filenames = []
    types = []
    backends = []
    filters = []
    exposures = []
    nexp = []
    repositionings = []
    shutter_states = []
    readouts = []
    binnings = []
    frame_positions = []
    frame_sizes = []
    pm_ra_cosdecs = []
    pm_decs = []
    comments = []
    schs = []
    statuses = []
    messages = []
    sched_times = []

    # Process each block
    for block in blocks:
        # Handle None values for target coordinates
        if block.get('target') is None:
            target_coords = "0h0m0.0s -90d0m0.0s"
        else:
            target_coords = f"{block['target_ra']} {block['target_dec']}"

        # Handle proper motion values
        pm_ra = block.get('pm_ra_cosdec', 0)
        pm_dec = block.get('pm_dec', 0)
        
        # Convert proper motion values to arcsec/hour if they're not already in the right format
        if isinstance(pm_ra, (int, float)):
            pm_ra = pm_ra * u.arcsec / u.hour
        if isinstance(pm_dec, (int, float)):
            pm_dec = pm_dec * u.arcsec / u.hour

        # Append values to lists
        ids.append(block.get('ID', 0))
        names.append(block.get('name', ''))
        start_times.append(block.get('start_time', None))
        end_times.append(block.get('end_time', None))
        targets.append(target_coords)
        priorities.append(block.get('priority', 0))
        observers.append(block.get('observer', []))
        codes.append(block.get('code', ''))
        titles.append(block.get('title', ''))
        filenames.append(block.get('filename', ''))
        types.append(block.get('type', ''))
        backends.append(block.get('backend', 0))
        filters.append(block.get('filter', ''))
        exposures.append(block.get('exposure', 0))
        nexp.append(block.get('nexp', 0))
        repositionings.append(block.get('repositioning', [0, 0]))
        shutter_states.append(block.get('shutter_state', False))
        readouts.append(block.get('readout', 0))
        binnings.append(block.get('binning', [1, 1]))
        frame_positions.append(block.get('frame_position', [0, 0]))
        frame_sizes.append(block.get('frame_size', [0, 0]))
        pm_ra_cosdecs.append(pm_ra)
        pm_decs.append(pm_dec)
        comments.append(block.get('comment', ''))
        schs.append(block.get('sch', ''))
        statuses.append(block.get('status', ''))
        messages.append(block.get('message', ''))
        sched_times.append(block.get('sched_time', None))

    # Create the table
    result_table = Table()
    result_table['ID'] = ids
    result_table['name'] = names
    result_table['start_time'] = start_times
    result_table['end_time'] = end_times
    result_table['target'] = targets
    result_table['priority'] = priorities
    result_table['observer'] = observers
    result_table['code'] = codes
    result_table['title'] = titles
    result_table['filename'] = filenames
    result_table['type'] = types
    result_table['backend'] = backends
    result_table['filter'] = filters
    result_table['exposure'] = exposures
    result_table['nexp'] = nexp
    result_table['repositioning'] = repositionings
    result_table['shutter_state'] = shutter_states
    result_table['readout'] = readouts
    result_table['binning'] = binnings
    result_table['frame_position'] = frame_positions
    result_table['frame_size'] = frame_sizes
    result_table['pm_ra_cosdec'] = pm_ra_cosdecs
    result_table['pm_dec'] = pm_decs
    result_table['comment'] = comments
    result_table['sch'] = schs
    result_table['status'] = statuses
    result_table['message'] = messages
    result_table['sched_time'] = sched_times

    return result_table


def table_to_blocks(table):
    blocks = []
    for row in table:
        # parse the constraints
        constraints = []
        for constraint in row["constraints"]:
            try:
                if constraint["type"] == "TimeConstraint":
                    constraints.append(
                        astroplan.TimeConstraint(
                            min=astrotime.Time(constraint["min"]),
                            max=astrotime.Time(constraint["max"]),
                        )
                    )
                elif constraint["type"] == "AtNightConstraint":
                    constraints.append(
                        astroplan.AtNightConstraint(
                            max_solar_altitude=constraint["max_solar_altitude"] * u.deg
                        )
                    )
                elif constraint["type"] == "AltitudeConstraint":
                    constraints.append(
                        astroplan.AltitudeConstraint(
                            min=constraint["min"] * u.deg,
                            max=constraint["max"] * u.deg,
                            boolean_constraint=constraint["boolean_constraint"],
                        )
                    )
                elif constraint["type"] == "AirmassConstraint":
                    constraints.append(
                        astroplan.AirmassConstraint(
                            min=constraint["min"],
                            max=constraint["max"],
                            boolean_constraint=constraint["boolean_constraint"],
                        )
                    )
                elif constraint["type"] == "MoonSeparationConstraint":
                    constraints.append(
                        astroplan.MoonSeparationConstraint(
                            min=constraint["min"] * u.deg,
                            max=constraint["max"] * u.deg,
                        )
                    )
                else:
                    logger.warning("Only time constraints are currently supported")
                    continue
            except:
                constraints.append(None)

        if row["ID"] is None:
            row["ID"] = astrotime.Time.now().mjd

        blocks.append(
            astroplan.ObservingBlock(
                target=astroplan.FixedTarget(row["target"]),
                duration=row["exposure"] * row["nexp"] * u.second,
                priority=row["priority"],
                name=row["name"],
                configuration={
                    "observer": row["observer"],
                    "code": row["code"],
                    "title": row["title"],
                    "filename": row["filename"],
                    "type": row["type"],
                    "backend": row["backend"],
                    "filter": row["filter"],
                    "exposure": row["exposure"],
                    "nexp": row["nexp"],
                    "repositioning": row["repositioning"],
                    "shutter_state": row["shutter_state"],
                    "readout": row["readout"],
                    "binning": row["binning"],
                    "frame_position": row["frame_position"],
                    "frame_size": row["frame_size"],
                    "pm_ra_cosdec": row["pm_ra_cosdec"],
                    "pm_dec": row["pm_dec"],
                    "comment": row["comment"],
                    "sch": row["sch"],
                    "ID": row["ID"],
                    "status": row["status"],
                    "message": row["message"],
                    "sched_time": row["sched_time"],
                },
                constraints=constraints,
            )
        )

    return blocks


def validate(schedule_table, observatory=None):
    logger.info("Validating the schedule table")

    if observatory is None:
        logger.info(
            "No observatory was specified, so validation will only check for basic formatting errors."
        )

    convert_to_blocks = False
    if type(schedule_table) is list:
        logger.info("Converting list of blocks to astropy table")
        schedule_table = blocks_to_table(schedule_table)
        convert_to_blocks = True

    assert (
        type(schedule_table) is Table or type(schedule_table) is Table.Row
    ), "schedule_table must be an astropy table or row"

    if type(schedule_table) is Table.Row:
        schedule_table = Table(schedule_table)

    # Check for required columns
    required_columns = [
        "name",
        "start_time",
        "end_time",
        "target",
        "priority",
        "observer",
        "code",
        "title",
        "filename",
        "type",
        "backend",
        "filter",
        "exposure",
        "nexp",
        "repositioning",
        "shutter_state",
        "readout",
        "binning",
        "frame_position",
        "frame_size",
        "pm_ra_cosdec",
        "pm_dec",
        "comment",
        "sch",
        "ID",
        "status",
        "message",
        "sched_time",
        "constraints",
    ]
    for column in required_columns:
        if column not in schedule_table.columns:
            logger.error(f"Column {column} is missing")
            raise ValueError(f"Column {column} is missing")

    # Check dtypes
    for colname in schedule_table.colnames:
        column = schedule_table[colname]
        match colname:
            case (
                "name"
                | "observer"
                | "code"
                | "title"
                | "filename"
                | "type"
                | "filter"
                | "comment"
                | "sch"
                | "status"
                | "message"
            ):
                if not np.issubdtype(column.dtype, np.dtype("U")):
                    logger.error(
                        f"Column '{column.name}' must be of type str, not {column.dtype}"
                    )
                    raise ValueError(
                        f"Column '{column.name}' must be of type str, not {column.dtype}"
                    )
            case "start_time" | "end_time":
                if type(column) is not astrotime.Time:
                    logger.error(
                        f"Column '{column.name}' must be of type astropy.time.Time, not {type(column)}"
                    )
                    raise ValueError(
                        f"Column '{column.name}' must be of type astropy.time.Time, not {type(column)}"
                    )
            case "target":
                if type(column) is not coord.SkyCoord:
                    logger.error(
                        f"Column '{column.name}' must be of type astropy.coordinates.SkyCoord, not {type(column)}"
                    )
                    raise ValueError(
                        f"Column '{column.name}' must be of type astropy.coordinates.SkyCoord, not {type(column)}"
                    )
            # case (
            #     "priority"
            #     | "nexp"
            #     | "readout"
            #     | "frame_position"
            #     | "frame_size"
            #     | "binning"
            #     | "repositioning"
            # ):
            #     if not np.issubdtype(column.dtype, np.dtype("int64")):
            #         logger.error(
            #             f"Column '{column.name}' must be of type int64, not {column.dtype}"
            #         )
            #         raise ValueError(
            #             f"Column '{column.name}' must be of type int64, not {column.dtype}"
            #         )
            case "exposure" | "pm_ra_cosdec" | "pm_dec":
                if not np.issubdtype(column.dtype, np.floating):
                    logger.error(
                        f"Column '{column.name}' must be of a float type, not {column.dtype}"
                    )
                    raise ValueError(
                        f"Column '{column.name}' must be of a float type, not {column.dtype}"
                    )
            case "shutter_state":
                if column.dtype != bool:
                    logger.error(
                        f"Column '{column.name}' must be of type bool, not {column.dtype}"
                    )
                    raise ValueError(
                        f"Column '{column.name}' must be of type bool, not {column.dtype}"
                    )

    # Obs-specific validation
    if observatory is not None:
        logger.info("Performing observatory-specific validation")
        for row in schedule_table:
            logger.info(f"Validating row {row.index}")

            if row["name"] == "TransitionBlock" or row["name"] == "EmptyBlock":
                logger.info(f"Skipping validation of {row['name']}")
                continue

            # Logging to debug and verify the input values
            logger.info(f"Target object: {row['target']}, Type: {type(row['target'])}")
            logger.info(
                f"Start time: {row['start_time']}, Type: {type(row['start_time'])}"
            )

            # Check if target is observable at start time
            altaz_obj = observatory.get_object_altaz(
                obj=row["target"],
                t=row["start_time"],
            )
            logger.info(f"AltAz Object: {altaz_obj}")
            if altaz_obj.alt < observatory.min_altitude:
                logger.error("Target is not observable at start time")
                row["status"] = "I"  # Invalid
                row["message"] = "Target is not observable at start time"
                continue

            # Check if source is observable at end time
            altaz_obj = observatory.get_object_altaz(
                obj=row["target"],
                t=row["end_time"],
            )
            if altaz_obj.alt < observatory.min_altitude:
                logger.error("Target is not observable at end time")
                row["status"] = "I"  # Invalid
                row["message"] = "Target is not observable at end time"
                continue

            # Check filter
            if len(observatory.filters) == 0:
                logger.info("No filters available, no filter check performed")
            elif row["filter"] not in observatory.filters:
                logger.error("Requested filter is not available")
                row["status"] = "I"  # Invalid
                row["message"] = "Requested filter is not available"
                continue

            # Check exposure time
            try:
                current_cam_state = observatory.camera.Connected
                observatory.camera.Connected = True
                if row["exposure"].to(u.second).value > observatory.camera.ExposureMax:
                    logger.error("Exposure time exceeds maximum")
                    row["status"] = "I"  # Invalid
                    row["message"] = "Exposure time exceeds maximum"
                    continue
                elif (
                    row["exposure"].to(u.second).value < observatory.camera.ExposureMin
                ):
                    logger.error("Exposure time is below minimum")
                    row["status"] = "I"  # Invalid
                    row["message"] = "Exposure time is below minimum"
                    continue
                observatory.camera.Connected = current_cam_state
            except:
                logger.warning(
                    "Exposure time range check failed because the driver is not available"
                )

            # Check repositioning, frame position, and frame size
            try:
                current_cam_state = observatory.camera.Connected
                observatory.camera.Connected = True
                if (
                    (
                        row["repositioning"][0] > observatory.camera.CameraXSize
                        or row["repositioning"][1] > observatory.camera.CameraYSize
                    )
                    and row["repositioning"] != [0, 0]
                    and row["repositioning"] != [None, None]
                ):
                    logger.error("Repositioning coordinates exceed camera size")
                    row["status"] = "I"  # Invalid
                    row["message"] = "Repositioning coordinates exceed camera size"
                if (
                    (
                        row["frame_position"][0] + row["frame_size"][0]
                        > observatory.camera.CameraXSize
                        or row["frame_position"][1] + row["frame_size"][1]
                        > observatory.camera.CameraYSize
                    )
                    and row["frame_position"] != [0, 0]
                    and row["frame_size"] != [0, 0]
                ):
                    logger.error("Frame position and size exceed camera size")
                    row["status"] = "I"  # Invalid
                    row["message"] = "Frame position and size exceed camera size"
                observatory.camera.Connected = current_cam_state
            except:
                logger.warning(
                    "Repositioning check failed because the driver is not available"
                )

            # Check readout
            try:
                current_cam_state = observatory.camera.Connected
                observatory.camera.Connected = True
                if row["readout"] >= len(observatory.camera.ReadoutModes):
                    logger.error("Readout mode not available")
                    row["status"] = "I"  # Invalid
                    row["message"] = "Readout mode not available"
                observatory.camera.Connected = current_cam_state
            except:
                logger.warning(
                    "Readout mode check failed because the driver is not available"
                )

            # Check binning
            try:
                current_cam_state = observatory.camera.Connected
                observatory.camera.Connected = True
                if row["binning"][0] > observatory.camera.MaxBinX:
                    logger.error("Binning exceeds maximum in X")
                    row["status"] = "I"  # Invalid
                    row["message"] = "Binning exceeds maximum in X"
                if row["binning"][1] > observatory.camera.MaxBinY:
                    logger.error("Binning exceeds maximum in Y")
                    row["status"] = "I"
                    row["message"] = "Binning exceeds maximum in Y"
                if (
                    row["binning"][0] != row["binning"][1]
                    and not observatory.camera.CanAsymmetricBin
                ):
                    logger.error("Binning must be square")
                    row["status"] = "I"
                    row["message"] = "Binning must be square"
                observatory.camera.Connected = current_cam_state
            except:
                logger.warning(
                    "Binning check failed because the driver is not available"
                )

    if convert_to_blocks:
        logger.info("Converting astropy table back to list of blocks")
        schedule_table = table_to_blocks(schedule_table)
    return schedule_table


def _mask_expander(arr, mask):
    return np.array([[mask[i]] * len(arr[i]) for i in range(len(arr))]).ravel()
