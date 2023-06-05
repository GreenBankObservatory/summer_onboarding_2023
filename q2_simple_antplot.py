import os
from pathlib import Path
import re

import numpy as np
from astropy.table import Table, Column, vstack, hstack
import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
from tqdm import tqdm

# Table.read, use table object to select columns by names
# vstack method
# use Figure object instead of plt.scatter (module level) for better parallel stuff
# type hinting

# don't need ** if only one directory down
# return single table for data
# stack tables
# have table containing column for project name
# access tables with masks
# string split, join
# regex - if not [result], raise ValueError
# pdb, ipdb for debugging from console
# find way to open debug upon error?


def get_raj_decj_data(filepath: str):
    """Returns the RA and DEC data from a given FITS file"""
    table = Table.read(filepath, hdu=2)
    return hstack([table["RAJ2000"], table["DECJ2000"]])


def print_values(filepath: str):
    """Prints the highest and lowest RA/DEC values from a FITS file"""
    table = get_raj_decj_data(filepath)
    ra = table["RAJ2000"]
    dec = table["DECJ2000"]
    print(
        "The lowest and highest values of the RAJ2000 data are"
        f" {ra.min()} and {ra.max()}. The lowest and highest values"
        f" of the DECJ2000 data are {dec.min()} and {dec.max()}."
    )


def plot_one_file(filepath: str):
    """Generates scatterplot of RA vs. DEC data for one file"""
    data = get_raj_decj_data(filepath)
    plt.scatter(data["DECJ2000"], data["RAJ2000"])
    plt.title("RAJ2000 vs. DECJ2000")
    plt.xlabel("deg")
    plt.ylabel("deg")
    plt.show()


def stack_tables(directory: str):
    """Returns combined tables of RA and DEC data for all files in a directory"""
    stacked_tables = []
    for filename in Path(directory).glob("*"):
        table = get_raj_decj_data(filename)
        stacked_tables.append(table)
    return vstack(stacked_tables)


def expanded_stack_tables(rootdir: str):
    """Combines RA, DEC data for files in subdirectories within a directory"""
    stacked_tables = []
    paths = list(Path(rootdir).glob("*/Antenna/*.fits"))
    for path in tqdm(paths):
        print(f"Loading {path}")
        table = get_raj_decj_data(path)
        stacked_tables.append(table)
    return vstack(stacked_tables)


def separated_stack_tables_dict(rootdir: str):
    """Returns a dictionary of RA, DEC data corresponding to different projects. Uses regular expression"""
    tables = {}
    for subdir in tqdm(os.listdir(rootdir)):
        project_match = re.match(r"[a-zA-Z0-9]+[_][\d]+[_]", str(subdir))
        if not project_match:
            raise ValueError("File does not match naming format")
        else:
            project = project_match.group(0)[:-1]
        paths = list(Path(os.path.join(rootdir, subdir)).glob("Antenna/*.fits"))
        for path in paths:
            table = get_raj_decj_data(path)
            if project not in list(tables):
                tables[project] = [table]
            else:
                tables[project].append(table)
    for project in tables:
        tables[project] = vstack(tables[project])
    return tables


def separated_stack_tables_dict_unnested(rootdir: str):
    """Same as separated_stack_tables_dict but un-nests loop through subdirs. Uses split/join"""
    tables = {}
    paths = list(Path(rootdir).glob("*/Antenna/*.fits"))
    for path in tqdm(paths):
        split = str(path).split("/")[-3].split("_")
        project = "_".join(split[:-1])
        table = get_raj_decj_data(path)
        if project not in list(tables):
            tables[project] = [table]
        else:
            tables[project].append(table)
    for project in tables:
        tables[project] = vstack(tables[project])
    return tables


def separated_stack_tables_table(rootdir: str):
    """Returns a table of RA, DEC, project name data"""
    stacked_table = []
    projects = []
    for subdir in tqdm(os.listdir(rootdir)):
        split = subdir.split("_")
        project = "_".join(split[:-1])
        paths = list(Path(os.path.join(rootdir, subdir)).glob("Antenna/*.fits"))
        for path in paths:
            table = Table.read(path, hdu=2)
            formatted_table = hstack(
                [
                    table["RAJ2000"],
                    table["DECJ2000"],
                    Column(name="Project", data=np.full(len(table), project)),
                ]
            )
            stacked_table.append(formatted_table)
            if project not in projects:
                projects.append(project)
    return vstack(stacked_table), projects


def separated_stack_tables_table_unnested(rootdir: str):
    """Same as separated_stack_tables_table but un-nests loop through subdirs"""
    stacked_table = []
    projects = []
    paths = list(Path(rootdir).glob("*/Antenna/*.fits"))
    for path in tqdm(paths):
        project_match = re.search(r"[a-zA-Z0-9]+[_][\d]+[_]", str(path))
        if not project_match:
            raise ValueError("File does not match naming format")
        else:
            project = project_match.group(0)[:-1]
        table = Table.read(path, hdu=2)
        formatted_table = hstack(
            [
                table["RAJ2000"],
                table["DECJ2000"],
                Column(name="Project", data=np.full(len(table), project)),
            ]
        )
        stacked_table.append(formatted_table)
        if project not in projects:
            projects.append(project)
    return vstack(stacked_table), projects


def plot_all_files(directory: str):
    """Generates scatterplot of RA vs. DEC data for all files in a directory"""
    table = stack_tables(directory)
    plt.scatter(table["DECJ2000"], table["RAJ2000"], s=0.5)
    plt.xlabel("deg")
    plt.ylabel("deg")
    plt.show()


def projection(directory: str):
    """Plots a Mollweide projection for stacked data in a directory"""
    table = stack_tables(directory)
    ra = coord.Angle(table["RAJ2000"])
    ra = ra.wrap_at(180 * u.degree)
    dec = coord.Angle(table["DECJ2000"])
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.scatter(ra.radian, dec.radian, s=0.5)
    ax.set_xticklabels(
        ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    )
    ax.grid(True)
    fig.show()


def expanded_projection(directory: str):
    """Plots a Mollweide projection for stacked antenna data in a directory"""
    table = expanded_stack_tables(directory)
    ra = coord.Angle(table["RAJ2000"])
    ra = ra.wrap_at(180 * u.degree)
    dec = coord.Angle(table["DECJ2000"])
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.scatter(ra.radian, dec.radian, s=0.5)
    ax.set_xticklabels(
        ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    )
    ax.grid(True)
    fig.show()


def color_projection_dict(directory: str):
    """Plots a Mollweide projection for stacked antenna data, with different colors for each project"""
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    tables = separated_stack_tables_dict(directory)
    for project in tables:
        ra = coord.Angle(tables[project]["RAJ2000"])
        ra = ra.wrap_at(180 * u.degree)
        dec = coord.Angle(tables[project]["DECJ2000"])
        ax.scatter(ra.radian, dec.radian, s=0.5)
    ax.set_xticklabels(
        ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    )
    ax.grid(True)
    fig.show()


def color_projection_dict_unnested(directory: str):
    """Same as color_projection_dict but uses separated_stack_tables_dict_unnested"""
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    tables = separated_stack_tables_dict_unnested(directory)
    for project in tables:
        ra = coord.Angle(tables[project]["RAJ2000"])
        ra = ra.wrap_at(180 * u.degree)
        dec = coord.Angle(tables[project]["DECJ2000"])
        ax.scatter(ra.radian, dec.radian, s=0.5)
    ax.set_xticklabels(
        ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    )
    ax.grid(True)
    fig.show()


def color_projection_table(directory: str):
    """Same as color_projection_dict but uses table stacking"""
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    table, projects = separated_stack_tables_table(directory)
    for project in projects:
        sliced = table[table["Project"] == project]
        ra = coord.Angle(sliced["RAJ2000"])
        ra = ra.wrap_at(180 * u.degree)
        dec = coord.Angle(sliced["DECJ2000"])
        ax.scatter(ra.radian, dec.radian, s=0.5)
    ax.set_xticklabels(
        ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    )
    ax.grid(True)
    fig.show()


def color_projection_table_unnested(directory: str):
    """Same as color_projection_table but uses separated_stack_tables_table_unnested"""
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    table, projects = separated_stack_tables_table_unnested(directory)
    for project in projects:
        sliced = table[table["Project"] == project]
        ra = coord.Angle(sliced["RAJ2000"])
        ra = ra.wrap_at(180 * u.degree)
        dec = coord.Angle(sliced["DECJ2000"])
        ax.scatter(ra.radian, dec.radian, s=0.5)
    ax.set_xticklabels(
        ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    )
    ax.grid(True)
    fig.show()


def main():
    # print_values(
    #     "/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/2014_11_20_01:05:05.fits"
    # )
    # plot_one_file(
    #     "/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/2014_11_20_01:05:05.fits"
    # )
    # plot_all_files("/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/")
    # projection("/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/")
    # expanded_projection("/home/scratch/tchamber/antenna_fits/")
    color_projection_table("/home/scratch/kwei/raw_data/")


if __name__ == "__main__":
    main()
