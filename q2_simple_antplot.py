import os
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack, hstack
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
# pdb, ipdb for debuggin from console
# find way to open debug upon error?


def get_raj_decj_data(filepath: str):
    """Returns the RA and DEC data from a given FITS file"""
    table = Table.read(filepath, hdu=2)
    # RAJ2000_data = table["RAJ2000"]
    # DECJ2000_data = table["DECJ2000"]
    # return RAJ2000_data, DECJ2000_data
    return hstack(table["RAJ2000"], table["DECJ2000"])


def print_values(filepath: str):
    """Simply prints the highest and lowest RA/DEC values from a FITS file"""
    RAJ2000_data, DECJ2000_data = get_raj_decj_data(filepath)
    print(
        "The lowest and highest values of the RAJ2000 data are"
        f" {RAJ2000_data.min()} and {RAJ2000_data.max()}. The lowest and highest values"
        f" of the DECJ2000 data are {DECJ2000_data.min()} and {DECJ2000_data.max()}."
    )


def plot_one_file(filepath: str):
    """Generates scatterplot of RA vs. DEC data for one file"""
    data = get_raj_decj_data(filepath)
    plt.scatter(data[1], data[0])
    plt.title("RAJ2000 vs. DECJ2000")
    plt.xlabel("deg")
    plt.ylabel("deg")
    plt.show()


def stack_tables(directory: str):
    """Returns combined tables of RA and DEC data for all files in a directory"""
    stacked_raj = []
    stacked_decj = []
    for filename in Path(directory).glob("*"):
        # filepath = os.path.join(directory, filename)
        raj, decj = get_raj_decj_data(filename)
        # stacked_raj.extend(data[0].tolist())
        # stacked_decj.extend(data[1].tolist())
        stacked_raj.append(raj)
        stacked_decj.append(decj)
    print(len(stacked_decj))
    return vstack(stacked_raj), vstack(stacked_decj)


def stack_tables2(directory: str):
    """Returns combined tables of RA and DEC data for all files in a directory"""
    stacked_raj = []
    stacked_decj = []
    for filename in Path(directory).glob("*"):
        # filepath = os.path.join(directory, filename)
        raj, decj = get_raj_decj_data(filename)
        stacked_raj.extend(raj)
        stacked_decj.extend(decj)
        # stacked_raj.append(raj)
        # stacked_decj.append(decj)
    return stacked_raj, stacked_decj


def expanded_stack_tables(rootdir: str):
    """Combines RA, DEC data for files in subdirectories within a directory"""
    stacked_raj = []
    stacked_decj = []
    paths = Path(rootdir).glob("*/Antenna/*.fits")
    for path in tqdm(paths):
        print(path)
        raj, decj = get_raj_decj_data(str(path))
        stacked_raj.append(raj)
        stacked_decj.append(decj)
    return vstack(stacked_raj), vstack(stacked_decj)


def expanded_stack_tables2(rootdir: str):
    """Combines RA, DEC data for files in subdirectories within a directory"""
    stacked_raj = []
    stacked_decj = []
    paths = Path(rootdir).glob("**/Antenna/*.fits")
    for path in tqdm(paths):
        raj, decj = get_raj_decj_data(path)
        stacked_raj.extend(raj)
        stacked_decj.extend(decj)
    return stacked_raj, stacked_decj


def separated_stack_tables(rootdir: str):
    """Returns a dictionary of RA, DEC data corresponding to different projects"""
    tables = {}
    for subdir in tqdm(os.listdir(rootdir)):
        project = str(subdir)[:11]
        # project_match = re.match("AGBT\\w{3}_\\w{3}", str(subdir))
        paths = Path(subdir).glob("*/Antenna/*.fits")
        for path in paths:
            raj, decj = get_raj_decj_data(str(path))
            tables[project]["RAJ"].append(raj)
            tables[project]["DECJ"].append(decj)
    return tables


def plot_all_files(directory: str):
    """Generates scatterplot of RA vs. DEC data for all files in a directory"""
    stacked_raj, stacked_decj = stack_tables(directory)
    plt.scatter(stacked_decj, stacked_raj, s=0.5)
    plt.xlabel("deg")
    plt.ylabel("deg")
    plt.show()


def projection(directory: str):
    """Plots a Mollweide projection for stacked data in a directory"""
    stacked_raj, stacked_decj = stack_tables(directory)
    ra = coord.Angle(stacked_raj * u.degree)
    ra = ra.wrap_at(180 * u.degree)
    dec = coord.Angle(stacked_decj * u.degree)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.scatter(ra.radian, dec.radian, s=0.5)
    ax.set_xticklabels(
        ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    )
    ax.grid(True)
    fig.show()


def colored_projection(directory: str):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="mollweide")
    tables = separated_stack_tables(directory)
    for project in tables:
        ra = coord.Angle(project["RAJ"] * u.degree)
        ra = ra.wrap_at(180 * u.degree)
        dec = coord.Angle(project["DECJ"] * u.degree)
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
    # plot_one_file('/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/
    # 2014_11_20_01:05:05.fits')
    # plot_all_files('/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/')
    projection("/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/")
    # projection("/home/scratch/tchamber/antenna_fits/")


if __name__ == "__main__":
    main()
