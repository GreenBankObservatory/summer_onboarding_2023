import os
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack
import astropy.coordinates as coord
import astropy.units as u
import matplotlib.pyplot as plt
from tqdm import tqdm

# Table.read, use table object to select columns by names
# vstack method
# use Figure object instead of plt.scatter (module level) for better parallel stuff
# type hinting


def get_raj_decj_data(filepath: str):
    """Returns the RA and DEC data from a given FITS file"""
    # hdu = fits.open(filepath)
    # data = hdu[2].data
    # RAJ2000_data = data.field("RAJ2000")
    # DECJ2000_data = data.field("DECJ2000")
    table = Table.read(filepath, hdu=2)
    RAJ2000_data = table["RAJ2000"]
    DECJ2000_data = table["DECJ2000"]
    # hdu.close()
    return RAJ2000_data, DECJ2000_data


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
    for filename in Path(directory).glob('*'):
        # filepath = os.path.join(directory, filename)
        raj, decj = get_raj_decj_data(filename)
        # stacked_raj.extend(data[0].tolist())
        # stacked_decj.extend(data[1].tolist())
        stacked_raj.append(raj)
        stacked_decj.append(decj)
    print(len(stacked_decj))
    return vstack(stacked_raj), vstack(stacked_decj)


def expanded_stack_tables(rootdir: str):
    stacked_raj = []
    stacked_decj = []
    paths = Path(rootdir).glob("**/*.fits")
    for path in tqdm(paths):
        print(path)
        if os.path.isfile(path) and 'RAJ2000' in Table.read(path, hdu=2).columns:
            raj, decj = get_raj_decj_data(str(path))
            stacked_raj.append(raj)
            stacked_decj.append(decj)
    return vstack(stacked_raj), vstack(stacked_decj)


def separated_stack_tables(rootdir: str):
    tables = {}
    for subdir in tqdm(os.listdir(rootdir)):
        project = str(subdir)[:11]
        # project_match = re.match("AGBT\\w{3}_\\w{3}", str(subdir))
        paths = Path(subdir).glob("**/*.fits")
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
    stacked_raj, stacked_decj = expanded_stack_tables(directory)
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
    color_number = 0
    for project in tables:
        ra = coord.Angle(project["RAJ"] * u.degree)
        ra = ra.wrap_at(180 * u.degree)
        dec = coord.Angle(project["DECJ"] * u.degree)
        ax.scatter(ra.radian, dec.radian, s=0.5, c=f"C{color_number}")
        color_number += 1
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
    # projection("/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/")
    projection('/home/scratch/tchamber/antenna_fits/')


if __name__ == "__main__":
    main()
