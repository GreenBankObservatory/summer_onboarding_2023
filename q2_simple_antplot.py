import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

def get_raj_decj_data(filepath):
    hdu = fits.open('/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/2014_11_20_01:05:05.fits')
    data = hdu[2].data
    RAJ2000_data = data.field('RAJ2000')
    DECJ2000_data = data.field('DECJ2000')
    hdu.close()
    return RAJ2000_data, DECJ2000_data

def print_values():
    RAJ2000_data, DECJ2000_data = get_raj_decj_data(filepath)
    print(f"The lowest and highest values of the RAJ2000 data are {RAJ2000_data.min()} and {RAJ2000_data.max()}. The lowest and highest values of the DECJ2000 data are {DECJ2000_data.min()} and {DECJ2000_data.max()}.")

def plot_one_file(filepath):
    data = get_raj_decj_data(filepath)
    plt.scatter(data[1], data[0])
    plt.title("RAJ2000 vs. DECJ2000")
    plt.xlabel("deg")
    plt.ylabel("deg")
    plt.show()

def stack_tables(directory):
    stacked_raj = []
    stacked_decj = []
    for filename in os.listdir(directory):
        filename = os.path.join(directory, filename)
        data = get_raj_decj_data(filepath)
        stacked_raj.append(data[0])
        stacked_decj.append(data[1])
    return stacked_raj, stacked_decj
  

def plot_all_files(directory):
    stacked_raj, stacked_decj = stack_tables(directory)
    plt.scatter(stacked_decj, stacked_raj, s=0.5)
    plt.xlabel("deg")
    plt.ylabel("deg")
    plt.show()

def projection(directory):
    stacked_raj, stacked_decj = stack_tables(directory)
    ra = coord.Angle(stacked_raj.filled(np.nan)*u.degree)
    ra = ra.wrap_at(180*u.degree)
    dec = coord.Angle(stacked_decj.filled(np.nan)*u.degree)
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.scatter(ra.radian, dec.radian)
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    ax.grid(True)

def main():
    # plot_one_file('/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/2014_11_20_01:05:05.fits')
    plot_all_files('/home/scratch/kwei/raw_data/AGBT13B_312_34/Antenna/')

if __name__ == "__main__":
    main()
