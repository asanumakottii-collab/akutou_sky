"""Read and analyze a saved Dataset without a new transport calculation."""
import argparse
from pathlib import Path
from akutou_sky import load_dataset, band_mean

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path,nargs='?',default=Path('results/improved/sky.nc'))
    args=parser.parse_args()
    sky=load_dataset(args.input)
    print(sky.total_XYZ.sel(tristimulus='Y',azimuth=0,elevation=15))
    print('XYZ mean over elevations 5..15 deg:')
    print(band_mean(sky.elevation.values,sky.total_XYZ.values,5,15))

if __name__=='__main__':
    main()
