"""Run: python examples/basic_usage.py --cache-dir model_cache --offline."""
import argparse
from pathlib import Path
from akutou_sky import ModelConfig, simulate, export_results

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir',type=Path)
    parser.add_argument('--offline',action='store_true')
    parser.add_argument('--output',type=Path,default=Path('results/library_example'))
    args=parser.parse_args()
    sky=simulate([-2,0,6,12],config=ModelConfig.from_preset('preview',aod550=.05,ozone_du=330),
                 elevations=[0,1,5,15,30,60,90],azimuths=[0,90,180],cache_dir=args.cache_dir,
                 allow_download=not args.offline,progress=print)
    export_results(sky,args.output)
    print(sky.total_XYZ.sel(tristimulus='Y'))
    print(args.output.resolve())

if __name__=='__main__':
    main()
