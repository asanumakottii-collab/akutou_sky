"""Command-line adapters for the library; computation itself remains in memory."""
import argparse
from dataclasses import fields
import json
from pathlib import Path
import re

from ._version import __version__
from .config import ModelConfig, PRESETS

def simulation_parser(*, output_root=None, tag='sky', cache_dir=None):
    parser = argparse.ArgumentParser(description='球面大気による朝焼け・夕焼け・薄明')
    _simulation_arguments(parser, output_root, tag, cache_dir)
    return parser

def _simulation_arguments(parser, output_root=None, tag='sky', cache_dir=None):
    parser.add_argument('--tag', default=tag, help='result directory name under --output-root')
    parser.add_argument('--output', type=Path, help='explicit result directory, overriding --tag')
    parser.add_argument('--output-root', type=Path, default=output_root or Path('results'))
    parser.add_argument('--preset', choices=PRESETS, default='standard')
    parser.add_argument('--config', type=Path, help='settings JSON; overrides preset')
    parser.add_argument('--d', type=float, nargs='+', default=[-5,-2,0,2,4,6,9,12,15,18])
    parser.add_argument('--elevations', type=float, nargs='+')
    parser.add_argument('--azimuths', type=float, nargs='+', default=[0])
    parser.add_argument('--background', type=Path)
    parser.add_argument('--cache-dir', type=Path, default=cache_dir)
    parser.add_argument('--offline', action='store_true', help='refuse atmospheric data downloads')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--no-export', action='store_true', help='save NetCDF and provenance only')
    parser.add_argument('--quiet', action='store_true')
    for field in fields(ModelConfig):
        option = '--'+field.name.replace('_','-')
        if field.name == 'solar_refraction':
            parser.add_argument(option, action=argparse.BooleanOptionalAction, default=None)
        else:
            help_text = 'HG only' if field.name in ('angstrom_exponent','aerosol_ssa','aerosol_g') else None
            parser.add_argument(option, type=type(field.default), default=None, help=help_text)
    parser.add_argument('--aod', dest='aod550', type=float, help=argparse.SUPPRESS)
    parser.add_argument('--step', dest='wavelength_step_nm', type=float, help=argparse.SUPPRESS)

def parser(*, output_root=None, tag='sky', cache_dir=None):
    top = argparse.ArgumentParser(description='akutou_sky · 朝焼け・夕焼け・薄明の分光モデル')
    top.add_argument('--version', action='version', version=f'akutou-sky {__version__}')
    sub = top.add_subparsers(dest='command', required=True)
    sim = sub.add_parser('simulate', help='compute spectra, colour and brightness')
    _simulation_arguments(sim, output_root, tag, cache_dir)
    exp = sub.add_parser('export', help='export a saved dataset as HTML/CSV/NetCDF')
    exp.add_argument('input', type=Path)
    exp.add_argument('--output','--out', dest='output', type=Path)
    exp.add_argument('--overwrite', action='store_true')
    val = sub.add_parser('validate', help='compare grid sensitivity; not observational accuracy')
    val.add_argument('--input', required=True, type=Path)
    val.add_argument('--output', type=Path, help='directory for validation.json and cached runs')
    val.add_argument('--checks', nargs='+', choices=['spectral','vertical','horizontal','angular','iterations','moments'],
                     default=['spectral','vertical','horizontal','angular','iterations'])
    val.add_argument('--cache-dir', type=Path, default=cache_dir)
    val.add_argument('--offline', action='store_true')
    val.add_argument('--quiet', action='store_true')
    plot = sub.add_parser('plot', help='save a colour/luminance PNG (requires plot extra)')
    plot.add_argument('input', type=Path)
    plot.add_argument('--output','--out', dest='output', required=True, type=Path)
    plot.add_argument('--azimuth', type=float, default=0.)
    plot.add_argument('--overwrite', action='store_true')
    return top

def simulate_command(args):
    from .api import simulate
    from .io import export_results, save_dataset
    from .provenance import manifest
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',args.tag):
        raise ValueError('--tag must contain only letters, digits, underscores and hyphens')
    settings = dict(PRESETS[args.preset])
    if args.config is not None:
        supplied = json.loads(args.config.read_text(encoding='utf-8'))
        if not isinstance(supplied,dict):
            raise ValueError('--config must contain a JSON object')
        settings.update(supplied)
    settings.update({field.name:getattr(args,field.name) for field in fields(ModelConfig)
                     if getattr(args,field.name) is not None})
    config = ModelConfig.from_dict(settings)
    out = args.output or args.output_root/args.tag
    if out.exists() and any(out.iterdir()) and not args.overwrite:
        raise FileExistsError(f'{out} already exists; choose a new directory or --overwrite')
    progress = None if args.quiet else lambda message: print(message,flush=True)
    dataset = simulate(args.d,config=config,elevations=args.elevations,azimuths=args.azimuths,
                       background=args.background,cache_dir=args.cache_dir,
                       allow_download=not args.offline,progress=progress)
    if args.no_export:
        raw = save_dataset(dataset,out/'sky.nc',overwrite=args.overwrite)
        (out/'manifest.json').write_text(json.dumps(manifest(dataset,raw),ensure_ascii=False,indent=2)+'\n')
        config.to_json(out/'requested_config.json')
        # An explicitly replaced run must not retain presentations of the old run.
        for name in ['colors.csv','spectra.csv','bands.csv','index.html','summary.json']:
            (out/name).unlink(missing_ok=True)
    else:
        export_results(dataset,out,overwrite=args.overwrite)
    for name in ['validation.json','preview.png']:
        (out/name).unlink(missing_ok=True)
    if not args.quiet:
        print(f'Saved: {out.resolve()}',flush=True)
    return out

def main(argv=None, *, default_output_root=None, default_tag='sky', default_cache=None):
    top = parser(output_root=default_output_root,tag=default_tag,cache_dir=default_cache)
    args = top.parse_args(argv)
    try:
        if args.command == 'simulate':
            simulate_command(args)
        elif args.command == 'export':
            from .io import load_dataset, export_results
            out = export_results(load_dataset(args.input),args.output or args.input.parent,overwrite=args.overwrite)
            print(f'Exported: {out.resolve()}')
        elif args.command == 'validate':
            from .validation import validate_sensitivity
            output = args.output or args.input.parent
            report = validate_sensitivity(args.input,checks=args.checks,
                output_dir=output,run_directory=output/'validation',cache_dir=args.cache_dir,
                allow_download=not args.offline,progress=None if args.quiet else lambda message:print(message,flush=True))
            if not args.quiet:
                print(json.dumps({'all_requested_checks_complete':report['all_requested_checks_complete'],
                                  'passes_all_sampled_checks':report['passes_all_sampled_checks']}))
            return 0 if report['passes_all_sampled_checks'] else 1
        elif args.command == 'plot':
            if args.output.exists() and not args.overwrite:
                raise FileExistsError(f'{args.output} exists; pass --overwrite')
            from . import plot_profiles
            print(plot_profiles(args.input,args.output,azimuth=args.azimuth))
    except (ValueError, OSError, ImportError) as exc:
        top.error(str(exc))
    return 0
