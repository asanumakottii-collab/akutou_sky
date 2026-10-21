"""akutou-sky command line: run, led, explore, check."""
import argparse
from pathlib import Path
import sys

from ._version import __version__


def _common(p):
    p.add_argument('--cache-dir', type=Path, default=None, help='atmospheric input data (default: see README)')
    p.add_argument('--offline', action='store_true', help='refuse to download missing input data')
    p.add_argument('--numerics', default='final', help='preview, final, or a JSON file')


def main(argv=None):
    parser = argparse.ArgumentParser(prog='akutou-sky', description='朝焼け・夕焼け・薄明の空の色と輝度')
    parser.add_argument('--version', action='version', version=__version__)
    sub = parser.add_subparsers(dest='command', required=True)

    run = sub.add_parser('run', help='compute one scene -> sky.nc, colors.csv, bands.csv')
    run.add_argument('--scene', default='standard', help='clear, standard, hazy, urban, rural, humid, dry, or JSON')
    run.add_argument('--d', type=float, nargs='+', default=[x / 2 for x in range(-12, 25)],
                     help='solar depressions in degrees (default -6..12 every 0.5)')
    run.add_argument('--elevations', type=float, nargs='+')
    run.add_argument('--azimuths', type=float, nargs='+', default=[0., 90., 180.])
    run.add_argument('--output', type=Path, default=Path('results/standard'))
    _common(run)

    led = sub.add_parser('led', help='sky.nc -> led_table.csv and akutou_sky_table.h for the ESP32')
    led.add_argument('input', type=Path, help='sky.nc from `akutou-sky run`')
    led.add_argument('--azimuth', type=float, default=0., help='viewing azimuth from the Sun (0 = sunward)')
    led.add_argument('--exponent', type=float, default=.5, help='luminance compression (1 = proportional)')
    led.add_argument('--reference-luminance', type=float, help='cd/m2 shown at full duty (default: brightest band)')
    led.add_argument('--primaries', type=Path, help='measured LED primaries JSON (default: sRGB)')
    led.add_argument('--output', type=Path, help='directory (default: next to the input)')

    explore = sub.add_parser('explore', help='four-sky dust/fog comparison page')
    explore.add_argument('--output', type=Path, default=Path('results/explorer'))
    _common(explore)

    check = sub.add_parser('check', help='numerical convergence and measured-illuminance comparison')
    check.add_argument('--scene', action='append', help='scene(s) to check (default: standard)')
    check.add_argument('--only', choices=['convergence', 'illuminance'])
    check.add_argument('--controls', nargs='+', help='subset of numerical controls to refine')
    check.add_argument('--output', type=Path, default=Path('results/checks'))
    _common(check)

    args = parser.parse_args(argv)
    say = lambda message: print(message, flush=True)
    try:
        if args.command == 'run':
            from .io import export
            from .solver import render
            ds = render(args.d, args.scene, args.numerics, args.elevations, args.azimuths,
                        cache_dir=args.cache_dir, allow_download=not args.offline, progress=say)
            say(f'saved {export(ds, args.output)}')
        elif args.command == 'led':
            from .color import add_color
            from .io import load
            from .led import led_table, load_primaries, write_csv, write_header
            ds = load(args.input)
            if 'XYZ' not in ds:
                ds = add_color(ds)
            table = led_table(ds, azimuth=args.azimuth, exponent=args.exponent,
                              reference_luminance=args.reference_luminance, primaries=load_primaries(args.primaries))
            out = args.output or args.input.parent
            out.mkdir(parents=True, exist_ok=True)
            write_csv(table, out / 'led_table.csv')
            import json
            name = json.loads(ds.attrs.get('scene_json', '{}'))
            write_header(table, out / 'akutou_sky_table.h', version=__version__,
                         scene_name=args.input.parent.name + f' (aod550={name.get("aod550")})')
            say(f'saved {out / "led_table.csv"} and {out / "akutou_sky_table.h"}')
        elif args.command == 'explore':
            from .explorer import build_explorer
            say(f'saved {build_explorer(args.output, args.numerics, cache_dir=args.cache_dir, allow_download=not args.offline, progress=say)}')
        elif args.command == 'check':
            from .check import convergence, illuminance, save_report
            ok = True
            for scene in args.scene or ['standard']:
                if args.only in (None, 'illuminance'):
                    report = illuminance(scene, args.numerics, cache_dir=args.cache_dir, allow_download=not args.offline)
                    save_report(report, args.output / f'illuminance_{Path(scene).stem}.json')
                    for row in report['rows']:
                        say(f'{scene} d={row["depression"]:+5.1f}  model {row["model_lux"]:10.4g} lx  '
                            f'USNO {row["usno_lux"]:10.4g} lx  ratio {row["ratio"]:.2f}')
                if args.only in (None, 'convergence'):
                    report = convergence(scene, args.numerics, args.controls, cache_dir=args.cache_dir,
                                         allow_download=not args.offline, progress=say)
                    save_report(report, args.output / f'convergence_{Path(scene).stem}.json')
                    ok &= report['all_pass']
            return 0 if ok else 1
    except (ValueError, FileNotFoundError, FileExistsError) as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
