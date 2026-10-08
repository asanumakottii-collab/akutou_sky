"""One-variable-at-a-time numerical sensitivity checks; not observational validation."""
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
import xarray as xr
from .model import compute_model
from .provenance import sha256
from .config import ModelConfig, source_grid
from .colorimetry import xyz_from_radiance

NUMERICAL_FIELDS = {'wavelength_step_nm','nodes','nsza','ground_step_m','middle_step_m','upper_step_m',
                    'horizontal_sampling',
                    'source_ground_step_m','source_middle_step_m','source_upper_step_m','max_source_samples','max_view_source_product','max_atmosphere_source_product',
                    'iterations','relative_tolerance','absolute_tolerance','legendre_moments','wavelength_chunk','num_threads'}

def comparison(reference, candidate):
    for key in ['depression','azimuth','elevation']:
        if not np.array_equal(reference[key], candidate[key]):
            raise ValueError(f'Mismatched {key} coordinates; comparison refused')
    cr = json.loads(reference.attrs['config_json'])
    cc = json.loads(candidate.attrs['config_json'])
    physical_changes = [key for key in cr if key not in NUMERICAL_FIELDS and cr[key] != cc[key]]
    if physical_changes:
        raise ValueError(f'Convergence comparison changes physics: {physical_changes}')
    def values(ds):
        xyz = xyz_from_radiance(ds.solar_stokes_radiance.sel(stokes='I').values,
                               ds.wavelength.values, ds.wavelength_bounds.values)
        den = xyz.sum(axis=-1, keepdims=True)
        xy = np.divide(xyz[...,:2], den, out=np.full_like(xyz[...,:2],np.nan), where=den>0)
        return xyz[...,1], xy
    y1, xy1 = values(reference)
    y2, xy2 = values(candidate)
    relative = np.abs(y2-y1)/np.maximum(y1,1e-12)
    color_distance = np.linalg.norm(xy2-xy1,axis=-1)
    result = {'changed_settings':{k:[cr[k],cc[k]] for k in cr if cr[k]!=cc[k]},
              'reference_floor_cd_m2':1e-12,'per_depression':[]}
    for i,d in enumerate(reference.depression.values):
        row={'depression_deg':float(d)}
        for name,mask in [('all',np.ones(reference.sizes['elevation'],bool)),
                          ('above_5deg',reference.elevation.values>=5),
                          ('below_5deg',reference.elevation.values<5)]:
            if mask.any():
                row[name]={'max_relative_Y':float(relative[i,:,mask].max()),
                           'median_relative_Y':float(np.median(relative[i,:,mask])),
                           'max_delta_xy':float(np.nanmax(color_distance[i,:,mask]))}
        result['per_depression'].append(row)
    result['max_relative_Y'] = float(relative.max())
    result['max_delta_xy'] = float(np.nanmax(color_distance))
    result['passes_5percent_Y_and_0_005_xy'] = bool(result['max_relative_Y']<=.05 and result['max_delta_xy']<=.005)
    return result

DEFAULT_CHECKS = ('spectral','vertical','horizontal','angular','iterations')
SUPPORTED_CHECKS = (*DEFAULT_CHECKS, 'moments')

def validate_sensitivity(input, *, checks=DEFAULT_CHECKS, output_dir=None,
                         run_directory=None, cache_dir=None, allow_download=True,
                         progress=None):
    """Compare refined grids on sampled rays; return a numerical-sensitivity report.

    Input may be a dataset or a NetCDF path. No result files are written when
    output_dir is None. If supplied, output_dir gets validation.json and a runs/
    cache (or the explicit run_directory). This does not validate observations.
    """
    from .io import load_dataset
    checks = tuple(checks)
    if not checks or len(set(checks)) != len(checks) or any(name not in SUPPORTED_CHECKS for name in checks):
        raise ValueError(f'checks must be unique names from {SUPPORTED_CHECKS}')
    input_path = None if isinstance(input, xr.Dataset) else Path(input)
    full = input if input_path is None else load_dataset(input_path)
    # Subset existing reference coordinates; no interpolation or fabricated samples.
    def nearest(axis, targets):
        values=full[axis].values
        return np.unique([values[np.argmin(abs(values-v))] for v in targets])
    reference=full.sel(depression=nearest('depression',[0,6,12,18]),
                       elevation=nearest('elevation',[0,1,5,15,45,90]))
    config=ModelConfig.from_dict(json.loads(reference.attrs['config_json']))
    next_nodes=next((v for v in [14,26,38,50,74,86,110,146,170,194,230,266,302] if v>config.nodes),None)
    variants={
        'spectral': {'wavelength_step_nm':config.wavelength_step_nm/2},
        'vertical': {'ground_step_m':config.ground_step_m/2,'middle_step_m':config.middle_step_m/2,
                     'upper_step_m':config.upper_step_m/2,
                     'source_ground_step_m':config.source_ground_step_m/2,
                     'source_middle_step_m':config.source_middle_step_m/2,
                     'source_upper_step_m':config.source_upper_step_m/2},
        'horizontal': {'nsza':2*config.nsza-1},
        'angular': {'nodes':next_nodes},
        'iterations': {'iterations':2*config.iterations,'relative_tolerance':config.relative_tolerance/10,
                       'absolute_tolerance':config.absolute_tolerance/10},
        'moments': {'legendre_moments':2*config.legendre_moments},
    }
    # Refine as far as the configured memory budget allows; report exact changes.
    affordable_columns = config.max_source_samples//(len(source_grid(config))*config.nodes*config.num_stokes)
    variants['horizontal']['nsza'] = min(2*config.nsza-1,affordable_columns)
    for factor in [2.,1.5,1.25,1.1]:
        vertical = dict(variants['vertical'])
        for key in ['source_ground_step_m','source_middle_step_m','source_upper_step_m',
                    'ground_step_m','middle_step_m','upper_step_m']:
            vertical[key] = getattr(config,key)/factor
        try:
            replace(config,**vertical).validate()
            variants['vertical']=vertical
            break
        except ValueError:
            continue
    output_dir = None if output_dir is None else Path(output_dir)
    if run_directory is not None and output_dir is None:
        raise ValueError('run_directory requires output_dir')
    out = None if output_dir is None else Path(run_directory) if run_directory is not None else output_dir/'runs'
    if out is not None:
        out.mkdir(parents=True,exist_ok=True)
        output_dir.mkdir(parents=True,exist_ok=True)
    report={'scope':'Sampled angles only, numerical sensitivity, not an error bound or observational validation',
            'input': '<in-memory Dataset>' if input_path is None else str(input_path.resolve()),
            'input_sha256':None if input_path is None else sha256(input_path),'reference_config':config.to_dict(),
            'sample_coordinates':{k:reference[k].values.tolist() for k in ['depression','azimuth','elevation']},'checks':{}}
    def save_report():
        report['all_requested_checks_complete'] = (len(report['checks'])==len(checks)
            and all(v.get('status')=='completed' for v in report['checks'].values()))
        report['passes_all_sampled_checks'] = (report['all_requested_checks_complete']
            and all(v['passes_5percent_Y_and_0_005_xy'] for v in report['checks'].values()))
        if output_dir is not None:
            (output_dir/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    for name in checks:
        try:
            if name=='angular' and next_nodes is None:
                raise ValueError('No supported higher angular rule available')
            if name=='horizontal' and variants[name]['nsza']<=config.nsza:
                raise ValueError('No room for a finer horizontal grid within the resource budget')
            variant=replace(config,**variants[name]).validate()
        except ValueError as exc:
            report['checks'][name]={'status':'skipped','reason':str(exc),
                                    'passes_5percent_Y_and_0_005_xy':False}
            save_report()
            continue
        path=None if out is None else out/f'{name}.nc'
        candidate=None
        if path is not None and path.exists():
            with xr.open_dataset(path) as cached:
                if (cached.attrs.get('transport_source_sha256') == reference.attrs.get('transport_source_sha256')
                    and cached.attrs['physical_inputs_json'] == reference.attrs['physical_inputs_json']
                    and json.loads(cached.attrs['config_json'])==variant.to_dict() and all(
                    np.array_equal(cached[k],reference[k]) for k in ['depression','azimuth','elevation'])):
                    candidate=cached.load()
        if candidate is None:
            candidate=compute_model(reference.depression.values,variant,reference.elevation.values,reference.azimuth.values,
                                    progress,cache_dir=cache_dir,allow_download=allow_download)
            if path is not None:
                from .io import save_dataset
                save_dataset(candidate,path,overwrite=True)
        report['checks'][name]=comparison(reference,candidate)
        report['checks'][name]['status']='completed'
        save_report()
        if progress:
            progress(f"{name}: max relative Y={report['checks'][name]['max_relative_Y']:.6g}, delta xy={report['checks'][name]['max_delta_xy']:.6g}")

    return report
