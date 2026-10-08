"""Recompute observer single scatter while preserving a verified split diffuse field.

This migration requires the archived source hash that produced the input. It
writes a new directory; the original datasets are never changed. A separate
fresh full-transport regression must be run before publishing the directory.
"""
from dataclasses import replace
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil

from akutou_sky import ModelConfig,load_dataset,save_dataset
from akutou_sky.model import replace_single_scattering
from akutou_sky.colorimetry import add_colorimetry
from akutou_sky.data import input_paths,resolve_cache_dir
from akutou_sky.explorer import _source_hashes,build_explorer
from akutou_sky._version import __version__

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--original-model-source',type=Path,required=True)
    p.add_argument('--cache-dir',type=Path,default=Path('model_cache'))
    p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    expected=hashlib.sha256(args.original_model_source.read_bytes()).hexdigest()
    manifest=json.loads((args.input/'manifest.json').read_text())
    assert expected==manifest['source_hashes']['model.py']
    assert args.input.resolve()!=args.output.resolve()
    if args.output.exists() and any(args.output.iterdir()) and not args.resume:
        raise FileExistsError('Output exists; use --resume to verify or replace migration results')
    cache=resolve_cache_dir(args.cache_dir)
    inputs={str(f.relative_to(cache)):hashlib.sha256(f.read_bytes()).hexdigest()
            for f in input_paths(cache,allow_download=False)}
    signature=json.dumps(dict(sources=_source_hashes(),inputs=inputs,
        backend=importlib.metadata.version('sasktran2'),library=__version__),sort_keys=True)
    old_cache={};new_cache={}
    def refine(path,destination,*,case=False):
        original=load_dataset(path)
        assert original.attrs['transport_source_sha256']==expected,path
        cfg=ModelConfig.from_dict(json.loads(original.attrs['config_json']))
        assert cfg.single_scatter_step_m==0 and cfg.delta_m_method=='split'
        if destination.exists() and args.resume:
            existing=load_dataset(destination)
            if (existing.attrs.get('transport_source_sha256')==_source_hashes()['model.py'] and
                existing.attrs.get('diffuse_component_reused_from_file_sha256')==hashlib.sha256(path.read_bytes()).hexdigest() and
                json.loads(existing.attrs['config_json'])==replace(cfg,single_scatter_step_m=500.).to_dict()):
                return
        updated=replace_single_scattering(original,replace(cfg,single_scatter_step_m=500.),
            cache_dir=cache,allow_download=False,coarse_cache=old_cache,fine_cache=new_cache)
        updated.attrs['diffuse_component_reused_from_file_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
        updated.attrs['diffuse_component_reused_from_file']=str(path.resolve())
        if case:updated.attrs['explorer_signature']=signature
        save_dataset(add_colorimetry(updated),destination,overwrite=args.resume)
    for i,item in enumerate(manifest['calculations'],1):
        source=args.input/item['path']
        assert hashlib.sha256(source.read_bytes()).hexdigest()==item['sha256'],source
        refine(source,args.output/item['path'],case=True)
        print(f'{i}/64 refined single scatter: {item["path"]}',flush=True)
    old_cache.clear();new_cache.clear()
    build_explorer(args.output,cache_dir=cache,allow_download=False,resume=True,progress=print)
    for name in ('urban','rural','dry','humid'):
        for source in sorted((args.input/'sensitivity'/name/'runs').glob('*.nc')):
            relative=source.relative_to(args.input)
            refine(source,args.output/relative)
            print(f'Refined sensitivity component: {relative}',flush=True)
    for source,name in [('docs/fog-and-dust.md','README.md'),('docs/numerical-refinement.md','numerical-refinement.md')]:
        shutil.copy2(source,args.output/name)
    print('Component replacement complete; fresh full-solver regression and validation remain required.',flush=True)

if __name__=='__main__':main()
