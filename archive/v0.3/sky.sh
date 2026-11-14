#!/bin/sh
set -eu
science_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
science_cache=${XDG_CACHE_HOME:-"$HOME/Library/Caches"}/akutou-science
science_python=${AKUTOU_PYTHON:-"$science_cache/python313/bin/python"}
if [ ! -x "$science_python" ]; then
    printf '%s\n' 'First run: sh setup_environment.sh' >&2
    exit 1
fi
export PYTHONPYCACHEPREFIX="$science_cache/pycache"
exec "$science_python" -m akutou_sky simulate \
    --output-root "$science_dir/results" --tag improved \
    --cache-dir "$science_dir/model_cache" "$@"
