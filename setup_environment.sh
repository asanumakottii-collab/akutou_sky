#!/bin/sh
# Keep packages outside cloud-synced Documents folders. Does not change system Python.
set -eu
science_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
science_cache=${XDG_CACHE_HOME:-"$HOME/Library/Caches"}/akutou-science
science_env="$science_cache/python313"
science_python=${AKUTOU_PYTHON:-python3.13}
if [ ! -x "$science_env/bin/python" ]; then
    "$science_python" -m venv "$science_env"
fi
"$science_env/bin/python" -m pip install -r "$science_dir/requirements.txt"
"$science_env/bin/python" -m pip install --no-deps -e "$science_dir[plot]"
printf '%s\n' "$science_env/bin/python"
