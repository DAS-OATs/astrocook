sphinx-apidoc -f -e --no-toc -o \
    docs/api astrocook \
    astrocook/core/atomic_data.py \
    astrocook/core/photometry.py \
    astrocook/core/session_manager.py \
    astrocook/core/spectrum_operations.py \
    astrocook/core/system_list_migration.py \
    astrocook/core/utils.py \
    astrocook/gui/*.py \
    astrocook/io/*.py

rm docs/api/astrocook.rst docs/api/astrocook.settings.rst docs/api/modules.rst 2> /dev/null