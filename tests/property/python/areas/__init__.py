"""Property-test area modules for the SageMath runner.

Each ``<area>.py`` in this package exports a ``FUNCTIONS`` dict mapping the
``function`` names used in ``tests/property/cases/<area>.cases.json`` to
callables.  ``runner.py`` discovers them by globbing this directory, so adding
an area never requires editing a shared file.

Modules whose name starts with ``_`` are helpers, not areas, and are skipped by
discovery.
"""
