try:
    from importlib.metadata import version as _version
    __version__ = _version("zamanda-yolculuk-zt005-2")
except Exception:  # source checkout without pip install -e .
    __version__ = "14.0.0"
