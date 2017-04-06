"""
=====================================
Neural decoding (:mod:`fklab.decode`)
=====================================

.. currentmodule:: fklab.decode

Tools for neural decoding.

.. automodule:: fklab.decode.compressed_kde

"""

from .compressed_kde import *

__all__ = [_s for _s in dir() if not _s.startswith('_')]
