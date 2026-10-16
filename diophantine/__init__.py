from .core import solve
from .parser import ParseError, parse_equation
from .result import Result, ParamFamily, RecursiveFamily, GeneratorFamily
from . import numtheory

__version__ = "1.0.0"

__all__ = ["solve", "ParseError", "parse_equation", "Result", "ParamFamily",
           "RecursiveFamily", "GeneratorFamily", "numtheory", "__version__"]
