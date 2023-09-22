import matplotlib as mpl

from cafeplotter.cafeparser import CafeParser
from cafeplotter.treeplotter import TreePlotter

__version__ = "0.2.0"

__all__ = [
    "CafeParser",
    "TreePlotter",
]


# Setting matplotlib rc(runtime configuration) parameters
# https://matplotlib.org/stable/tutorials/introductory/customizing.html
mpl_rc_params = {
    # SVG
    "svg.fonttype": "none",
}
mpl.rcParams.update(mpl_rc_params)
