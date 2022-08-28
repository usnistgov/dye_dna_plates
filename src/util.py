concentration_label = r"$C_i\;\left[\dfrac{\mu\mathrm{mol}}{\mathrm{L}}\right]$" 
import os
BASE_DIR = os.path.abspath(os.path.dirname(__file__))


def figure_name_to_abspath(fname: str) -> str:
    """Figure name to absolute path

    Parameters
    ----------
    fname : str
        name of figure

    Returns
    -------
    str
        absolute path to name of figure 
    """
    return os.path.join(BASE_DIR, "..", "out", fname)