import pandas as pd
import typing
import os
import numpy as np
from matplotlib.colors import Normalize, Colormap
from src.wells import ORDERED_WELLS

BASE_DIR = os.path.abspath(os.path.dirname(__file__))
PATH_TO_DATA = os.path.join(BASE_DIR, "..", "data")
cmap_temperature = Colormap("coolwarm")
cmap_temperature_norm = Normalize(vmin=283., vmax=353.)
F_REF = 1e6
C_REF = 1e-6


def excel_to_data(f_name: str, channel="GREEN"):
    """
    Convert "Raw Data" sheet of excel file to pandas dataframe.
    Uses Equation (9) of manuscript
    to calculate temperature associated with each cycle.

    Parameters
    ----------
    f_name : str
        name of excel file
    channel : str, optional
        name of channel to investigate, defaults to :code:`"GREEN"`

    Returns
    -------
    pd.DataFrame
        Formatted data frame, with wells sorted from A1, A2 ... H11, H12

    """
    sheet = pd.ExcelFile(f_name)
    df = sheet.parse(sheet_name="Raw Data", header=7)
    df["T [K]"] = 353. - 0.5*df["Cycle"]
    df = df.pivot(index="T [K]", columns=['Well'], values=channel).sort_index(ascending=False)
    return df[ORDERED_WELLS]


def get_C(file_name):
    """
    Get total dye concentration 
    associated with each well.

    Parameters
    ----------
    file_name : str
        CSV file containing values of
        concentrations of dye in units of :math:`\\mu` mol/L.
        Formatted like a 96-well plate:

        +-----+---+---+---+---+---+---+---+---+---+---+---+---+
        |  Row| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11| 12|
        +=====+===+===+===+===+===+===+===+===+===+===+===+===+
        | A   | . | . | . | . | . | . | . | . | . | . | . | . |
        +-----+---+---+---+---+---+---+---+---+---+---+---+---+
        | B   | . | . | . | . | . | . | . | . | . | . | . | . |
        +-----+---+---+---+---+---+---+---+---+---+---+---+---+
        | C   | . | . | . | . | . | . | . | . | . | . | . | . |
        +-----+---+---+---+---+---+---+---+---+---+---+---+---+
        | D   | . | . | . | . | . | . | . | . | . | . | . | . |
        +-----+---+---+---+---+---+---+---+---+---+---+---+---+
        | E   | . | . | . | . | . | . | . | . | . | . | . | . |
        +-----+---+---+---+---+---+---+---+---+---+---+---+---+
        | F   | . | . | . | . | . | . | . | . | . | . | . | . |
        +-----+---+---+---+---+---+---+---+---+---+---+---+---+
        | G   | . | . | . | . | . | . | . | . | . | . | . | . |
        +-----+---+---+---+---+---+---+---+---+---+---+---+---+
        | H   | . | . | . | . | . | . | . | . | . | . | . | . |
        +-----+---+---+---+---+---+---+---+---+---+---+---+---+


    Returns
    -------
    Dictionary
        Mapping of well name ("A1",...) to dye concentration [units of mol/L]

    """
    df = pd.read_csv(file_name)
    c = {}
    for col in map(str, range(1, 13)):
        for i, row in enumerate(df["Row"]):
            c[row + col] = df[col][i] * 1e-6

    return c


class RawData:
    """Stores raw data.

    Attributes
    ----------
    F : np.ndarray
        Fluorescence data, :math:`F_d^{t,l}` (see Eq. 11)
    C : np.ndarray
        Dye concentrations, :math:`C` (see Eq. 10)
    B : float
        DNA concentration, :math:`B_d` in mol/L
    t : str
        Type of DNA, :math:`t`, :code:`"SS"` or :code:`"DS"` or :code:`"None"`.
    l : str
        Replicate name, :math:`l` is :code:`"A"`, :code:`"B"`, or code:`"C"`.
    N : int
        number of nucleobases per single strand

    """
    def __init__(self, fluorescence_file_name, B_d, t, l, N,
                       dye_conc_file_name="dye_conc_uM.csv"):
        """Scale the data before interpolating/solving optimization problem.

        Parameters
        ----------
        fluorescence_file_name : str
            name of fluorescence file within data folder
        B_d : float
            Total concentration of nucleobases in mol/L :math:`B_d`
        dye_conc_file_name : str, optional
            name of dye concentration file name within data folder, defaults to :code:`"dye_conc_uM.csv"`
        t : str
            Type of DNA, :code:`"SS"`, :code:`"DS"`, or :code:`"None"`.
        l : str
            Replicate name, :math:`\\ell` is A, B, or C
        N : int
            number of nucleobases per single strand

        """
        df = excel_to_data(os.path.join(PATH_TO_DATA, fluorescence_file_name))
        c_total = get_C(file_name=os.path.join(PATH_TO_DATA, dye_conc_file_name))

        self.T = df.index.to_numpy()
        self.C = np.array([c_total[key] for key in df.columns])
        self.F = df.to_numpy()
        self.B = B_d
        self.N = N

        self.t = t
        assert ((t == "SS") or (t == "DS") or (t == "None")), "Incorrect input DNA type"
        self.l = l
        assert (l in list("ABC")), "Incorrect replicate name"


class CombinedData:
    """Combined data from several replicate plates. 
    See changes to :math:`F` and :math:`C` in 
    *Raw Data and Scaling* portion of *Results and Discussion.*

    Attributes
    ----------
    F : np.ndarray
        Fluorescence data, :math:`\\mathbf{F}_{\\mathbf{D}}^t` in Equation (13a)
    C : np.ndarray
        Dye concentrations, :math:`\\mathbf{C}_{\\mathbf{D}}^t` in Equation (13a)
    D : float
        DNA concentration, :math:`\\mathbf{D}` in Equation (13a)
    t : str
        Type of DNA, :math:`t`, :code:`"SS"` or :code:`"DS"`
    N : int
        number of nucleobases per single strand, :math:`N`
    M_tls : np.array
        :math:`\\mathbf{M}^\\mathrm{TLS}` in Equation (18a), set externally, defaults to np.array([])
    C_hat : np.array
        :math:`\\widehat{\\mathbf{C}}` in Equation (18a), set externally, defaults to np.array([])
    V_M : np.array
        :math:`V(\\mathbf{M})` (see section S1.3 in supporting material), defaults to np.array([])
    V_C : np.array
        :math:`V(\\mathbf{C})` (see section S1.3 in supporting material), defaults to np.array([])
    M_std : np.array
        :math:`\\sqrt{V(\\mathbf{M})}`, defaults to np.array([])
    C_std : np.array
        :math:`\\sqrt{V(\\mathbf{C})}`, defaults to np.array([])


    """
    def __init__(self, *replicates: typing.Tuple[RawData]) -> None:
        """Combine replicate :math:`F_d^{t,l}` and :math:`C`

        Parameters
        ----------
        replicates : typing.List[Data]
            list of plates to gather together as replicates
        """
        t = None
        T = None
        B = None
        N = None
        F = []
        C = []
        for cls in replicates:
            if t is not None:
                assert np.max(np.abs(T - cls.T)) < 1e-14, "Inconsistent temperatures"
                assert t == cls.t, "Cannot combine different DNA types"
                assert N == cls.N, "Cannot combine with different DNA lengths"
                assert np.abs(B - cls.B) < 1e-14, "Cannot combine different DNA concentrations"
            else:
                t = cls.t
                B = cls.B
                T = cls.T[:]
                N = cls.N

            F.append(cls.F)
            C.append(cls.C)
            
        self.F = np.hstack(F) / F_REF
        self.C = np.hstack(C) / C_REF
        self.t = t
        self.N = N
        self.D = B / C_REF / self.N
        self.T = T

        # default params to be changed upon solution of (22)
        self.M_tls = np.array([])
        self.C_hat = np.array([])
        self.V_C = np.array([])
        self.V_M = np.array([])
        self.M_std = np.array([])
        self.C_std = np.array([])
    
    def make_subset(self, F_min):
        """Make subset of data as described in *Raw Data and Scaling* portion of manuscript
        
        Note
        ----
        :math:`\\mathbf{F}` and :math:`\\mathbf{C}` are overwritten.
        
        """
        wells = self.C*C_REF > 0.5e-6

        self.F = self.F[:, wells]
        self.C = self.C[wells]

        if self.t == 'DS':
            wells = self.C*C_REF <= 2.5e-6
            self.F = self.F[:, wells]
            self.C = self.C[wells]
        
        i_T = 0
        while np.mean(self.F[i_T, :]) < F_min:
            i_T += 1
        
        self.F = self.F[i_T:, :]
        self.T = self.T[i_T:]

