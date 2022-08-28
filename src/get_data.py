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
    Uses Equation (12) of manuscript 
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
        CSV file formatted like a 96-well plate. The top left corner looks like

        +-----+---+---+
        |  Row| 1 | 2 |
        +=====+===+===+
        | A   | . | . |
        +-----+---+---+
        | B   | . | . |
        +-----+---+---+
        | C   | . | . |
        +-----+---+---+

        The values are concentrations of dye in units of mol/L



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
        Fluorescence data, :math:`F_k^{t\\ell}` (see Equation 14 of main text)
    C : np.ndarray
        Dye concentrations, :math:`C` (see Equation 13 of main text)
    D : float
        DNA concentration, :math:`D_k` in mol/L
    t : str
        Type of DNA, :math:`t`, :code:`"SS"` or :code:`"DS"` or :code:`"None"`.
    l : str
        Replicate name, :math:`\\ell` is A, B, or C

    """
    def __init__(self, fluorescence_file_name, D_k, t, l,
                       dye_conc_file_name="dye_conc_uM.csv"):
        """Scale the data before interpolating/solving optimization problem.

        Parameters
        ----------
        fluorescence_file_name : str
            name of fluorescence file within data folder
        D_k : float
            Total concentration of DNA in mol/L :math:`D_k`
        dye_conc_file_name : str, optional
            name of dye concentration file name within data folder, defaults to :code:`"dye_conc_uM.csv"`
        t : str
            Type of DNA, :code:`"SS"`, :code:`"DS"`, or :code:`"None"`.
        l : str
            Replicate name, :math:`\\ell` is A, B, or C

        """
        df = excel_to_data(os.path.join(PATH_TO_DATA, fluorescence_file_name))
        c_total = get_C(file_name=os.path.join(PATH_TO_DATA, dye_conc_file_name))

        self.T = df.index.to_numpy()
        self.C = np.array([c_total[key] for key in df.columns])
        self.F = df.to_numpy()
        self.D = D_k

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
        Fluorescence data, :math:`\\mathbf{F}_k^t` in Equation (16a)
    C : np.ndarray
        Dye concentrations, :math:`\\mathbf{C}_k^t` in Equation (16a)
    D : float
        DNA concentration, :math:`\\mathbf{D}_k` in Equation (16a)
    t : str
        Type of DNA, :math:`t`, :code:`"SS"` or :code:`"DS"`
    M_tls : np.array
        :math:`\\mathbf{M}^\\mathrm{TLS}`, set externally, defaults to np.array([])
    C_hat : np.array
        :math:`\\widehat{\\mathbf{C}}`, set externally, defaults to np.array([])
    V_M : np.array
        :math:`V(\\mathbf{M})` (see Section S1.2), defaults to np.array([])
    V_C : np.array
        :math:`V(\\mathbf{C})` (see Section S1.2), defaults to np.array([])
    M_std : np.array
        :math:`\\sqrt{V(\\mathbf{M})}` (see Section S1.2), defaults to np.array([])
    C_std : np.array
        :math:`\\sqrt{V(\\mathbf{C})}` (see Section S1.2), defaults to np.array([])


    """
    def __init__(self, *replicates: typing.Tuple[RawData]) -> None:
        """Combine replicate :math:`F` and :math:`C`

        Parameters
        ----------
        replicates : typing.List[Data]
            list of plates to gather together as replicates
        """
        t = None
        T = None
        D = None
        F = []
        C = []
        for cls in replicates:
            if t is not None:
                assert np.max(np.abs(T - cls.T)) < 1e-14, "Inconsistent temperatures"
                assert t == cls.t, "Cannot combine different DNA types"
                assert np.abs(D - cls.D) < 1e-14, "Cannot combine different DNA concentrations"
            else:
                t = cls.t
                D = cls.D
                T = cls.T[:]

            F.append(cls.F)
            C.append(cls.C)
            
        self.F = np.hstack(F) / F_REF
        self.C = np.hstack(C) / C_REF
        self.t = t
        self.D = D / C_REF
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
        while (np.mean(self.F[i_T, :]) < F_min):
            i_T += 1
        
        self.F = self.F[i_T:, :]
        self.T = self.T[i_T:]



pluto_plateA_rot = RawData(fluorescence_file_name="6-2022_intercalate-pluto-inter_PlateA_Rot_Pluto_6-14-22_data.xls",
                        D_k=2e-6, t="SS", l="A",
                        dye_conc_file_name="dye_conc_uM_rotated.csv")