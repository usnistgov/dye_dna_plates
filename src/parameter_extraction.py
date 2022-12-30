import numpy as np
from src.get_data import CombinedData, C_REF
R_kJ_molK = 8.314 / 1000.


def calculate_relative_brightness(f_SS, f_DS):
    """Calculate relative brightness, Equation (28).

    Parameters
    ----------
    f_SS : np.array
        Molar fluorescence of single-stranded DNA, :math:`\\mathbf{f}^\\mathrm{SS}`.
    f_DS : np.array
        Molar fluorescence of double-stranded DNA, :math:`\\mathbf{f}^\\mathrm{DS}`.

    Returns
    -------
    np.array
        Relative brightness, :math:`\\mathbf{f}^\\mathrm{DS}_j/\\mathbf{f}^\\mathrm{SS}_j` for each :math:`j` associated with SS.

    """
    m = f_SS.shape[0]
    return f_DS[-m:]/f_SS


def calculate_relative_brightness_err(SS_M1, SS_M2, DS_M1, DS_M2, SS_V_M1, SS_V_M2, DS_V_M1, DS_V_M2) -> np.array:
    """Estimate error in relative brightness, Equation (28).

    Parameters
    ----------
    SS_M1 : np.array
        :math:`\\mathbf{M}_1^\\mathrm{SS}`
    SS_M2 : np.array
        :math:`\\mathbf{M}_2^\\mathrm{SS}`
    DS_M1 : np.array
        :math:`\\mathbf{M}_1^\\mathrm{DS}`
    DS_M2 : np.array
        :math:`\\mathbf{M}_2^\\mathrm{DS}`
    SS_V_M1 : np.array
        :math:`V(\\mathbf{M}_1^\\mathrm{SS})`
        _description_
    SS_V_M2 : np.array
        :math:`V(\\mathbf{M}_2^\\mathrm{SS})`
    DS_V_M1 : np.array
        :math:`V(\\mathbf{M}_1^\\mathrm{DS})`
    DS_V_M2 : np.array
        :math:`V(\\mathbf{M}_2^\\mathrm{DS})`

    Returns
    -------
    np.array
        Estimate of error in relative brightness
    """
    m = SS_M1.shape[0]
    M1D, M2D = DS_M1[-m:], DS_M2[-m:]
    V1D, V2D = DS_V_M1[-m:], DS_V_M2[-m:]

    Va = 2*SS_V_M1 + SS_V_M2
    Vb = (V1D + V2D + (M1D + M2D)*(M1D + M2D))*(2*V1D + V2D + (2*M1D - M2D)*(2*M1D - M2D)) \
        - (M1D + M2D)*(M1D + M2D)*(2*M1D - M2D)*(2*M1D - M2D) \
        + (2*V1D + V2D + (2*M1D + M2D)*(2*M1D + M2D))*(V1D + V2D + (M1D - M2D)*(M1D - M2D)) \
        - (2*M1D + M2D)*(2*M1D + M2D)*(M1D - M2D)*(M1D - M2D)
    Ea = (2*SS_M1 - SS_M2)
    Eb = (M1D + M2D)*(2*M1D - M2D) - (2*M1D + M2D)*(M1D - M2D)
    E_num = Ea*Eb
    V_numerator = (Va + Ea*Ea)*(Vb + Eb*Eb) - Ea*Ea*Eb*Eb
    Va = 2*V1D + V2D
    Vb = (SS_V_M1 + SS_V_M2 + (SS_M1 + SS_M2)*(SS_M1 + SS_M2))*(2*SS_V_M1 + SS_V_M2 + (2*SS_M1 - SS_M2)*(2*SS_M1 - SS_M2)) \
        - (SS_M1 + SS_M2)*(SS_M1 + SS_M2)*(2*SS_M1 - SS_M2)*(2*SS_M1 - SS_M2) \
        + (2*SS_V_M1 + SS_V_M2 + (2*SS_M1 + SS_M2)*(2*SS_M1 + SS_M2))*(SS_V_M1 + SS_V_M2 + (SS_M1 - SS_M2)*(SS_M1 - SS_M2)) \
        - (2*SS_M1 + SS_M2)*(2*SS_M1 + SS_M2)*(SS_M1 - SS_M2)*(SS_M1 - SS_M2)
    Ea = (2*M1D - M2D)
    Eb = (SS_M1 + SS_M2)*(2*SS_M1 - SS_M2) - (2*SS_M1 + SS_M2)*(SS_M1 - SS_M2)
    V_denominator = (Va + Ea*Ea)*(Vb + Eb*Eb) - Ea*Ea*Eb*Eb
    E_den = Ea*Eb
    return np.sqrt((E_den*E_den*V_numerator + E_num*E_num*V_denominator)/E_den/E_den/E_den/E_den)


class Parameters:
    """Stores multiple instances of CombinedData for one DNA type
    
    Attributes
    ----------
    N : int
        number of nucleobases, :math:`N`
    M1 : np.array
        :math:`\\mathbf{M}^\\mathrm{TLS}` associated with :math:`\\mathbf{D}=1`
    M2 : np.array
        :math:`\\mathbf{M}^\\mathrm{TLS}` associated with :math:`\\mathbf{D}=2`. 
        Several high temperatures are removed to reflect smaller temperature range associated with :math:`\\mathbf{D}=1`
    C1 : np.array
        :math:`\\widehat{\\mathbf{C}}` associated with :math:`\\mathbf{D}=1`
    C2 : np.array
        :math:`\\widehat{\\mathbf{C}}` associated with :math:`\\mathbf{D}=2`
    r : np.array
        :math:`r`, as defined in Equation (S7).
    V_C1 : np.array
        :math:`V(\\mathbf{C})` associated with :math:`\\mathbf{D}=1`
    V_C2 : np.array
        :math:`V(\\mathbf{C})` associated with :math:`\\mathbf{D}=2`
    V_M1 : np.array
        :math:`V(\\mathbf{M})` associated with :math:`\\mathbf{D}=1`
    V_M2 : np.array
        :math:`V(mathbf{M})` associated with :math:`\\mathbf{D}=2`. 
        Several high temperatures are removed to reflect smaller temperature range associated with :math:`\\mathbf{D}=1`
    dT : float
        Change in temperature from one cycle to next, :math:`\\Delta T`
    
    """
    def __init__(self, cls1: CombinedData, cls2: CombinedData):
        """Initialize data

        Note
        ----
        Since different temperature ranges for each, need to make subset
        of dataset that has more temperatures.
        Dataset with lower DNA concentration :code:`cls1` always has less temperatures.

        Parameters
        ----------
        cls1 : CombinedData
            Data of DNA type at :math:`\\mathbf{D}=1`
        cls2 : CombinedData
            Data of DNA type at :math:`\\mathbf{D}=2`
        """
        assert np.abs(cls1.D - 1.) < 1e-12, "Wrong imput concentration"
        assert np.abs(cls2.D - 2.) < 1e-12, "Wrong imput concentration"

        self.C1 = cls1.C_hat
        self.C2 = cls2.C_hat

        m = cls1.T.shape[0]
        self.M1 = cls1.M_tls
        self.M2 = cls2.M_tls[-m:]
        self.T = cls1.T

        self.V_C1 = cls1.V_C
        self.V_C2 = cls2.V_C
        self.V_M1 = cls1.V_M
        self.V_M2 = cls2.V_M[-m:]

        self.r = self.M1 / self.M2
        self.dT = self.T[1] - self.T[0]

        self.N = cls1.N
        assert self.N == cls2.N, "Cannot combine different DNA lengths"

    def get_phi_1(self) -> np.array:
        """Get :math:`\\varphi_{1}`, vectorized version of Equation (S6a)

        Returns
        -------
        np.array
            
        """
        return 2*self.r - 1
    
    def get_std_phi_1(self) -> np.array:
        """Get estimate of standard deviation in :math:`\\varphi_{1}`

        Returns
        -------
        np.array

            .. math::
                \\frac{2}{\\mathbf{M}_{2j}}\\sqrt{V(\\mathbf{M}_{1j}) + r_j^2V(\\mathbf{M}_{2j})}

        """
        return 2*np.sqrt(
            self.V_M1 + self.r*self.r*self.V_M2
        )/self.M2
    
    def get_phi_2(self):
        """Get :math:`\\varphi_{2}`, vectorized version of Equation (S6b)

        Returns
        -------
        np.array
            
        """
        return 2 - 1/self.r
    
    def get_std_phi_2(self):
        """Get estimate of standard deviation in :math:`\\varphi_{2}`

        Returns
        -------
        np.array

            .. math::
                \\frac{1}{\\mathbf{M}_{1j}}\\sqrt{V(\\mathbf{M}_{2j}) + r_j^{-2}V(\\mathbf{M}_{1j})}

        """
        return np.sqrt(
            self.V_M2 + self.V_M1/self.r/self.r
        )/self.M1
    
    def get_theta_b_all_1(self, n):
        return n*self.C1[-1]*(2*self.r - 1)/self.N

    def get_psi_j1(self, j: int):
        return self.C1*(2*self.r[j] - 1)
    
    def get_std_psi_j1(self, j: int):
        return np.sqrt(
            (2*self.r[j] - 1)*(2*self.r[j] - 1)*self.V_C1
            + 4*self.C1*self.C1/self.M2[j]/self.M2[j]*self.V_M1[j]
            + 4*self.C1*self.C1*self.r[j]*self.r[j]/self.M2[j]/self.M2[j]*self.V_M2[j]
        )
    
    def get_psi_j2(self, j: int):
        return self.C2*(1 - 1/2/self.r[j])
    
    def get_std_psi_j2(self, j: int):
        return np.sqrt(
            (1 - 1/2/self.r[j])*(1 - 1/2/self.r[j])*self.V_C2
            + (self.C2*self.C2/4/self.M1[j]/self.M1[j])*self.V_M2[j]
            + self.C2*self.C2/self.r[j]/self.r[j]/self.M1[j]/self.M1[j]/4*self.V_M1[j]
        )
    
    def get_K(self) -> np.array:
        """Get :math:`\\mathbf{K}` from vectorized version of Equation (24)

        Returns
        -------
        np.array
            :math:`\\mathbf{K}`
        """
        return (2*self.M1 - self.M2)/(self.M2 - self.M1)/2
    
    def get_K_std(self) -> np.array:
        """Get standard deviation estimate of :math:`\\mathbf{K}`

        Returns
        -------
        np.array

            .. math::
                \\sqrt{
                    \\Delta \\mathbf{M}^2\\left(4V(\\mathbf{M}_1)+2V(\\mathbf{M}_2)\\right)
                    + \\frac{\\Delta H^2}{8\\Delta \\mathbf{M}^4}\\left(V(\\mathbf{M}_1)+V(\\mathbf{M}_2)\\right)
                }
            
            where :math:`\\Delta H := 2\\mathbf{M}_1 - \\mathbf{M}_2`, 
            :math:`\\Delta \\mathbf{M} := \\mathbf{M}_2 - \\mathbf{M}_1`
        
        
        """
        DM = self.M2 - self.M1
        DH = 2*self.M1 - self.M2
        return np.sqrt((DM*DM*(4*self.V_M1+ 2*self.V_M2) + DH*DH*(self.V_M1+ self.V_M2))/8/DM/DM/DM/DM)
    
    def get_f(self) -> np.array:
        """Get :math:`\\mathbf{f}` from vectorized version of Equation (27)

        Returns
        -------
        np.array
            :math:`\\mathbf{f}`
        """
        return self.M1 + self.M2 - (2*self.M1 + self.M2)/(2*self.M1 - self.M2)*(self.M1 - self.M2)

    def get_f_std(self) -> np.array:
        """Get standard deviation estimate of :math:`\\mathbf{f}`

        Returns
        -------
        np.array

            .. math::
                \\sqrt{
                    V(\\mathbf{M}_1) + V(\\mathbf{M}_2)
                     + \\frac{E_B^2V_A + E_A^2 V_B}{E_B^4}
                }
            
            where 

            .. math::
                V_A := (V_B + H_+^2)(V(\\mathbf{M}_1) + V(\\mathbf{M}_2) + \\Delta \\mathbf{M}^2) - E_A^2

                V_B := 2V(\\mathbf{M}_1) + V(\\mathbf{M}_2)

                E_A := \\Delta \\mathbf{M} H_{+}

                E_B := 2\\mathbf{M}_1 - \\mathbf{M}_2

                \\Delta \\mathbf{M} := \\mathbf{M}_2 - \\mathbf{M}_1

                H_+ := 2\\mathbf{M}_1 + \\mathbf{M}_2
            
        """
        HH =  2*self.M1 + self.M2
        DM = self.M2 - self.M1
        E_B = 2*self.M1 - self.M2 
        V_B = 2*self.V_M1 + self.V_M2
        E_A = DM*HH
        V_A = (V_B + HH*HH)*(self.V_M1 + self.V_M2 + DM*DM) - E_A*E_A

        return np.sqrt(
            self.V_M1 + self.V_M2 + (E_B*E_B*V_A + E_A*E_A*V_B)/(E_B*E_B*E_B*E_B)
        )
    
    def get_dg(self) -> np.array:
        """Get free energy of dye binding, :math:`\\Delta g`, vectorized version of Equation (29).

        Returns
        -------
        np.array

        """
        K = self.get_K()
        return -R_kJ_molK * self.T * np.log(K/C_REF/self.N)
    
    def get_dh(self) -> np.array:
        """Get differential enthalpy of binding, :math:`\\Delta h` as the vectorized version of Equation (30).

        Returns
        -------
        np.array

        """
        K = self.get_K()
        return R_kJ_molK * self.T[2:]*self.T[:-2] / 2 / self.dT * np.log(K[2:] / K[:-2])
    
    def get_dg_std(self) -> np.array:
        """Get estimate of standard deviation in :math:`\\Delta g`.

        Returns
        -------
        np.array

            .. math::
                \\frac{RT_j}{2\\Delta \\mathbf{M}_j\\left(2\\mathbf{M}_{j1}-\\mathbf{M}_{j2}\\right)} \\sqrt{
                    \\mathbf{M}_{j2}^2 V(\\mathbf{M}_{j1})
                    +
                    \\mathbf{M}_{j1}^2 V(\\mathbf{M}_{j2})
                }
            
            where :math:`\\Delta \\mathbf{M}_j = \\mathbf{M}_{j2} - \\mathbf{M}_{j1}`.

        """
        return R_kJ_molK * self.T / 2. /(2.*self.M1 - self.M2)/(self.M2 - self.M1) * np.sqrt(
            self.M2*self.M2 * self.V_M1
            +
            self.M1*self.M1*self.V_M2
        )
    
    def get_dh_std(self) -> np.array:
        """Get estimate of error in :math:`\\Delta h_j`

        Returns
        -------
        np.array

        """
        K = self.get_K()

        DM = self.M2 - self.M1
        # partial derivatives, ignoring constant prefactor
        dhj__M1jp1 = self.M2[2:]/2./K[2:]/DM[2:]/DM[2:]
        dhj__M2jp1 = self.M1[2:]/2./K[2:]/DM[2:]/DM[2:]
        dhj__M1jm1 = self.M2[:-2]/2./K[:-2]/DM[:-2]/DM[:-2]
        dhj__M2jm1 = self.M1[:-2]/2./K[:-2]/DM[:-2]/DM[:-2]

        return R_kJ_molK * self.T[2:]*self.T[:-2] / 2 / self.dT * np.sqrt(
            dhj__M1jp1*dhj__M1jp1*self.V_M1[2:]
            + dhj__M2jp1*dhj__M2jp1*self.V_M2[2:]
            + dhj__M1jm1*dhj__M1jm1*self.V_M1[:-2]
            + dhj__M2jm1*dhj__M2jm1*self.V_M2[:-2]
        )
        