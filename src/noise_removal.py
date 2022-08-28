import numpy as np


def compute_M_LS(F, C):
    """Calculate :math:`\\mathbf{M}` by least-squares approximation 

    Returns
    -------
    np.array
        :math:`\\mathbf{M}^\\mathrm{LS}`, see Equation (21)

    """
    return F@C/np.inner(C, C)


def compute_c_plus(F, C, M_minus, rho_squared):
    """Compute updated guess of concentrations, :math:`\\mathbf{c}_+`

    Parameters
    ----------
    F : np.ndarray
        Fluorescence :math:`\\mathbf{F}`
    C : np.array
        Dye Concentration :math:`\\mathbf{C}`
    M_minus : np.array
        Guess for :math:`\\mathbf{M}`, :math:`\\mathbf{M}_{-}`
    rho_squared : float
        Weight, :math:`\\rho^2`

    Returns
    -------
    np.ndarray
        :math:`\\mathbf{c}_+` by Equation (S3a)
    """
    return (C*rho_squared + F.T@M_minus)/(rho_squared + np.inner(M_minus, M_minus))


def compute_M_plus(F, c_plus):
    """Get updated guess for M

    Parameters
    ----------
    F : np.ndarray
        Fluorescence matrix :math:`\\mathbf{F}`
    c_plus : np.array
        Concentration matrix updated :math:`\\mathbf{c}_{+}`

    Returns
    -------
    np.array 
        :math:`\\mathbf{M}_{+}` by Equation (S3b)
    """
    return F@c_plus/np.inner(c_plus, c_plus)


def predictor_corrector(F, C, rho_squared, maxiter=100000, print_iter=True):
    """Solve Equation (22) with predictor--corrector algorithm

    Parameters
    ----------
    F : np.ndarray
        Fluorescence data :math:`\\mathbf{F}`
    C : np.ndarray
        Dye concentration data :math:`\\mathbf{C}`
    rho_squared : float
        Weighting factor for concentrations, :math:`\\rho^2` in Equation (22a)
    maxiter : int, optional
        maximum iterations allowed, by default 100000
    print_iter : bool, optional
        whether or not to print the total number of iterations performed, by default True

    Returns
    -------
    (M, c) : tuple(np.array, np.array) 
        solution, :math:`(\\mathbf{M}^\\mathrm{TLS},\\widehat{\\mathbf{C}})`
    
    """
    m, n = F.shape
    if not C.shape == (n,):
        C = C.reshape(n)
    
    M_minus = compute_M_plus(F, C)

    k = 0
    while (k == 0 or dM_2 >= 1e-8) and k < maxiter:
        c_plus = compute_c_plus(F, C, M_minus, rho_squared)
        M_plus = compute_M_plus(F, c_plus)

        error = np.abs(M_plus - M_minus)
        dM_2 = np.sqrt(np.dot(error, error))

        M_minus = M_plus.copy()
        k = k+1
    
    if k == maxiter:
        print("Did not converge!, exiting!")
        quit()
    if print_iter:
        print("Total number of iterations was %i" % k)
    
    return M_plus, c_plus

