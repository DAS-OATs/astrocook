from .message import *
from .vars import *
import ast
from astropy import constants as ac
from astropy.table import Table
from copy import deepcopy as dc
import cProfile
#from decorator import decorator
import json
import logging
from matplotlib import pyplot as plt
import numpy as np
import pstats
from scipy.ndimage import median_filter
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
from scipy.special import wofz
#from lmfit.lineshapes import gaussian as gauss

prefix = 'functions'

def _gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def _fadd(a, u, deriv=False):
    """ @brief Real part of the Faddeeva function Re(F)
    @param a First abstract variable
    @param u Second abstrac variable
    @return Re(F(a, u))
    """
    D = wofz(u + 1j * a)
    F = np.real(D)
    if deriv:
        dF_da = 2 * (F*a + np.imag(D)*u - 1/np.sqrt(np.pi))
        dF_du = 2 * (np.imag(D)*a - F*u)
        return F, (dF_da, dF_du)
    else:
        return F

def _voigt_par_convert(x, z, N, b, btur, trans):
    if trans == 'unknown':
        xem = z*au.nm
        xobs = z*au.nm
    else:
        xem = xem_d[trans]
        xobs = xem*(1+z)
    fosc = fosc_d[trans]
    gamma = gamma_d[trans]/au.s
    b_qs = np.sqrt(b**2 + btur**2)
    atom = fosc * ac.e.esu**2 / (ac.m_e * ac.c)
    tau0 = np.sqrt(np.pi) * atom * N * xem / b_qs
    a = 0.25 * gamma * xem / (np.pi * b_qs)
    u = ac.c/b_qs * ((x/xobs).to(au.dimensionless_unscaled) - 1)
    return tau0, a, u

def _voigt_par_convert_new(x, z, N, b, btur, trans, deriv=False):
    if trans == 'unknown':
        xem = z #*au.nm
        xobs = z #*au.nm
    else:
        xem = xem_d[trans].value
        xobs = xem*(1+z)
    fosc = fosc_d[trans]
    gamma = gamma_d[trans] #/au.s
    b_qs = np.sqrt(b**2 + btur**2)
    atom = fosc *  844.7972564303736 #* au.Fr**2 * au.s / (au.kg * au.m)

    tau0 = np.sqrt(np.pi) * atom * N * xem / b_qs
    a = 0.25 * gamma * xem / (np.pi * b_qs)
    u = 299792458/b_qs * (x/xobs - 1)

    #print('before', a)
    if deriv:
        dtau0_dlogN = np.sqrt(np.pi) * atom * xem / b_qs * N * np.log(10)
        dtau0_db = -np.sqrt(np.pi) * atom * N * xem/(b_qs**3) * b
        dtau0_dbtur = -np.sqrt(np.pi) * atom * N * xem/(b_qs**3) * btur

        da_db = 0.25 * gamma * xem / (np.pi * b_qs**3) * b
        da_dbtur = 0.25 * gamma * xem / (np.pi * b_qs**3) * btur

        du_dz = -299792458/b_qs * x*xem/(xobs**2)
        du_db = -299792458/(b_qs**3) * b * (x/xobs - 1)
        du_dbtur = -299792458/(b_qs**3) * btur * (x/xobs - 1)


    tau0_f = 1e-17
    a_f = 1e-12
    u_f = 1e-3

    tau0_f_mod = 1.35e-1
    a_f_mod = 2e-8
    u_f_mod = 6e-2

    tau0 = tau0 * tau0_f
    #print(N, tau0)
    #print(N+10, np.sqrt(np.pi) * atom * (N+10) * xem / b_qs * tau0_f)
    #print((np.sqrt(np.pi) * atom * (N+10) * xem / b_qs * tau0_f-tau0)/10)
    a = a * a_f
    u = u * u_f
    #print(a)
    if deriv:
        dtau0_dlogN = dtau0_dlogN * tau0_f #* 60
        dtau0_db = dtau0_db * tau0_f #* 6e1
        dtau0_dbtur = dtau0_dbtur * tau0_f

        da_db = da_db * a_f #* 2e4
        da_dbtur = da_dbtur * a_f

        du_dz = du_dz * u_f * 6e1
        du_db = du_db * u_f * 5.5e1
        du_dbtur = du_dbtur * u_f

        return tau0, a, u, \
            (dtau0_dlogN, dtau0_db, dtau0_dbtur), \
            (da_db, da_dbtur),  (du_dz, du_db, du_dbtur)
    else:
        return tau0, a, u


def lines_voigt_jac(x0, x, series='CIV', resol=70000, spec=None, apply_bounds_transformation=True):
    for i in range(0, len(x0), 3):
        z, logN, b = x0[i], x0[i+1], x0[i+2]
        btur = 0
        x = x * au.nm
        z = z * au.dimensionless_unscaled
        N = 10**logN / au.cm**2
        b = b * au.km/au.s
        btur = btur * au.km/au.s

        dI_dz = np.zeros(len(x))
        dI_dlogN = np.zeros(len(x))
        #dI_dlogN_new = np.zeros(len(x))
        dI_db = np.zeros(len(x))
        dI_dbtur = np.zeros(len(x))
        #print('jac', z, N, b)

        #model = lines_voigt(x, z, logN, b, btur, series)
        for t in trans_parse(series):
            tau0, a, u, \
                (dtau0_dlogN, dtau0_db, dtau0_dbtur), \
                (da_db, da_dbtur),  (du_dz, du_db, du_dbtur) \
                = _voigt_par_convert_new(x.value, z.value, N.value, b.value,
                                         btur.value, t, deriv=True)
            F, (dF_da, dF_du) = _fadd(a, u, deriv=True)

            dI_dtau0 = -F * np.exp(-tau0 * F)
            #print('jac', a)
            dI_dF = -tau0 * np.exp(-tau0 * F)

            dI_dz += dI_dF*dF_du*du_dz
            dI_dlogN += dI_dtau0*dtau0_dlogN
            #dI_dlogN_new += F*dtau0_dlogN
            dI_db += dI_dtau0*dtau0_db + dI_dF*dF_da*da_db + dI_dF*dF_du*du_db
            dI_dbtur += dI_dtau0*dtau0_dbtur + dI_dF*dF_da*da_dbtur + dI_dF*dF_du*du_dbtur

        #dI_dlogN_new = -model*dI_dlogN_new
        dI_dz = convolve_simple(dI_dz, psf_gauss(x.value, resol, spec))
        dI_dlogN = convolve_simple(dI_dlogN, psf_gauss(x.value, resol, spec))
        dI_db = convolve_simple(dI_db, psf_gauss(x.value, resol, spec))

        """
    print(np.array([dI_dz, dI_dlogN, dI_db]).T)
    return np.array([dI_dz, dI_dlogN, dI_db]).T#, dI_dbtur])
        """
        if i==0:
            J = np.array([dI_dz, dI_dlogN, dI_db])
        else:
            #print(J)
            J = np.append(J, [dI_dz], axis=0)
            J = np.append(J, [dI_dlogN], axis=0)
            J = np.append(J, [dI_db], axis=0)
            #print(J)
    #print(J)
    #J = np.reshape(J, (len(J)//len(x0), len(x0)))
    #print(J)
    return J.T
    #"""

def lines_voigt(x, z, logN, b, btur, series='Ly_a'):
    """ @brief Voigt function (real part of the Faddeeva function, after a
    change of variables)

    @param x Wavelength domain (in nm)
    @param z Redshift
    @param N Column density (in cm^-2)
    @param b Doppler broadening (in km s^-1)
    @param btur Turbulent broadening (in km s^-1)
    @param series Series of ionic transition
    @param xem Wavelength of the line (in nm)
    @param tab Table with the Faddeeva function
    @return Voigt function over x
    """

    #x = x[0] * au.nm
    x = x * au.nm
    z = z * au.dimensionless_unscaled
    N = 10**logN / au.cm**2
    b = b * au.km/au.s
    btur = btur * au.km/au.s
    model = np.ones(np.size(np.array(x)))
    #print('voigt', z, N, b)
    for t in trans_parse(series):
        tau0, a, u = _voigt_par_convert_new(x.value, z.value, N.value, b.value,
                                            btur.value, t)
        dtau0, da, du = _voigt_par_convert_new(x.value, z.value, N.value+10, b.value,
                                              btur.value, t)
        F = _fadd(a, u)
        dF = _fadd(da, du)
        #print(tau0, F, np.exp(-tau0 * F))
        #print(dtau0, dF, np.exp(-dtau0 * dF))
        #print((np.exp(-dtau0 * dF)-np.exp(-tau0 * F))/10)
        #print('voigt', a)
        model *= np.array(np.exp(-tau0 * F))

    return model

def zero(x):
    return 0*x


def adj_gauss(x, z, ampl, sigma, series='Ly_a'):
    model = np.ones(len(x))
    #for t in series_d[series]:
    for t in trans_parse(series):
        c = (1+z)*xem_d[t].to(au.nm).value
        model += ampl*np.exp(-(0.5 * (x-c) / sigma)**2)
    return model

def convolve(data, psf):
    s = 0
    l = 0
    for i, k in enumerate(psf):
        #print(i, k)
        s += l
        l = len(k)
        k_arr = k[np.where(k>0)]
        k_arr = k_arr/np.sum(k_arr)
        #data_arr = data
        data_arr = data[s:s+l]
        pad_l = len(k_arr)#*2
        pad = np.ones(pad_l)

        temp_arr = np.concatenate((pad*data_arr[0], data_arr, pad*data_arr[-1]))
        conv = np.convolve(temp_arr, k_arr, mode='valid')[pad_l//2+1:][:l]
        #print(len(temp_arr), len(conv))
        #plt.plot(range(len(conv)),temp_arr[1:-1], c='red', alpha=0.5)
        #plt.plot(range(len(conv)),conv, c='black', alpha=0.5)
        #plt.show()
        if i == 0:
            ret = conv
        else:
            ret = np.append(ret, conv)
    #print([len(p) for p in psf], len(ret))
    return ret


def convolve_simple(dat, kernel):
    """simple convolution of two arrays"""
    npts = len(dat) #max(len(dat), len(kernel))
    pad = np.ones(npts)
    #tmp = np.concatenate((pad*dat[0], dat, pad*dat[-1]))
    tmp = np.pad(dat, (npts, npts), 'edge')
    out = np.convolve(tmp, kernel/np.sum(kernel), mode='valid')
    noff = int((len(out) - npts) * 0.5)
    #ret = (out[noff:])[:npts]
    ret = out[noff:noff+npts]
    #print(len(dat), len(kernel), len(ret))
    return ret


def create_xmin_xmax(x):
    mean = 0.5*(x[1:]+x[:-1])
    xmin = np.append(x[0], mean)
    xmax = np.append(mean, x[-1])
    return xmin, xmax


def detect_local_minima(arr):
    #https://stackoverflow.com/questions/3986345/how-to-find-the-local-minima-of-a-smooth-multidimensional-array-in-numpy-efficie
    # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    # apply the local minimum filter; all locations of minimum value
    # in their neighborhood are set to 1
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
    local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
    # local_min is a mask that contains the peaks we are
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    #
    # we create the mask of the background
    background = (arr==0)
    #
    # a little technicality: we must erode the background in order to
    # successfully subtract it from local_min, otherwise a line will
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    #
    # we obtain the final mask, containing only peaks,
    # by removing the background from the local_min mask
    detected_minima = local_min #- eroded_background
    #return np.where(detected_minima)
    return detected_minima

def expr_check(node):
    if isinstance(node, list):
        iter = node
    elif isinstance(node, ast.List):
        iter = node.elts
    else:
        return expr_eval(node)

    ret = []
    for i in iter:
        if isinstance(node, ast.Num):
            ret.append(float(i.n))
        elif isinstance(node, ast.Constant):
            ret.append(float(i.n))
        else:
            ret.append(expr_eval(i))
    try:
        ret = np.ravel(np.asarray(ret, dtype='float64'))
    except:
        ret = tuple([np.asarray(r, dtype='float64') for r in ret])
    return ret

def expr_eval(node):
    if isinstance(node, ast.Num): # <number>
        return node.n

    elif isinstance(node, ast.NameConstant): # <number>
        return node.value

    elif isinstance(node, ast.BinOp): # <left> <operator> <right>
        return py_ops[type(node.op)](expr_check(node.left),
                                     expr_check(node.right))

    elif isinstance(node, ast.Call):
        #print(node.args)
        #print(type(node.args))
        try:
            ret = getattr(np, expr_eval(node.func))(expr_check(node.args))
        except:
            ret = getattr(np, expr_eval(node.func))(*expr_check(node.args))
        return ret

    elif isinstance(node, ast.Compare):
        left = expr_check(node.left)
        for i, (o,c) in enumerate(zip(node.ops, node.comparators)):
            if i == 0:
                cond = py_ops[type(o)](left, expr_eval(c))
            else:
                cond = np.logical_and(cond, py_ops[type(o)](left, expr_eval(c)))
            left = expr_eval(c)
        return cond

    elif isinstance(node, ast.Name):
        return node.id

    elif isinstance(node, ast.UnaryOp): # <operator> <operand> e.g., -1
        return py_ops[type(node.op)](expr_eval(node.operand))
    else:
        #raise TypeError(node)
        return expr_check(node)

def lines_voigt_N_tot(x, z, N_tot, N_other, b, btur, series='Ly_a'):
#def lines_voigt_N_tot(x, z, logN_tot, logN_other, b, btur, series='Ly_a'):
    logN = np.log10(N_tot-N_other)
    if logN == -np.inf:
        logN = pars_std_d['logN']
    return lines_voigt(x, z, logN, b, btur, series)


def log2_range(start, end, step):
    start = np.log2(start)
    end = np.log2(end)
    start_r = np.floor if step < 0 else np.ceil
    end_r = np.ceil if step < 0 else np.floor
    if start % 1 > 0:
        start = np.round(start)
        logging.warning("I'm only using integer powers of 2. I changed 'start' "
                        "to %1.0f." % 2**start)
    if end % 1 > 0:
        end = np.round(end)
        logging.warning("I'm only using integer powers of 2. I changed 'end' "
                        "to %1.0f." % 2**end)
    log2 = np.arange(start, end+step, step)
    return np.power(2, log2)


def meta_parse(meta):
    s = ""
    for m in meta:
        if m not in forbidden_keywords and m[:5] not in forbidden_keywords:
            s += "%s: %s / %s \n" % (m, meta[m], meta.comments[m])
    return s[:-2]


def psf_gauss_wrong(x, #center, resol):
              resol, reg=None):
    """ @brief Gaussian PSF

    The function returns a gaussian array for each element of a selected region
    in the wavelength domain

    @param x Wavelength domain (in nm)
    @param c_min Starting pixel of the region
    @param c_max Ending pixel of the region
    @param center Center wavelength of the region
    @param resol Resolution
    @return Gaussian PSF over x
    """

    _, inds, _ = np.intersect1d(x, reg, return_indices=True)
    c = np.nanmedian(reg)
    sigma = c / resol * 4.246609001e-1
    psf = np.exp(-0.5*((x[inds]-c) / sigma)**2)
    #psf[np.where(psf < 1e-4)] = 0.0
    #psf = np.zeros(len(x))
    #psf[len(x)//2] = 1
    ret = [np.array(psf)]
    #plt.plot(x, psf)
    ret = psf
    return ret

def psf_gauss(x, resol, spec=None):
    c = x[len(x)//2]
    #resol = np.interp(c, spec.x, spec.t['resol'])
    #resol = 40000
    sigma = c / resol * 4.246609001e-1
    psf = np.exp(-0.5*((spec.x.to(xunit_def).value-c) / sigma)**2)
    psf = psf[np.where(psf > 1e-6)]
    #xout = spec.x.to(xunit_def).value[np.where(psf > 1e-6)]
    #psf[np.where(psf < 1e-4)] = 0.0
    #psf = np.zeros(len(x))
    #psf[len(x)//2] = 1
    #ret = [np.array(psf)]
    #plt.plot(xout*10, psf)
    """
    ret = psf
    return ret
    """
    if len(psf)==0:
        #print(x, spec.x.to(xunit_def).value)
        return psf_gauss(spec.x.to(xunit_def).value, resol, spec)
    else:
        ret = psf
        return ret
    #profile.disable()
    #ps = pstats.Stats(profile)
    #ps.print_stats()

def resol_check(spec, resol, prefix=prefix):
    check = resol is not None, 'resol' in spec.t.colnames
    resol = resol if check[0] else None
    print(msg_resol(check, prefix))
    return np.logical_or(*check), resol



def running_mean(x, h=1):
    """ From https://stackoverflow.com/questions/13728392/moving-average-or-running-mean """

    n = 2*h+1
    cs = np.nancumsum(np.insert(x, 0, 0))
    norm = np.nancumsum(~np.isnan(np.insert(x, 0, 0)))
    rm = (cs[n:] - cs[:-n]) / (norm[n:]-norm[:-n])
    return np.concatenate((h*[rm[0]], rm, h*[rm[-1]]))


def running_median(x, h=1):
    return median_filter(x, h)


def running_rms(x, xm, h=1):
    """ From https://stackoverflow.com/questions/13728392/moving-average-or-running-mean """
    n = 2*h+1
    rs = np.nancumsum(np.insert((x-xm)**2, 0, 0))
    norm = np.nancumsum(~np.isnan(np.insert(x, 0, 0)))
    #rm = (cs[n:] - cs[:-n]) / float(n)
    rms = np.sqrt((rs[n:] - rs[:-n]) / (norm[n:]-norm[:-n]))
    return np.concatenate((h*[rms[0]], rms, h*[rms[-1]]))


def to_x(z, trans):
    if trans == 'unknown':
        return z.to(au.nm)
    else:
        return (1+z)*xem_d[trans].to(au.nm)

def to_z(x, trans):
    if trans == 'unknown':
        return x.to(au.nm).value
    else:
        return (x.to(au.nm)/xem_d[trans].to(au.nm)).value-1


def trans_parse(series):
    trans = []
    for s in series.replace(';',',').replace(':',',').split(','):
        if '_' in s:
            trans.append(s.strip())
        else:
            for t in series_d[s.strip()]:
                trans.append(t)
    return trans


def elem_expand(elem, sess_sel):
    return '\n'.join([str(sess_sel)+','+r for r in elem.split('\n')])
    #return '\n'.join([str(sess_sel)+','+r+',C'+str(i%10) \
    #                 for i,r in enumerate(elem.split('\n'))])

# Adapted from http://ginstrom.com/scribbles/2008/09/07/getting-the-selected-cells-from-a-wxpython-grid/
def corners_to_cells(top_lefts, bottom_rights):
    """
    Take lists of top left and bottom right corners, and
    return a list of all the cells in that range
    """
    cells = []
    for top_left, bottom_right in zip(top_lefts, bottom_rights):

        rows_start = top_left[0]
        rows_end = bottom_right[0]

        cols_start = top_left[1]
        cols_end = bottom_right[1]

        rows = range(rows_start, rows_end+1)
        cols = range(cols_start, cols_end+1)

        cells.extend([(row, col)
            for row in rows
            for col in cols])

    return cells

def get_selected_cells(grid):
    """
    Return the selected cells in the grid as a list of
    (row, col) pairs.
    We need to take care of three possibilities:
    1. Multiple cells were click-selected (GetSelectedCells)
    2. Multiple cells were drag selected (GetSelectionBlock…)
    3. A single cell only is selected (CursorRow/Col)
    """

    top_left = grid.GetSelectionBlockTopLeft()

    if top_left:
        bottom_right = grid.GetSelectionBlockBottomRight()
        return corners_to_cells(top_left, bottom_right)

    selection = list(grid.GetSelectedCells())

    row = grid.GetGridCursorRow()
    col = grid.GetGridCursorCol()
    from wx.grid import GridCellCoords
    if not selection:
        return [GridCellCoords(row, col)]

    return [GridCellCoords(row, col)]+selection


def str_to_dict(str):
    return json.loads(str)


def range_str_to_list(range):
    return [[float(i) for i in r.split('-')] for r in range.split(',')]


def class_find(obj, cl, up=[]):
    if hasattr(obj, '__dict__'):
        for i in obj.__dict__:
            #print(obj, i)
            if isinstance(obj.__dict__[i], cl):
                print(up, i, 'caught!')
            else:
                class_find(obj.__dict__[i], cl, up+[i])
    elif isinstance(obj, dict):
        for i in obj:
            if isinstance(obj[i], cl):
                print(up, i, 'caught!')

def class_mute(obj, cl):
    if hasattr(obj, '__dict__'):
        for i in obj.__dict__:
            if isinstance(obj.__dict__[i], cl):
                obj.__dict__[i] = str(cl)
            elif not isinstance(obj.__dict__[i], Table):
            #else:
                class_mute(obj.__dict__[i], cl)
    elif isinstance(obj, dict):
        for i in obj:
            if isinstance(obj[i], cl):
                obj[i] = str(cl)

def class_unmute(obj, cl, targ):
    if hasattr(obj, '__dict__'):
        for i in obj.__dict__:
            if obj.__dict__[i]==str(cl):
                obj.__dict__[i] = targ
            else:
                class_unmute(obj.__dict__[i], cl, targ)
    elif isinstance(obj, dict):
        for i in obj:
            if obj[i]==str(cl):
                obj[i] = targ

def x_convert(x, zem=0, xunit=au.km/au.s):
    xem = (1+zem) * 121.567*au.nm
    equiv = [(au.nm, au.km/au.s,
              lambda x: np.log(x/xem.value)*ac.c.to(au.km/au.s),
              lambda x: np.exp(x/ac.c.to(au.km/au.s).value)*xem.value)]
    return x.to(xunit, equivalencies=equiv)

import functools
import warnings

"""
@decorator
class arg_fix:
    #Decorator ensuring backward compatibility when an argument name is
    #modified in a function definition.
    #from https://gist.github.com/rfezzani/002181c8667ec4c671421a4d938167eb

    def __init__(self, arg_mapping):
        #Args:
        #    arg_mapping (dict): mapping between the function's old argument
        #        names and the new ones.
        self.arg_mapping = arg_mapping
        self.warning_msg = ("'%s' is a deprecated argument name " +
                            "for the function '%s', use '%s' instead.")

    def __call__(self, f):
        @functools.wraps(f)
        def fixed_f(*args, **kwargs):
            for old_arg, new_arg in self.arg_mapping.items():
                if old_arg in kwargs:
                    #  warn that the function interface has changed:
                    warnings.warn(self.warning_msg %
                        (old_arg, f.__name__, new_arg), DeprecationWarning)
                    # Substitute new_arg to old_arg
                    kwargs[new_arg] = kwargs.pop(old_arg)

            # Call the function with the fixed arguments
            return f(*args, **kwargs)
        return fixed_f
"""
