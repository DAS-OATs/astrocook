from .cookbook_continuum_old import CookbookContinuumOld
from .functions import running_mean
from .message import *
from .utils import parse_range
from .vars import qso_composite, xem_d

import astropy.constants as ac
import astropy.table as at
import astropy.units as au
from scipy.interpolate import interp1d
import logging

class CookbookContinuum(CookbookContinuumOld):

    def __init__(self):
        super(CookbookContinuum, self).__init__()


    def clip_flux(self, zem, ran='all', smooth_len_lya=5000, smooth_len_out=400,
                  kappa=2, template=False, fudge='auto', knots_dist=2000,
                  mode='update'):


        """ @brief Clip flux
        @details Estimate the continuum by clipping absorbers.
        @url continuum_cb.html#clip-flux
        @param zem Emission redshift
        @param ran Wavelength range (nm)
        @param smooth_len_lya Smoothing length in the Lyman alpha forest (km/s)
        @param smooth_len_out Smoothing length outside the Lyman alpha forest (km/s)
        @param kappa Number of sigma to reject absorber
        @param template Use a composite spectrum to reduce the Lyman alpha emission peak
        @param fudge Fudge factor to scale the continuum
        @param knots_dist Distance between knots (km/s)
        @param mode Update or replace
        @return 0
        """
        
        try:
            zem = float(zem)
            xmin, xmax = parse_range(ran)
            smooth_len_lya = float(smooth_len_lya)
            smooth_len_out = float(smooth_len_out)
            kappa = float(kappa)
            template = str(template) == 'True'
            fudge = None if fudge=='auto' else float(fudge)
            knots_dist = float(knots_dist)
            mode = str(mode)
        except:
            logging.error(msg_param_fail)
            return 0

        if mode not in ['update', 'replace']:
            logging.warning("I cannot understand the mode. I will use `update`.")


        maxiter = 1000

        spec = self.sess.spec

        if 'cont_no_telluric' in spec._t.colnames:
            logging.info("An existing continuum included a telluric model. "
                         "I'm updating the continuum without changing the "
                         "telluric model.")
            spec._t['telluric'] = spec._t['cont']/spec._t['cont_no_telluric']
            spec._t['cont_telluric'] = spec._t['cont']
            spec._t['cont'] = spec._t['cont_no_telluric']


        dv = spec._dv()

        prox = 5000 * au.km/au.s
        lya_obs = xem_d['Ly_a']*(1+zem)
        lya_prox = (lya_obs*(1-prox/ac.c)).value

        # Extract range total
        xsel_tot = np.logical_and(spec._t['x']>xmin, spec._t['x']<xmax)
        spec_t_tot = spec._t[xsel_tot]

        exclude_criteria = spec_t_tot['y'] <= 3. * spec_t_tot['dy']
        #xsel_rej = np.logical_and(xsel_tot, ~exclude_criteria)
        xsel_rej = ~exclude_criteria
        spec_t_rej = spec_t_tot[xsel_rej]

        # Normalize to template
        if template:
            x_template = (qso_composite['x']/10.)* (1.+zem)
            y_template = qso_composite['y']
            template_interpolation = interp1d(x_template, y_template)
            y_interp = template_interpolation(spec_t_rej['x'])
            y_interp_cont = template_interpolation(spec_t_tot['x'])
            spec_t_rej['y'] /= y_interp
            spec_t_rej['dy'] /= y_interp

        dv = dv[xsel_tot][xsel_rej]

        # Extract ranges
        xpad_before = np.logical_and(spec._t['x']>xmin, spec._t['x']<lya_prox)
        xsel_before = np.logical_and(spec_t_rej['x']>xmin, spec_t_rej['x']<lya_prox)
        xpad_after = np.logical_and(spec._t['x']>lya_prox, spec._t['x']<xmax)
        xsel_after = np.logical_and(spec_t_rej['x']>lya_prox, spec_t_rej['x']<xmax)

        intervals = []
        if np.sum(xsel_before)>0:
            xsels_before = np.where(xpad_before==1)[0][0]
            spec_t_before = spec_t_rej[xsel_before]
            intervals.append((xsel_before, xsels_before, spec_t_before, smooth_len_lya))

        if np.sum(xsel_after)>0:
            xsels_after = np.where(xpad_after==1)[0][0]
            spec_t_after = spec_t_rej[xsel_after]
            intervals.append((xsel_after, xsels_after, spec_t_after, smooth_len_out))

        if 'mask_unabs' not in spec._t.colnames:
            logging.info("I'm adding column `mask_unabs`.")
        if 'cont' not in spec._t.colnames:
            logging.info("I'm adding column `cont`.")
        if mode == 'replace' or 'mask_unabs' not in spec._t.colnames:
            spec._t['mask_unabs'] = at.Column(np.zeros(len(spec._t)), dtype=int)
        if mode == 'replace' or 'cont' not in spec._t.colnames:
            spec._t['cont'] = at.Column(np.array(None, ndmin=1), dtype=float)

        x_rm = np.array([])
        y_rm = np.array([])
        sel_tot = np.array([])
        for xsel, xsels, spec_t, smooth_len_local in intervals:

            # Compute running median
            hwindow = int(smooth_len_local/np.min(dv[xsel]))//8

            import matplotlib.pyplot as plt
            plt.plot(spec_t['x'], spec_t['y'])
            #plt.show()
            y_rm_interval = running_mean(spec_t['y'], h=hwindow)
            # Clip flux
            sum_sel = len(spec_t)
            sel = np.zeros(len(spec_t), dtype=bool)
            for i in range(maxiter):
                sel[~sel] = spec_t['y'][~sel]-y_rm_interval<-kappa*spec_t['dy'][~sel]
                selw = np.where(sel==1)[0]#+xsels
                x_rm_interval = spec_t['x'][~sel]
                y_rm_interval = running_mean((spec_t['y'])[~sel], h=hwindow)

                if i == maxiter-1 and sum_sel != np.sum(sel):
                    logging.warning("Clipping not converged after {} iterations! "\
                                    .format(i+1))
                if sum_sel != np.sum(sel):
                    sum_sel = np.sum(sel)

                else:
                    logging.info("Clipping converged after %i iterations." % (i+1))
                    break
            x_rm = np.append(x_rm, x_rm_interval)
            y_rm = np.append(y_rm, y_rm_interval)
            sel_tot = np.append(sel_tot,sel)

        _, _, mask_unabs = np.intersect1d(x_rm, spec._t['x'], return_indices=True)
        spec._t['mask_unabs'][mask_unabs] = np.ones(len(mask_unabs))

        # Compute fudge if `auto` and apply it

        #"""
        if fudge is None:
            diff = running_mean(spec_t_rej['y'][np.logical_not(sel_tot)]-y_rm, hwindow)
            fudge = 1+diff/y_rm
        y_rm *= fudge
        #"""
        cont_temp = np.interp(spec._t['x'][xsel_tot], x_rm, y_rm)*spec_t_tot['y'].unit

        if template: cont_temp*=y_interp_cont

        # Update old continuum
        if mode=='update' and 'cont' in spec._t.colnames:
            cont_old = spec._t['cont']
        else:
            cont_old = np.array([np.nan]*len(spec._t))
        cont_old[xsel_tot] = cont_temp
        spec._t['cont_old'] = cont_old
        spec._gauss_convolve(std=smooth_len_out, input_col='cont_old',
                             output_col='cont')
        spec._t.remove_column('cont_old')

        # Propagate to lines
        lines = self.sess.lines
        if lines is not None:
            lines._t['cont'] = np.interp(lines._t['x'], spec._t['x'],
                                         spec._t['cont'])

        # Extract nodes
        self.nodes_extract(delta_x=knots_dist, mode='cont')

 
        if 'cont_no_telluric' in spec._t.colnames:
            spec._t['cont_no_telluric'] = spec._t['cont']
            spec._t['cont'] = spec._t['cont']*spec._t['telluric']
            spec._t.remove_columns(['cont_telluric', 'telluric'])

        return 0


    def fit_pl(self):
        """@brief Fit power law 🚧
        @details 🚧
        @url continuum_cb.html#fit-power-law
        """

        return 0


    def update_deabs(self):
        """@brief Update after de-absorbing 🚧
        @details 🚧
        @url continuum_cb.html#update-after-de-absorbing
        """

        return 0


    def correct_lya(self, zem, input_col='y', mode='basic', logN_thres=100,
                 percentile=100):
        """ @brief Correct for Ly-a opacity
        @details Correct the spectrum flux for the effective Lyman-alpha opacity.
        @url continuum_cb.html#correct-for-ly-a-opacity
        @param zem Emisson redshift
        @param input_col Column to correct
        @param mode Correction mode ('basic' or 'inoue')
        @param logN_col Threshold for logarithmic column density
        @param percentile Percentile to compute the threshold from the column
        density distribution (only if logN_col is None)
        @return 0
        """

        try:
            zem = float(zem)
            logN_thres = None if logN_thres in [None, 'None'] else float(logN_thres)
        except:
            logging.error(msg_param_fail)
            return 0

        spec = self.sess.spec
        systs = self.sess.systs

        if logN_thres is None:
            try:
                #logN_thres = np.median(self.sess.systs._t['logN'])
                logN_thres = np.percentile(self.sess.systs._t['logN'], percentile)
                logging.info("I estimated the threshold column density from "
                             "the system table: logN_thres = %2.2f." \
                             % logN_thres)
            except:
                logN_thres = 100
                logging.warning("No systems: I couldn't estimate the threshold "
                                "column density. I am using logN_thres = "\
                                "%2.2f." % logN_thres)
        if mode == 'inoue':
            inoue_all = getattr(spec, '_lya_corr_inoue')(zem, input_col,
                                                         apply=False)
            basic_all = getattr(spec, '_lya_corr_basic')(zem, 100, input_col,
                                                         apply=False)
        else:
            inoue_all = np.zeros(len(spec.x))
            basic_all = np.zeros(len(spec.x))
        basic = getattr(spec, '_lya_corr_basic')(zem, logN_thres, input_col,
                                                 apply=False)


        self._lya_corr = 1+(inoue_all-1)*(basic-1)/(basic_all-1)
        """
        plt.plot(spec.x, inoue_all)
        plt.plot(spec.x, basic_all)
        plt.plot(spec.x, basic)
        plt.plot(spec.x, self._lya_corr, linewidth=3)
        """
        #plt.show()
        taucorr = dc(spec._t[input_col])
        taucorr *= self._lya_corr
        spec._t[input_col+'_taucorr'] = taucorr

        w = spec.x.value > (1+zem)*xem_d['Ly_lim'].value
        logging.info("Mean correction bluewards from Lyman limit: %3.2f." \
                     % np.mean(self._lya_corr[w]))
        return 0
