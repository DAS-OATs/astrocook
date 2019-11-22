
    def add_fit(self, series='CIV', z=1.6971, logN=13, b=10, resol=70000,
                chi2r_thres=None, maxfev=100, append=True):
        """ @brief Add and fit a Voigt model for a system.
        @param series Series of transitions
        @param z Guess redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param maxfev Maximum number of function evaluation
        @param append Append system to existing system list
        @return 0
        """

        z = float(z)
        logN = float(logN)
        b = float(b)
        resol = float(resol)
        if chi2r_thres == None or chi2r_thres == 'None':
            chi2r_thres = np.inf
        else:
            chi2r_thres = float(chi2r_thres)
        maxfev = int(maxfev)

        systs = SystList(sess=self)
        #mods = ModelList()
        if append and self.systs != None:
            self.systs._append(systs)
            #self.systs._mods._append(systs._mods)
        else:
            self.systs = systs

        #"""
        self.systs._add_fit(series, z, logN, b, resol, chi2r_thres,
                            fit_kws={'maxfev': maxfev})
        """
        mod = SystModel2(self.systs)

        self.systs._t.add_row(['voigt_func', series, z, z, logN, b, None, self.systs._id])
        mod.new_voigt(series, z, logN, b, resol)
        mod.fit(fit_kws={'maxfev': maxfev})
        self.systs._fit_save(mod)
        print(self.systs._t)
        #"""

        if chi2r_thres != np.inf:
            z_rem = self.systs._clean(chi2r_thres)
        else:
            z_rem = []

        self._update_spec()

        return 0

    def add_fit_from_lines(self, series='CIV', z_start=1.71, z_end=1.18,
                           dz=1e-4, logN=14, b=10, resol=70000,
                           chi2r_thres=None, maxfev=100, append=True):
        """ @brief Add and fit Voigt models to a line list, given a redshift
        range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param dz Threshold for redshift coincidence
        @param N Guess column density
        @param b Guess doppler broadening
        @param resol Resolution
        @param chi2r_thres Reduced chi2 threshold to accept the fitted model
        @param maxfev Maximum number of function evaluation
        @param append Append systems to existing system list
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        dz = float(dz)
        logN = float(logN)
        b = float(b)
        resol = float(resol)
        if chi2r_thres == None or chi2r_thres == 'None':
            chi2r_thres = np.inf
        else:
            chi2r_thres = float(chi2r_thres)
        maxfev = int(maxfev)

        z_range = self.lines._syst_cand(series, z_start, z_end, dz)

        self.systs = SystList(sess=self)
        #mods = ModelList()
        self.systs._add_fit(series, z_range, logN, b, resol, chi2r_thres,
                            fit_kws={'maxfev': maxfev}, verb=True)
        self._update_spec()

        return 0


# syst_list 2019-03-11

    def _domain(self, ys, thres=1e-3): #series, z, zmin, zmax):
        """ @brief Define domain for fitting. """

        spec = self._spec
        #m = mod.eval(x=self._xs, params=pars)
        #c = np.where(m<1-thres)
        c = np.where(ys<1-thres)[0]
        print(c)
        print(np.ediff1d(c))
        print(np.where(np.ediff1d(c)>1))

        #"""
        xt = np.array(self._xs[c])
        xc = np.split(xt, np.where(np.ediff1d(c)>1)[0])

        if 'deabs' in spec._t.colnames:# and 1 == 0:
            yc = np.array(spec._t['deabs'][c]/spec._t['cont'][c])
        else:
            yc = np.array(spec.y[c]/spec._t['cont'][c])
        wc = np.array(spec._t['cont'][c]/spec.dy[c])
        """
        xc = self._xs
        yc = np.full(len(xc), np.nan)
        wc = np.full(len(xc), np.nan)
        if 'deabs' in spec._t.colnames:# and 1 == 0:
            yc[c] = np.array(spec._t['deabs'][c]/spec._t['cont'][c])
        else:
            yc[c] = np.array(spec.y[c]/spec._t['cont'][c])
        wc[c] = np.array(spec._t['cont'][c]/spec.dy[c])
        #"""

        return xc, yc, wc

    def _fit(self, mod, pars, xc, yc, wc):


        #pars.pretty_print()
        fit = mod.fit(yc, pars, x=xc[0], weights=wc)#, nan_policy='omit')
        mod = fit
        plt.plot(self._xs, mod.eval(x=[self._xs], params=fit.params))
        #plt.plot(self._xs, mod._lines[0].eval(x=self._xs, params=fit.params))
        #plt.plot(self._xs, mod._psf.eval(x=self._xs, params=fit.params))
        plt.plot(xc, yc)
        plt.show()
        pars = fit.params
        #pars.pretty_print()
        self._ys = mod.eval(x=[self._xs], params=pars)

        return mod, pars

    def _group(self, mod, pars, ys, thres=1.e-3):
        """ @brief Define group of systems. A group is the set of systems with
        at least one line within the footprint of the system described by the
        model in input.
        """

        group = [len(self._t)-1]
        c = 0
        for i, s in enumerate(self._t[:-1]):  # All systems except the last
            yn = s['ys']
            if np.amin(np.maximum(ys, yn)) < 1-thres:
                group.append(i)
                m = ModelLines(lines_voigt, c+1, s['series'], s['z'], s['zmin'],
                               s['zmax'], **s['pars'])
                mod *= m
                pars.update(m._pars)

                c += 1

        group = np.array(group)
        ys = mod.eval(x=[self._xs], params=pars)
        return mod, pars, ys, group

    def _single_std(self, series='Ly_a', z=2.0, N=1e13, b=10, btur=0,
                    resol=35000, psf=False, eval=True, add=True):

        # Create model
        func = 'voigt'
        dz = 0.0
        z_tol = 1e-4
        zmin = z-z_tol
        zmax = z+z_tol
        pars = {'N': N, 'b': b, 'b_tur': btur, 'resol': resol}
        dpars = {'N': None, 'b': None, 'b_tur': None, 'resol': None}
        chi2r = None

        """
        lines = ModelLines(lines_voigt, 0, series, z, zmin, zmax, **pars)
        if psf:
            psf = ModelPSF(psf_gauss, 0, series, z, zmin, zmax, **pars)
            mod = LMComposite(lines, psf, convolve)
            mod._pars = lines._pars
            mod._pars.update(psf._pars)
        else:
            mod = lines
        pars = mod._pars
        pars.pretty_print()
        """
        mod = SystModelStd(0, series, z, zmin, zmax, **pars)
        #"""
        if eval:
            ys = mod.eval(x=[self._xs], params=mod._pars)
        else:
            ys = None

        plt.plot(self._xs, mod.eval(x=[self._xs], params=mod._pars))
        plt.plot(self._xs, mod._lines[0].eval(x=[self._xs], params=mod._pars))
        plt.plot(self._xs, mod._psf.eval(x=[self._xs], params=mod._pars)[0])
        plt.show()

        if add:
            self._t.add_row([series, func, z, dz, zmin, zmax, mod, pars, dpars,
                             chi2r, ys])

        return mod, mod._pars, ys

    def _update_spec(self, fit=None):#, fit):
        """ @brief Update spectrum after fitting """

        spec = self._spec
        y = spec.y
        if 'model' not in spec._t.colnames:
            print(prefix, "I'm adding column 'model'.")
            spec._t['model'] = np.empty(len(spec.x), dtype=float)*y.unit
        if 'deabs' not in spec._t.colnames:
            spec._t['deabs'] = y

        s = self._s
        cont = spec._t['cont']
        model = spec._t['model']
        deabs = spec._t['deabs']

        model[s] = cont[s]
        for i, r in enumerate(self._t):
            #m = Model(series=r['series'], z=r['z'], zmin=r['zmin'],
            #          zmax=r['zmax'], pars=r['pars'])
            m = ModelLines(lines_voigt, i, r['series'], r['z'], r['zmin'],
                           r['zmax'], **r['pars'])
            #model[s] = m._mod.eval(x=self._xs, params=m._mod._pars) * model[s]
            model[s] = m.eval(x=[self._xs], params=m._pars) * model[s]
        deabs[s] = cont[s] + y[s] - model[s]


    def _update_systs(self, mod, pars, ys, group=None):
        """ @brief Update system list after fitting """

        #pars = fit.params
        if group is None:
            group = [len(self._t)-1]
        for i, w in enumerate(group):
            self._t[w]['mod'] = mod
            self._t[w]['z'] = pars['lines_voigt_'+str(i)+'_z'].value \
                         #*au.dimensionless_unscaled
            self._t[w]['dz'] = pars['lines_voigt_'+str(i)+'_z'].stderr\
                          #*au.dimensionless_unscaled
            #"""
            self._t[w]['pars'] = {'N': pars['lines_voigt_'+str(i)+'_N'].value,
                         'b': pars['lines_voigt_'+str(i)+'_b'].value,
                         'btur': pars['lines_voigt_'+str(i)+'_btur'].value}
                         #'resol': pars['psf_gauss_'+str(i)+'_resol'].value}
                         #'ampl': pars['adj_gauss_'+str(i)+'_ampl'].value}
            self._t[w]['dpars'] = {'N': pars['lines_voigt_'+str(i)+'_N'].stderr,
                          'b': pars['lines_voigt_'+str(i)+'_b'].stderr,
                          'btur': pars['lines_voigt_'+str(i)+'_btur'].stderr}
                          #'resol': pars['psf_gauss_'+str(i)+'_resol'].stderr}
                          #'ampl': pars['adj_gauss_'+str(i)+'_ampl'].stderr}
            #"""
            self._t[w]['chi2r'] = mod.redchi

    def fit_from_lines(self, series='Ly_a', z_start=2.5, z_end=2.0, N=1e14,
                       b=10, thres=1e-3):
        """ @brief Fit Voigt models to a line list, given a redshift range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param N Guess column density
        @param b Guess doppler broadening
        @param thres Threshold for grouping
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        N = float(N)
        b = float(b)
        thres = float(thres)

        z_lines = [[(x.to(au.nm)/xem_d[t].to(au.nm))-1. \
                    for t in series_d[series]] for x in self._lines.x]
        z_lines = np.ravel(z_lines)
        if z_end < z_start:
            z_range = z_lines[np.logical_and(z_lines<z_start, z_lines>z_end)]\
                          [::-1]
        else:
            z_range = z_lines[np.logical_and(z_lines>z_start, z_lines<z_end)]

        #for z in np.arange(z_start, z_end, z_step):
        for z in z_range:
            self.fit(series=series, z=z, N=N, b=b, group_thres=thres,
                     update=False)
            print(prefix, "I've fitted a system at redshift %2.4f..." % z,
                  end='\r')
        print(prefix, "I've fitted all systems between redshift %2.4f and "\
              "%2.4f." % (z_start, z_end))

        self._update_spec()

        return 0

    def fit_from_deabs(self, series='Ly_a', z_start=2.5, z_end=2.0, N=1e12,
                       b=10, thres=1e-3):
        """ @brief Fit Voigt models to spectrum residuals, given a redshift
        range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param N Guess column density
        @param b Guess doppler broadening
        @param thres Threshold for grouping
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        N = float(N)
        b = float(b)
        thres = float(thres)

        spec_deabs = dc(self._spec)
        spec_deabs.gauss_convolve(input_col='deabs')
        lines_deabs = spec_deabs.peaks_find()

        z_lines = [[(x.to(au.nm)/xem_d[t].to(au.nm))-1. \
                    for t in series_d[series]] for x in lines_deabs.x]
        z_lines = np.ravel(z_lines)
        if z_end < z_start:
            z_range = z_lines[np.logical_and(z_lines<z_start, z_lines>z_end)]\
                          [::-1]
        else:
            z_range = z_lines[np.logical_and(z_lines>z_start, z_lines<z_end)]

        #for z in np.arange(z_start, z_end, z_step):
        for z in z_range:
            self.fit(series=series, z=z, N=N, b=b, group_thres=thres,
                     update=False)
            print(prefix, "I've fitted a system at redshift %2.4f..." % z,
                  end='\r')
        print(prefix, "I've fitted all systems between redshift %2.4f and "\
              "%2.4f." % (z_start, z_end))

        self._update_spec()

        return 0

    def fit_range(self, series='Ly_a', z_start=2.5, z_end=2.0, z_step=1e-3,
                  N=1e13, b=10, thres=1e-3):
        """ @brief Fit Voigt models at constant steps in redshift range.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param z_step Redshift step
        @param N Guess column density
        @param b Guess doppler broadening
        @param thres Threshold for grouping
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        z_step = float(z_step)
        N = float(N)
        b = float(b)
        thres = float(thres)

        for z in np.arange(z_start, z_end, z_step):
            self.fit(series=series, z=z, N=N, b=b, group_thres=thres,
                     update=False)
            print(prefix, "I've fitted a system at redshift %2.4f..." % z,
                  end='\r')
        print(prefix, "I've fitted all systems between redshift %2.4f and "\
              "%2.4f." % (z_start, z_end))

        self._update_spec()

        return 0

    def fit_slide(self, series='CIV', z_start=1.13, z_end=1.71, z_step=5e-4,
                  logN=14, b=10.0):
        """ @brief Slide a set of Voigt models across a spectrum and fit them
        where they suit the spectrum.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param z_step Redshift step
        @param logN Column density (logarithmic)
        @param b Doppler parameter
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        z_step = float(z_step)
        N = 10**float(logN)
        b = float(b)

        z_range = np.arange(z_start, z_end, z_step)

        xs = self._xs
        #z_arr = []

        count = 0

        # Candidate model
        mod, pars, _ = self._single_std(series, 0., N, b, eval=False, add=False)

        # Null model

        chi2_arr = []
        for z in z_range:
            print(prefix, "I'm scanning the spectrum: now at redshift %2.4f..." \
                  % z, end='\r')
            self._spec._shift_rf(z)
            self._xs = np.array(self._spec._safe(self._spec.x).to(au.nm))
            ys = mod.eval(x=self._xs, params=pars)
            xc, yc, wc = self._domain(ys)
            ym = mod.eval(x=xc, params=pars)
            chi2 = np.sum(((ym-yc)/wc)**2)
            chi2_arr.append(chi2)


        chi2_sm = running_mean(chi2_arr, 1)
        #plt.plot(z_range, chi2_arr/chi2_sm)
        #plt.show()
        chi2_found = sigma_clip(chi2_arr/chi2_sm, sigma=5)
        z_found = z_range[chi2_found.mask]

        for z_cen in z_found:
            print(prefix, "I'm checking a candidate at redshift %2.4f...      " \
                  % z_cen)#, end='\r')
            for z in np.arange(z_cen-1e-4, z_cen+1e-4, 1e-5):
                self._spec._shift_rf(z)

                self._xs = np.array(self._spec._safe(self._spec.x).to(au.nm))
                ys = mod.eval(x=self._xs, params=pars)
                xc, yc, wc = self._domain(ys)
                ym = mod.eval(x=xc, params=pars)
                chi2 = np.sum(((ym-yc)*wc)**2)/len(ym)

                ym_0 = np.ones(len(xc))
                chi2_0 = np.sum(((ym_0-yc)*wc)**2)/len(ym)

                ym_1 = np.append(ym[:len(ym)//2], np.ones(len(xc)-len(ym)//2))
                chi2_1 = np.sum(((ym_1-yc)*wc)**2)/len(ym)

                ym_2 = np.append(np.ones(len(ym)//2), ym[len(ym)//2:])
                chi2_2 = np.sum(((ym_2-yc)*wc)**2)/len(ym)

                #print("%2.3e %2.3e %2.3e %2.3e" % (chi2, chi2_0, chi2_1, chi2_2))
                if chi2 < chi2_0 and chi2 < chi2_1 and chi2 < chi2_2 :
                    self._spec._shift_rf(0)
                    self._xs = np.array(self._spec._safe(self._spec.x).to(au.nm))
                    mod_f, pars_f, ys_f = self._single_std(series, z, N, b)
                    mod_f, pars_f, ys_f, group_f = self._group(mod_f, pars_f,
                                                               ys_f)
                    xc_f, yc_f, wc_f = self._domain(ys_f)
                    mod_f, pars_f = self._fit(mod_f, pars_f, xc_f, yc_f, wc_f)
                    self._update_systs(mod_f, pars_f, group_f)
                    print(prefix, "I've fitted a system at redshift %2.4f."\
                          "        " % z)
                    #plt.plot(xc, yc)

                    try:
                        self._update_spec(mod_f)
                    except:
                        pass
                    plt.plot(xc, ym)
                    plt.plot(xc, ym_0)
                    plt.plot(xc, ym_1)
                    plt.plot(xc, ym_2)
                    plt.show()

        self._spec._shift_rf(0)

        return 0


    def fit(self, series='Ly_a', z=2.0, N=1e13, b=10, btur=0, resol=35000,
            ampl=0.0, group_thres=1e-3, domain_thres=1e-3, update=True):
        """ @brief Create and fit a Voigt model for a system.
        @param series Series of transitions
        @param z Redshift
        @param N Guess column density
        @param b Guess Doppler broadening
        @param btur Guess turbulence broadening
        @param resol Resolution
        @param ampl Amplitude of the continuum adjustment
        @param group_thres Threshold for grouping
        @param domain_thres Threshold for fitting domain
        @param update Flag to update the spectrum
        @return 0
        """

        z = float(z)
        N = float(N)
        b = float(b)
        btur = float(btur)
        resol = float(resol)
        ampl = float(ampl)
        group_thres = float(group_thres)
        domain_thres = float(domain_thres)

        mod, pars, ys = self._single_std(series, z, N, b, btur, resol, psf=True)
        mod, pars, ys, group = self._group(mod, pars, ys, group_thres)
        #mod, pars, ys, group = self._psf(mod, pars, ys)
        xc, yc, wc = self._domain(ys, domain_thres)
        mod, pars = self._fit(mod, pars, xc, yc, wc)
        self._update_systs(mod, pars, group)
        if update:
            self._update_spec(mod)

        return 0


###### OLD

    # Method of SystList
    def fit_slide(self, series='Ly_a', z_start=2.5, z_end=2.0, z_step=-1e-3,
                  logN_start=16, logN_end=13, logN_step=-0.5,
                  b_start=5., b_end=105., b_step=10.):
        """ @brief Slide a set of Voigt models across a spectrum and fit them
        where they suit the spectrum.
        @param series Series of transitions
        @param z_start Start redshift
        @param z_end End redshift
        @param z_step Redshift step
        @param logN_start Start column density (logarithmic)
        @param logN_end End column density (logarithmic)
        @param logN_step Column density step (logarithmic)
        @param b_start Start Doppler parameter
        @param b_end End Doppler parameter
        @param b_step Doppler parameter step
        @return 0
        """

        z_start = float(z_start)
        z_end = float(z_end)
        z_step = float(z_step)
        logN_start = float(logN_start)
        logN_end = float(logN_end)
        logN_step = float(logN_step)
        b_start = float(b_start)
        b_end = float(b_end)
        b_step = float(b_step)

        z_range = np.arange(z_start, z_end, z_step)
        N_range = 10**np.arange(logN_start, logN_end, logN_step)
        b_range = np.arange(b_start, b_end, b_step)

        xs = self._xs
        #z_arr = []

        count = 0
        for z in z_range:
            chi2_arr = np.empty((len(N_range),len(b_range)))
            chi2_n_arr = np.empty((len(N_range),len(b_range)))
            print(prefix, "I'm scanning the spectrum: now at redshift %2.4f" \
                  % z, end='\r')

            for i, N in enumerate(N_range):

                for j, b in enumerate(b_range):

                    mod, pars, ys = self._single_std(series, z, N, b, add=False)
                    xc, yc, wc = self._domain(ys)
                    ym = mod.eval(x=xc, params=pars)
                    chi2 = np.sum(((ym-yc)/wc)**2)
                    chi2_arr[i, j] = chi2

                    # Null model
                    mod_n, pars_n, ys_n = self._single_std(series, z, 0., add=False)
                    xc_n, yc_n, wc_n = self._domain(ys_n)
                    ym_n = mod_n.eval(x=xc_n, params=pars_n)
                    chi2_n = np.sum(((ym_n-yc_n)/wc_n)**2)
                    chi2_n_arr[i, j] = chi2_n

            #print(chi2_arr)
            #print(chi2_n_arr)
            min_chi2 = np.amin(chi2_arr)
            argmin_chi2 = np.unravel_index(np.argmin(chi2_arr, axis=None),
                                           chi2_arr.shape)
            #print(min_chi2)
            #print(chi2_n_arr[argmin_chi2])

            if (min_chi2<chi2_n_arr[argmin_chi2]):
                count += 1
                N = N_range[argmin_chi2[0]]
                b = b_range[argmin_chi2[1]]
                print(z, chi2, chi2_n, N, b)
                #plt.plot(xc, yc)
                #plt.plot(xc, ym)
                #plt.show()
                mod_f, pars_f, ys_f = self._single_std(series, z, N, b)
                #pars_f.pretty_print()
                mod_f, pars_f, ys_f, group_f = self._group(mod_f, pars_f, ys_f)
                #pars_f.pretty_print()
                xc_f, yc_f, wc_f = self._domain(ys_f)
                mod_f, pars_f = self._fit(mod_f, pars_f, xc_f, yc_f, wc_f)
                #pars_f.pretty_print()
                #z_arr.append(z)
                self._update_systs(mod_f, pars_f, group_f)
        print(prefix, "I've scanned the spectrum between redshift %2.4f and "\
              "%2.4f and fitted %i systems." % (z_start, z_end, count))
            #print(z_arr)

        #plt.plot(z_range, chi2_arr)
        #plt.plot(z_range, chi2_n_arr)
        #plt.show()
        try:
            self._update_spec(mod_f)
        except:
            pass


        """
        for N in N_range:
            chi2_arr = []
            chi2_n_arr = []
            count = 0
            for z in z_range:
                print(prefix, "I'm scanning the spectrum: now at redshift %2.4f" \
                    % z, end='\r')
                comp, pars = self._model_voigt(series, z, N, b=20, add=False)
                #comp, pars, group = self._group(comp, pars)
                xc, yc, wc = self._domain(comp, pars)
                ym = comp.eval(x=xc, params=pars)
                chi2 = np.sum(((ym-yc)/wc)**2)
                chi2_arr.append(chi2)

                # Null model
                comp_n, pars_n = self._model_voigt(series, z, 0., b=20, add=False)
                #comp_n, pars_n, group_n = self._group(comp_n, pars_n)
                xc_n, yc_n, wc_n = self._domain(comp_n, pars_n)
                ym_n = comp_n.eval(x=xc_n, params=pars_n)
                chi2_n = np.sum(((ym_n-yc_n)/wc_n)**2)
                chi2_n_arr.append(chi2_n)

                if (chi2<0.5*chi2_n):
                    count += 1
                    print(z, chi2, chi2_n)
                    plt.plot(xc, yc)
                    plt.plot(xc, ym)
                    plt.show()
                    comp_f, pars_f = self._model_voigt(series, z, N)
                    #pars_f.pretty_print()
                    comp_f, pars_f, group_f = self._group(comp_f, pars_f)
                    #pars_f.pretty_print()
                    xc_f, yc_f, wc_f = self._domain(comp_f, pars_f)
                    comp_f, pars_f = self._fit_voigt(comp_f, pars_f, xc_f, yc_f, wc_f)
                    #pars_f.pretty_print()
                    #z_arr.append(z)
                    self._update_systs(comp_f, group_f)
            print(prefix, "I've scanned the spectrum between redshift %2.4f and "\
                  "%2.4f for column density %3.2e and fitted %i systems."
                  % (z_start, z_end, N, count))
            #print(z_arr)

            plt.plot(z_range, chi2_arr)
            plt.plot(z_range, chi2_n_arr)
            plt.show()
            try:
                self._update_spec(comp_f)
            except:
                pass
        """

        return 0


class ModelAdj(LMModel):
    """ Class for continuum models

    A ModelAdj is a model for the local continuum."""

    def __init__(self, call, count,
                 **kwargs):
        func = call._adj_func
        func_name = func.__name__
        if func_name != 'adj_gauss':
            print(prefix, "Only 'adj_gauss' function supported for lines.")
            return None
        self._prefix = func_name+'_'+str(count)+'_'
        super(ModelAdj, self).__init__(func, prefix=self._prefix)
        getattr(self, '_pars_'+func_name)(call)

    def _pars_adj_gauss(self, call):
        d = call._adj_d
        self._pars = self.make_params()
        self._pars.add_many(
            (self._prefix+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            (self._prefix+'ampl', d['ampl'], d['ampl_vary'], d['ampl_min'],
             d['ampl_max'], d['ampl_expr']),
            (self._prefix+'sigma', d['sigma'], d['sigma_vary'], d['sigma_min'],
             d['sigma_max'], d['sigma_expr']))

class Model(object):
    """ Class for models

    A Model is a combination of Lmfit Models for instrument PSF,
    continuum adjustment, and system profile."""

    def __init__(self,
                 series='Ly_a',
                 lines_func=lines_voigt,
                 psf_func=psf_gauss,
                 adj_func=adj_gauss,
                 z=None,
                 zmin=None,
                 zmax=None,
                 pars=None,
                 count=0):
        self._series = series
        self._lines_func = lines_func
        self._psf_func = psf_func
        self._adj_func = adj_func
        self._z = z
        self._zmin = zmin
        self._zmax = zmax
        self._pars = pars
        self._count = count
        self._set_defaults()
        self._create_comp()

    def _create_comp(self):
        #adj = ModelAdj(self, self._count, **self._adj_d)
        lines = ModelLines(self, self._count, self._series, **self._lines_d)#._lines_func, series=self._series)
        psf = ModelPSF(self, self._count, **self._psf_d)
        #adj_lines = adj * lines
        adj_lines = lines
        self._comp = LMComposite(adj_lines, psf, convolve)
        #self._comp._pars = adj._pars
        self._comp._pars = lines._pars
        #self._comp._pars.update(adj._pars)
        self._comp._pars.update(psf._pars)

    def _set_defaults(self):

        # Line defaults
        """
        self._adj_d = adj_gauss_d
        self._adj_d['z'] = self._z
        self._adj_d['z_min'] = self._zmin
        self._adj_d['z_max'] = self._zmax
        """
        self._lines_d = lines_voigt_d
        self._lines_d['z'] = self._z
        self._lines_d['z_min'] = self._zmin
        self._lines_d['z_max'] = self._zmax
        self._psf_d = psf_gauss_d
        self._psf_d['z'] = self._z
        self._psf_d['zmin'] = self._zmin
        self._psf_d['zmax'] = self._zmax
        if self._pars is not None:
            for p in self._pars:
                """
                if p in self._adj_d:
                    self._adj_d[p] = self._pars[p]
                """
                if p in self._lines_d:
                    self._lines_d[p] = self._pars[p]
                if p in self._psf_d:
                    self._psf_d[p] = self._pars[p]





class ModelPSF(LMModel):
    """ Class for psf models

    A ModelPSF is a model for the instrumental PSF."""

    def __init__(self, call, count,
                 **kwargs):
        func = call._psf_func
        func_name = func.__name__
        if func_name != 'psf_gauss':
            print(prefix, "Only 'psf_gauss' function supported for PSF.")
            return None
        self._prefix = func_name+'_'+str(count)+'_'
        super(ModelPSF, self).__init__(func, prefix=self._prefix)#, **kwargs)
        getattr(self, '_pars_'+func_name)(call)

    def _pars_psf_gauss(self, call):
        d = call._psf_d
        self._pars = self.make_params()
        self._pars.add_many(
            (self._prefix+'z', d['z'], d['z_vary'], d['z_min'], d['z_max'],
             d['z_expr']),
            (self._prefix+'resol', d['resol'], d['resol_vary'], d['resol_min'],
             d['resol_max'], d['resol_expr']))
