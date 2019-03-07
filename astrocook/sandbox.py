
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
