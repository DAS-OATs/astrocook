from astropy.io import ascii
import datetime as dt
import logging

class Workflow(object):


    def __init__(self, gui, cb):
        self._gui = gui
        self._cb = cb

    def civ_full(self):
        cb = self._cb
        path = '/Users/guido/GitHub/test_data/CIV/'
        targ_list = ascii.read(path+'targets.csv')

        logging.info("I'm launching civ_full...")
        for l in targ_list:
            start = dt.datetime.now()
            logging.info("I'm starting with target %s..." % l['name'])

            zem = l['z']
            xmin = l['lambdamin']
            xmax = l['lambdamax']
            self._gui._panel_sess._on_open(path+'reduced/'+l['name']+'.fits')

            sess_start = self._gui._sess_sel
            if sess_start.spec.meta['object'] == 'J2123-0050':
                sess = sess_start.cb.region_extract(xmin=xmin, xmax=xmax)
            else:
                sess = sess_start

            """
            self._cb.gauss_convolve(std=10)
            self._cb.peaks_find(kappa=3.0)
            sess.lines._t.remove_rows(sess.lines.y == 0)
            if np.mean(sess.spec._t['y'])<1 and np.std(sess.spec._t['y'])<1:
                sess.spec._t['cont'] = [1] * len(sess.spec._t)*sess.spec.y.unit
            if 'cont' not in sess.spec._t.colnames:
                self._cb.nodes_extract(delta_x=1000)
                self._cb.nodes_interp()

            sess_reg = sess.cb.region_extract(xmin=xmin, xmax=xmax)
            self._gui._panel_sess._on_add(sess_reg, open=False)
            sess_center = dc(sess_reg)
            sess_center.add_syst_from_lines(z_end=20, maxfev=10)#series='unknown')

            sess_reg.lines.t['x'] = (1+sess_center.systs.t['z'])\
                                    *xem_d['Ly_a'].to(sess_reg.spec.x.unit)
            sess_reg.lines.t['logN'] = sess_center.systs.t['logN']

            #sess_reg.add_syst_from_lines(series='SiII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='SiIV', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            sess_reg.add_syst_from_lines(series='CIV', logN=None, b=20.0,
                                         dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='FeII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)
            #sess_reg.add_syst_from_lines(series='MgII', logN=None, b=20.0,
            #                             dz=5e-5, z_end=zem, maxfev=10)

            sess_reg.add_syst_from_resids(chi2r_thres=2.0, logN=13.0, b=10.0,
                                          maxfev=10)

            sess_reg.compl_syst(n=10)#, z_start=2.128, z_end=2.1372)
            sess_reg.add_syst_slide(col='deabs')#, z_start=1.6, z_end=1.61)
            sess_reg.syst_merge()
            self._gui._refresh()
            sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_'
                          +datetime.date.today().isoformat()+'.xxx')
            sess_reg.save('/data/cupani/CIV/analyzed/'+t+'_latest.xxx')
            time_end = datetime.datetime.now()
            print("%s; computation time: %s" \
                  % (datetime.datetime.now(), time_end-time_start))
            """
