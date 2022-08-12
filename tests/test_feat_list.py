import astropy.table as at
from copy import deepcopy as dc
import numpy as np
import pytest
import wx

app = wx.App(False)
from astrocook.gui import GUI
gui = GUI(['tests/test_feat_list.acs'])

from astrocook.feat_list import Feat, FeatList, FeatTable
from astrocook.spectrum import Spectrum
spec = gui._sess_sel.spec
chunk = spec._t[78:105]
systs = gui._sess_sel.systs
feats = FeatList()
feat = Feat(chunk)
featt = FeatTable()


class TestFeat:

    def test__ccf_compute(self):
        xmean = 579
        deltav = 10
        feat._ccf_compute(xmean, deltav)
        assert feat._ccf_xmean == xmean
        assert feat._ccf_deltav == deltav


    def test__ew_compute(self):
        feat._ew_compute()
        assert feat._ew.value == pytest.approx(0.017371824670513115)
        assert feat._dew.value == pytest.approx(0.0002772926552396542)
        assert str(feat._ew.unit) == 'nm'
        assert str(feat._dew.unit) == 'nm'


    def test__fwhm_compute(self):
        feat._fwhm_compute()
        assert feat._fwhm.value == pytest.approx(0.05015866415169512)
        assert str(feat._fwhm.unit) == 'nm'


    def test__logN_compute(self):
        feat._logN_compute()
        assert not hasattr(feat, '_logN')

        feat_copy = dc(feat)
        feat_copy._systs_join(systs)
        feat_copy._logN_compute()
        assert feat_copy._logN == pytest.approx(13.103592543222868)


    def test__snr_compute(self):
        feat._snr_compute()
        assert feat._snr == pytest.approx(60.96067414495566)


    def test__systs_check(self):
        assert feat._systs_check() == 0

        feat_copy = dc(feat)
        feat_copy._systs_join(systs)
        assert feat_copy._systs_check() == 1


    def test__systs_join(self):
        feat_copy = dc(feat)
        feat_copy._systs_join(systs)
        assert feat_copy._systs_orig == systs
        assert feat_copy._systs == {k: systs._d[k] for k in (40, 41, 42) if k in systs._d}
        assert feat_copy._trans == {40: 'CIV_1548', 41: 'CIV_1548', 42: 'CIV_1548'}


    def test__systs_logN_tot(self):
        feat_copy = dc(feat)
        feat_copy._systs_join(systs)
        feat_copy._systs_logN_tot()
        assert feat_copy._logN_tot is None
        assert feat_copy._dlogN_tot is None

        feats_copy = dc(feats)
        systs_copy = dc(systs)
        feats_copy.create(spec, systs_copy, thres=1e-2)
        feats_copy._logN_tot(systs_copy, set_specs=False)
        feats_copy._l[0]._logN_tot_par.value = 14
        feats_copy._l[0]._logN_tot_par.stderr = 1
        feats_copy._l[0]._systs_logN_tot()
        assert feats_copy._l[0]._logN_tot == pytest.approx(14)
        assert feats_copy._l[0]._dlogN_tot == pytest.approx(1)

        feats_copy._l[0]._N_tot_par = dc(feats_copy._l[0]._logN_tot_par)
        feats_copy._l[0]._systs_logN_tot()
        assert feats_copy._l[0]._logN_tot == pytest.approx(1.146128035678238)
        assert feats_copy._l[0]._dlogN_tot == pytest.approx(0.031073953374422217)
        feats_copy._l[0]._N_tot_par.stderr = None
        feats_copy._l[0]._systs_logN_tot()
        assert feat_copy._dlogN_tot is None


    def test__systs_stats(self):
        self.test__logN_compute()
        self.test__xz_compute()
        self.test__ew_compute()
        self.test__fwhm_compute()
        self.test__snr_compute()


    def test__xz_compute(self):
        feat._xz_compute()
        assert not hasattr(feat, '_z')
        assert not hasattr(feat, '_x')

        feat_copy = dc(feat)
        systs_copy = dc(systs)
        feat_copy._systs_join(systs_copy)
        feat_copy._systs = {}
        feat_copy._xz_compute()
        assert not hasattr(feat_copy, '_z')
        assert not hasattr(feat_copy, '_x')

        feat_copy = dc(feat)
        systs_copy = dc(systs)
        feat_copy._systs_join(systs_copy)
        for s in feat_copy._systs.values():
            s._pars['logN'] = 0
        feat_copy._xz_compute()
        assert feat_copy._z == pytest.approx(2.735670073315127)
        assert feat_copy._x.value == pytest.approx(578.3579350186773)
        assert str(feat_copy._x.unit) == 'nm'

        feat_copy = dc(feat)
        feat_copy._systs_join(systs)
        feat_copy._xz_compute()
        assert feat_copy._z == pytest.approx(2.735755719507012)
        assert feat_copy._x.value == pytest.approx(578.3711947963634)
        assert str(feat_copy._x.unit) == 'nm'


class TestFeatList():


    #def test_t(self):
    #    assert feats.t == feats._t

    def test__add(self):
        feats_copy = dc(feats)
        feats_copy._add(chunk, systs)
        assert str(feats_copy._xunit) == 'nm'
        assert str(feats_copy._yunit) == 'erg / (Angstrom cm2 s)'
        assert len(feats_copy._l) == 1


    def test__check_attr(self):
        assert feats._check_attr('none') is None
        assert feats._check_attr('_l') == feats._l


    #def test__load(self):

    def test__logN_tot(self):
        feats_copy = dc(feats)
        systs_copy = dc(systs)
        feats_copy.create(spec, systs_copy, thres=1e-2)
        feats_copy._logN_tot(systs_copy)
        assert systs_copy._N_tot_specs == {
            45: ('lines_voigt_45_', 16672607372083.424,
                 '+10**lines_voigt_44_logN+10**lines_voigt_43_logN+10**lines_'\
                 'voigt_42_logN+10**lines_voigt_41_logN+10**lines_voigt_40_logN'),
            48: ('lines_voigt_48_', 1825561191171.5479, '')}

        feats_copy = dc(feats)
        systs_copy = dc(systs)
        feats_copy.create(spec, systs_copy, thres=1e-2)
        feats_copy._logN_tot(systs_copy, sel='0:1')
        assert systs_copy._N_tot_specs == {
            45: ('lines_voigt_45_', 16672607372083.424,
                 '+10**lines_voigt_44_logN+10**lines_voigt_43_logN+10**lines_'\
                 'voigt_42_logN+10**lines_voigt_41_logN+10**lines_voigt_40_logN')}

        feats_copy = dc(feats)
        systs_copy = dc(systs)
        feats_copy.create(spec, systs_copy, thres=1e-2)
        feats_copy._logN_tot(systs_copy, sel=0)
        assert systs_copy._N_tot_specs == {}

        feats_copy = dc(feats)
        systs_copy = dc(systs)
        feats_copy.create(spec, systs_copy, thres=1e-2)
        feats_copy._l[0]._systs = {}
        feats_copy._logN_tot(systs_copy)
        assert systs_copy._N_tot_specs == {}

        feats_copy = dc(feats)
        systs_copy = dc(systs)
        feats_copy.create(spec, systs_copy, thres=1e-2)
        feats_copy._l[0]._systs = {40: feats_copy._l[0]._systs[40]}
        feats_copy._logN_tot(systs_copy)
        assert systs_copy._N_tot_specs == {
            48: ('lines_voigt_48_', 1825561191171.5479, '')}

        feats_copy = dc(feats)
        systs_copy = dc(systs)
        feats_copy.create(spec, systs_copy, thres=1e-2)
        s = feats_copy._l[0]._systs
        s = {40: feats_copy._l[0]._systs[40]}
        s[max(s.keys())]._mod._pars.add("dummy_par_0", 1, True, 0, 2, None)
        feats_copy._logN_tot(systs_copy)
        assert systs_copy._N_tot_specs == {
            45: ('lines_voigt_45_', 16672607372083.424,
                 '+10**lines_voigt_44_logN+10**lines_voigt_43_logN+10**lines_'\
                 'voigt_42_logN+10**lines_voigt_41_logN+10**lines_voigt_40_logN'),
            48: ('lines_voigt_48_', 1825561191171.5479, '')}

        feats_copy = dc(feats)
        systs_copy = dc(systs)
        feats_copy.create(spec, systs_copy, thres=1e-2)
        feats_copy._logN_tot(systs_copy, set_specs=False)
        assert systs_copy._N_tot_specs == {}


    def test__maxs_from_spec(self):
        feats_copy = dc(feats)
        feats_copy._maxs_from_spec(spec)
        assert (feats_copy._maxs == np.asarray([0, 168, 273, -1])).all

        feats_copy = dc(feats)
        feats_copy._maxs_from_spec(spec, height=1)
        assert (feats_copy._maxs == np.asarray([0, -1])).all

        feats_copy = dc(feats)
        feats_copy._maxs_from_spec(spec, prominence=1)
        assert (feats_copy._maxs == np.asarray([0, -1])).all


    #def test__save(self):


    def test__select_isolated(self):
        feats._select_isolated()


    def test__systs_update(self):
        feats_copy = dc(feats)
        feats_copy.create(spec, systs, thres=1e-2)
        feats_copy._systs_update(systs)
        systs_dicts = [{k: systs._d[k] for k in (40,41,42,43,44,45) \
                        if k in systs._d},
                       {k: systs._d[k] for k in (46,47,48) if k in systs._d},
                       {k: systs._d[k] for k in (40,41,42,43,44,45,46,47,48) \
                        if k in systs._d}]
        trans_dicts = [{40: 'CIV_1548', 41: 'CIV_1548', 42: 'CIV_1548', \
                        43: 'CIV_1548', 44: 'CIV_1548', 45: 'CIV_1548'},
                       {46: 'CIV_1548', 47: 'CIV_1548', 48: 'CIV_1548'},
                       {40: 'CIV_1550', 41: 'CIV_1550', 42: 'CIV_1550', \
                        43: 'CIV_1550', 44: 'CIV_1550', 45: 'CIV_1550', \
                        46: 'CIV_1550', 47: 'CIV_1550', 48: 'CIV_1550'},]

        logNs = [13.222003522968974, 12.598483068140343, 13.31470544234731]
        zs = [2.736307425275211, 2.7386654915260893, 2.7384028275005026]
        xs = [578.4566101040783, 578.8216868642659, 579.7444075234056]
        ews = [0.022028967649447246, 0.0055587149113473285, 0.014897938451949692]
        dews = [0.0003734132873864521, 0.0003185964270165629, 0.0005032446499793494]
        fwhms = [0.05015866415169512, 0.08109204866002528, 0.05024205441884533]
        snrs = [62.81903479311659, 65.88475838204226, 62.986046902136756]

        for i, f in enumerate(feats_copy._l):
            fc = dc(f)
            fc._systs_join(systs)
            assert fc._systs_orig == systs
            assert fc._systs == systs_dicts[i]
            assert fc._trans == trans_dicts[i]
            assert fc._logN == pytest.approx(logNs[i])
            assert fc._x.value == pytest.approx(xs[i])
            assert fc._ew.value == pytest.approx(ews[i])
            assert fc._dew.value == pytest.approx(dews[i])
            assert fc._fwhm.value == pytest.approx(fwhms[i])
            assert fc._snr == pytest.approx(snrs[i])
            assert str(fc._x.unit) == 'nm'
            assert str(fc._ew.unit) == 'nm'
            assert str(fc._dew.unit) == 'nm'
            assert str(fc._fwhm.unit) == 'nm'


    def test__table_update(self):
        feats_copy = dc(feats)
        feats_copy.create(spec, systs, thres=1e-2)
        feats_copy._table_update()
        table = at.Table()
        x = [578.301345706367, 578.4788414597525, 578.7722146785759,
             578.8919227872219, 579.266653288144, 579.8504794237512]
        m = [42987.088838146, 42638.41469010227, 42111.62443671379,
             41854.64352380562, 40694.39651763134, 38941.287827606626]
        c = [43166.573550840425, 42836.98033069381, 42200.45411706748,
             41903.61326477187, 40843.362917297294, 39024.6883197286]
        table['x'] = at.Column(np.array(x, ndmin=1), dtype=float, unit='nm')
        table['model'] = at.Column(np.array(m, ndmin=1), dtype=float, unit='erg / (Angstrom cm2 s)')
        table['cont'] = at.Column(np.array(c, ndmin=1), dtype=float, unit='erg / (Angstrom cm2 s)')
        assert (feats_copy._t['x'] == table['x']).any
        assert (feats_copy._t['model'] == table['model']).any
        assert (feats_copy._t['cont'] == table['cont']).any


    def test__z_lock(self):
        feats_copy = dc(feats)
        feats_copy.create(spec, systs, thres=1e-2)
        feats_copy._table_update()
        feats_copy._systs_update(systs)
        feats_copy._z_lock()
        constr = {'lines_voigt_40_b': (40, 'b', None),
                  'lines_voigt_41_z': (41, 'z', 'lines_voigt_40_z+0.00005949597072'),
                  'lines_voigt_41_b': (41, 'b', None),
                  'lines_voigt_42_z': (42, 'z', 'lines_voigt_40_z+0.00021225759506'),
                  'lines_voigt_42_b': (42, 'b', None),
                  'lines_voigt_43_z': (43, 'z', 'lines_voigt_40_z+0.00041601458070'),
                  'lines_voigt_43_b': (43, 'b', None),
                  'lines_voigt_44_z': (44, 'z', 'lines_voigt_40_z+0.00060045339307'),
                  'lines_voigt_44_b': (44, 'b', None),
                  'lines_voigt_45_z': (45, 'z', 'lines_voigt_40_z+0.00074250385575'),
                  'lines_voigt_45_b': (45, 'b', None),
                  'lines_voigt_46_b': (46, 'b', None),
                  'lines_voigt_47_z': (47, 'z', 'lines_voigt_46_z+-0.00006131361621'),
                  'lines_voigt_47_b': (47, 'b', None)}
        assert feats_copy._l[0]._systs_orig._constr == constr


    def test_create(self):
        self.test__maxs_from_spec()
        self.test__table_update()

        feats_copy = dc(feats)
        feats_copy.create(spec, systs, thres=1)
        table = at.Table()
        x = []
        m = []
        c = []
        table['x'] = at.Column(np.array(x, ndmin=1), dtype=float, unit='nm')
        table['model'] = at.Column(np.array(m, ndmin=1), dtype=float, unit='erg / (Angstrom cm2 s)')
        table['cont'] = at.Column(np.array(c, ndmin=1), dtype=float, unit='erg / (Angstrom cm2 s)')
        assert (feats_copy._t['x'] == table['x']).any
        assert (feats_copy._t['model'] == table['model']).any
        assert (feats_copy._t['cont'] == table['cont']).any
