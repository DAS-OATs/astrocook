from .functions import *
from .message import *
from .vars import *
from astropy import units as au
import logging

class CookbookEdit(object):
    """ Cookbook of utilities for editing data
    """

    def __init__(self):
        super(CookbookEdit, self).__init__()


    def _struct_parse(self, struct, length=2):
        sess_list = self.sess._gui._sess_list

        parse = struct.split(',')

        if len(parse) < length:
            logging.error("I can't parse the structure.")
            return None

        # Session
        sessn = parse[0]
        try:
            sessn = int(sessn)
            parse[0] = sessn
        except:
            logging.error(msg_param_fail)
            return None

        if sessn > len(sess_list):
            logging.error("I can't find session %s." % sessn)
            return None

        # Attribute
        attrn = parse[1]
        sess = sess_list[sessn]
        if not hasattr(sess, attrn):
            logging.error(msg_attr_miss(attrn))
            return None

        attr = getattr(sess, attrn)
        if attr is None:
            logging.warning("Attribute %s is None." % attrn)
            return attrn, attr, parse


        if length==3:
            # Column
            coln = parse[2]
            if coln not in getattr(sess, attrn)._t.colnames:
                logging.error(msg_col_miss(coln))
                return None
            col = getattr(sess, attrn)._t[coln]
            return coln, col, parse
        else:
            return attrn, attr, parse


    def shift_bary(self, v=None):
        """ @brief Shift to barycentric frame
        @details Shift x axis to the barycentric frame of the solar system.
        @param v Velocity in the barycentric frame (km/s)
        @return 0
        """

        try:
            v = float(v)
        except:
            try:
                v = self.sess.spec.meta['v_bary']
            except ValueError:
                logging.error(msg_param_fail)
                return 0

        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._shift_bary(v)
            except:
                logging.debug(msg_attr_miss(s))
        return 0


    def shift_from_rf(self, z=0):
        """ @brief Shift from rest frame
        @details Shift x axis from rest frame to the original frame.
        @param z Redshift to use for shifting
        @return 0
        """

        try:
            z = float(z)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        for s in self.sess.seq:
            try:
                z_to = z-getattr(self.sess, s)._rfz
                getattr(self.sess, s)._shift_rf(z_to)
            except:
                logging.debug(msg_attr_miss(s))
        return 0


    def shift_to_rf(self, z=0):
        """ @brief Shift to rest frame
        @details Shift x axis to the rest frame.
        @param z Redshift to use for shifting
        @return 0
        """

        try:
            z = float(z)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._shift_rf(z)
            except:
                logging.debug(msg_attr_miss(s))
        return 0


    def struct_import(self, struct='0,systs', mode='replace'):
        """ @brief Import a data structure from a session into the current one
        @details The structure to be imported is described by a string with the
        session number and the structure tag (spec, lines, systs), separated by
        a comma (e.g. 0,spec, meaning "spectrum from session 0"). The imported
        structure is either replaced or appended to the corresponding one in the
        current session.
        @param struct Structure
        @param mode Mode (replace or append)
        @return 0
        """

        parse = self._struct_parse(struct)
        if parse is None: return 0
        attrn, attr, _ = parse
        attr = dc(attr)
        if attrn == 'systs' \
            and 'cont' not in self.sess._gui._sess_sel.spec.t.colnames:
            logging.error("Attribute %s requires a continuum. Please try "
                          "Recipes > Guess continuum before." % attrn)
            return 0

        if mode=='replace':
            if attr is None:
                logging.warning("I'm replacing structure with None.")
                setattr(self.sess._gui._sess_sel, attrn, attr)
                return 0

            if attrn in ['lines', 'systs']:
                #spec = self.sess._gui._sess_sel.spec
                x = self.sess._gui._sess_sel.spec.x.to(au.nm)
                attr = attr._region_extract(np.min(x), np.max(x))
            if attr is None:
                logging.warning("I'm replacing structure with None")
                setattr(self.sess._gui._sess_sel, attrn, attr)
                return 0

                # Redefine regions from spectrum
            if attrn == 'systs':
                for m in attr._mods_t:
                    mod = m['mod']
                    mod._spec = self.sess._gui._sess_sel.spec
                    mod._xf, mod._yf, mod._wf = \
                        mod._make_regions(mod, mod._spec._safe(mod._spec.x)\
                                               .to(au.nm).value)
            setattr(self.sess._gui._sess_sel, attrn, attr)

        if mode=='append':
            if attr is None:
                logging.warning("I'm not appending None.")
                return 0
            attr_dc = dc(attr)
            if attrn == 'systs':
                id_max = np.max(getattr(self.sess._gui._sess_sel, attrn)._t['id'])
                attr_dc._t['id'] = attr_dc._t['id']+id_max
            #print(len(attr_dc._t))
            #print(len(np.unique(attr_dc._t['id'])))
            #print(len(attr_dc._t))
            #print(len(getattr(self.sess._gui._sess_sel, attrn)._t))
            getattr(self.sess._gui._sess_sel, attrn)._append(attr_dc)

        if attrn=='systs':
            self.sess._gui._sess_sel.cb._mods_recreate()
            self.sess._gui._sess_sel.cb._spec_update()

        return 0


    def struct_modify(self, col='', expr=''):
        """ @brief Modify a data structure using a binary operator
        @details Modify a data structure using a binary operator. An output
        column is computed from an expression with input columns as arguments.
        The expression must be parsable by AST, with columns described by a
        string with the session number, the structure tag (spec, lines, systs),
        and the column name separated by a comma (e.g. 0,spec,x, meaning "column
        x of spectrum from session 0"). Columns can be from different data
        structures only if they have the same length. If the output column
        already exists, it is overwritten.
        @param col Output column
        @param expr Expression
        @return 0
        """

        # Reversed to parse sessions with higher number first, and avoid overwriting
        sess_list = self.sess._gui._sess_list[::-1]
        for i, s in enumerate(sess_list):
            if s.spec is not None:
                for c in sorted(s.spec._t.colnames, key=len, reverse=True):
                    expr = expr.replace('%i,spec,%s' % (len(sess_list)-1-i, c),
                                        str(list(np.array(s.spec._t[c]))))
            if s.lines is not None:
                for c in sorted(s.lines._t.colnames, key=len, reverse=True):
                    expr = expr.replace('%i,lines,%s' % (len(sess_list)-1-i, c),
                                        str(list(np.array(s.lines._t[c]))))
            if s.systs is not None:
                for c in sorted(s.systs._t.colnames, key=len, reverse=True):
                    expr = expr.replace('%i,systs,%s' % (len(sess_list)-1-i, c),
                                        str(list(np.array(s.systs._t[c]))))

        #print(expr)
        #print(len(expr))
        out = expr_eval(ast.parse(expr, mode='eval').body)

        _, _, all_out = self._struct_parse(col, length=2)
        struct = getattr(self.sess._gui._sess_list[all_out[0]], all_out[1])
        if all_out[2] in struct._t.colnames: # and False:
            col_out = struct._t[all_out[2]]
            try:
                struct._t[all_out[2]] = expr_eval(ast.parse(expr, mode='eval').body) * col_out.unit
            except:
                struct._t[all_out[2]] = expr_eval(ast.parse(expr, mode='eval').body)
        else:
            struct._t[all_out[2]] = expr_eval(ast.parse(expr, mode='eval').body)

        return 0


    def x_convert(self, zem=0, xunit=au.km/au.s):
        """ @brief Convert x axis
        @details Convert the x axis to wavelength or velocity units. The x axis
        can be converted to any unit of wavelength or velocity (default: nm and
        km/s). The conversion applies to both the spectrum and the line list.
        When converting to and from velocity units, the zero point is set at
        (1+zem)λ_Lya (where λ_Lya = 121.567 nm is the rest-frame wavelength
        of the Lyman-alpha transition).
        @param zem Emission redshift
        @param xunit Unit of wavelength or velocity
        @return 0
        """

        try:
            zem = float(zem)
        except ValueError:
            logging.error(msg_param_fail)
            return 0

        xunit = au.Unit(xunit)
        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._x_convert(zem, xunit)
            except:
                logging.debug(msg_attr_miss(s))
        return 0


    def y_convert(self, yunit=au.electron/au.nm):
        """ @brief Convert y axis
        @details Convert the y axis to different units. The y axis can be
        expressed in different units depending on how it was calibrated
        (default: erg/(cm^2 s nm)). It can be converted to any unit of the same
        physical quantity. The conversion applies to both the spectrum and the
        line list.
        @param yunit Unit of flux density
        @return 0
        """

        yunit = au.Unit(yunit)

        for s in self.sess.seq:
            try:
                getattr(self.sess, s)._y_convert(yunit=yunit)
            except:
                logging.debug(msg_attr_miss(s))
        return 0
