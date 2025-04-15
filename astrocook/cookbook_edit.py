from .cookbook_edit_old import CookbookEditOld
from .functions import *
from .message import *
from .vars import *

from astropy import units as au
import logging

class CookbookEdit(CookbookEditOld):
    """ Cookbook of utilities for editing data
    """

    def __init__(self):
        super(CookbookEdit, self).__init__()


    def modify_columns(self):
        """@brief Modify columns 🚧
        @details 🚧
        @url edit_cb.html#modify-columns
        """

        return 0


    def import_systs(self, source=0, mode='replace'):
        """@brief Import system list
        @details Import a system list into the currently selected session.
        @url edit_cb.html#import-system-list
        @param source Source session
        @param mode Mode (replace or append)
        @return 0
        """

        try:
            source = str(source)
        except:
            logging.error(msg_param_fail)
            return 0

        struct = source+',systs'
        return self.struct_import(struct, mode)
    
    
    def import_telluric(self, source=0, col='telluric_model', merge_cont=True):
        """@brief Import telluric model
        @details Import a telluric model into the currently selected session 
        and optionally merge it into the continuum.
        @url edit_cb.html#import-telluric-model
        @param source Source session
        @param col Telluric model column
        @param merge_cont Merge telluric model into continuum
        @return 0
        """
        
        try:
            source = str(source)
            col = str(col)
            merge_cont = str(merge_cont) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0
        
        struct = source+',spec,'+col
        
        _, model, _ = self._struct_parse(struct, length=3)
        t = self.sess.spec._t
        
        
        if len(t)!=len(model):
            logging.error("I cannot import the telluric model. The spectrum "
                          "table has a different length.")
            return 0

        logging.info("Creating column `telluric_model`...")
        t['telluric_model'] = model
        if merge_cont and 'cont' in t.colnames:
            logging.info("Merging the telluric model with the continuum...")
            if 'cont_no_telluric' not in t.colnames:
                logging.info("Creating column `cont_no_telluric` to store the" 
                             "existing continuum...")
                t['cont_no_telluric'] = t['cont']
            t['cont'] *= t['telluric_model']
        elif merge_cont and 'cont' not in t.colnames:
            logging.warning("I cannot merge the model with the continuum: "
                            "continuum not found.")
        return 0
        

    def extract(self):
        """@brief Extract 🚧
        @details 🚧
        @url edit_cb.html#extract
        """

        return 0


    def mask(self, shift=0, tell=True, sky=True, cond=''):
        """@brief Mask
        @details Mask telluric lines, sky lines, or spectral regions defined by
        a condition.
        @url edit_cb.html#mask
        @param shift Shift to the barycentric frame (km/s)
        @param tell Mask telluric lines
        @param sky Mask sky lines
        @param cond Condition
        @return 0
        """

        try:
            shift = float(shift)
            tell = str(tell) == 'True'
            sky = str(sky) == 'True'
        except:
            logging.error(msg_param_fail)
            return 0

        if tell: self.telluric_mask(shift=shift)
        if sky: self.sky_mask(shift=shift)
        if cond!='':
            self.mask_cond(cond=cond, new_sess=False)
            self.sess.spec._t['mask'] = np.logical_not(self.sess.spec._t['mask'])
            mask = np.where(np.array(self.sess.spec._t['mask']==0, dtype=bool))
            self.sess.spec._t['y'][mask] = np.nan

        return 0

    def combine(self, name='*_combined', unique=True, _sel=''):
        """ @brief Combine two or more sessions
        @details Create a new session combining the spectra from two or more
        other sessions.

        The recipe collects all the bins from the original spectra and puts them
        all together in the new spectrum. The bins retain their original size
        (defined by `xmin` and `xmax`), so they may overlap in the final
        spectrum. By default, they are ordered by ascending `x`.

        All other structures from the original sessions (line lists, etc.) are
        not propagated to the new one.

        N.B. To select sessions, either click on the session window or provide
        a list through the hidden parameter `_sel`.
        @url edit_cb.html#combine
        @param name Name of the output session
        @return Combined session
        """

        unique = str(unique) == 'True'

        name_in = name
        #sel = self._tab._get_selected_items()
        sel = self.sess._gui._sess_item_sel
        sess_list = self.sess._gui._sess_list

        """
        if isinstance(_sel, list) and _sel != []:
            sel = _sel
        if isinstance(_sel, str) and _sel != '':
            try:
                sel = [int(s) \
                       for s in _sel.replace('[','').replace(']','').split(',')]
            except:
                pass
        if sel == []:
            sel = range(len(sess_list))
        self._gui._sess_item_sel = sel
        """
        sel = _sel

        struct_out = {}
        #for struct in sess_list[sel[0]].seq:
        #    struct_out[struct] = dc(getattr(sess_list[sel[0]], struct))


        if name_in[0] == '*':
            name = sess_list[sel[0]].name

        logging.info("Combining sessions %s..." % ', '.join(str(s) for s in sel))
        for s in sel[1:]:
            #spec._t = at.vstack([spec._t, self._gui._sess_list[s].spec._t])

            for struct in sess_list[s].seq:
                if getattr(sess_list[s], struct) != None:
                    if struct_out[struct] != None:
                        struct_out[struct]._append(
                            getattr(sess_list[s], struct), unique=unique)
                    else:
                        struct_out[struct] = dc(getattr(sess_list[s], struct))


            if name_in[0] == '*':
                name += '_' + sess_list[s].name

        struct_out['spec']._t.sort('x')
        if name_in[0] == '*':
            name += name_in[1:]
        from .session import Session
        sess = Session(gui=self.sess._gui, name=name, spec=struct_out['spec'],
                       nodes=struct_out['nodes'], lines=struct_out['lines'],
                       systs=struct_out['systs'])
        return sess

    def equalize(self, xmin, xmax, _sel='', cont=True):
        """ @brief Equalize two sessions
        @details Equalize the spectrum of two sessions, based on their flux
        ratio within a wavelength window.

        By default, the last-selected spectrum is equalized to the
        first-selected one (which is left unchanged). Equalization is done in
        place, without creating a new session.

        To compute the rescaling factor, the recipe takes the medians of the `y`
        columns of the two spectra between `xmin` and `xmax`. The `y` and `dy`
        columns of the second spectrum are then multiplied by $$
        \\textrm{med}($$`y`$$_1)/\\textrm{med}($$`y`$$_2)$$.

        N.B. To select sessions, either click on the session window or provide
        a list through the hidden parameter `_sel`.
        @param xmin Minimum wavelength (nm)
        @param xmax Maximum wavelength (nm)
        @return 0
        """

        try:
            xmin = float(xmin) * au.nm
            xmax = float(xmax) * au.nm
        except ValueError:
            logging.error(msg_param_fail)
            return None

        """
        sel = self._gui._sess_item_sel
        if isinstance(_sel, list) and _sel != []:
            sel = _sel
        if isinstance(_sel, str) and _sel != '':
            sel = [int(s) \
                for s in _sel.replace('[','').replace(']','').split(',')]
        self._gui._sess_item_sel = sel
        """
        sel = _sel
        logging.info("Equalizing session %i to session %i... "
                     % (sel[1], sel[0]))

        for i,s in enumerate(sel):
            sess = self.sess._gui._sess_list[s]
            w = np.where(np.logical_and(sess.spec.x>xmin, sess.spec.x<xmax))[0]
            if len(w)==0:
                logging.error("I can't use this wavelength range for "
                              "equalization. Please choose a range covered by "
                              "both sessions.")
                return(0)
            if i == 0:
                f = np.nanmedian(sess.spec.y[w]).value
                #print(np.median(sess.spec.y[w]))
            else:
                f = f/np.nanmedian(sess.spec.y[w]).value
                #print(np.median(sess.spec.y[w]), f)
                logging.info("Equalization factor: %3.4f." % f)
                sess.spec.y = f*sess.spec.y
                sess.spec.dy = f*sess.spec.dy
                if cont and 'cont' in sess.spec._t.colnames:
                    sess.spec._t['cont'] = f*sess.spec._t['cont']

        return 0