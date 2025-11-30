class CookbookViewOld(object):
    """ Cookbook of utilities for editing data
    """

    def __init__(self):
        super(CookbookViewOld, self).__init__()


    def z_ax(self, trans='Ly_a'):
        """ @brief Show redshift axis
        @details Show an additional axis on the plot with wavelength converted
        into redshift for a given transition
        @param trans Transition
        @return 0
        """
        self.sess._ztrans = trans

        return 0
