class CookbookTemplates(object):
    """ Cookbook of spectral templates
    """

    def __init__(self):
        super(CookbookTemplates, self).__init__()

    def bb(self, temp=6000, scale=1.0):
        """ @brief Blackbody template
        @details Create a blackbody template of a given temperature
        The wavelength range used for sliding is defined in velocity units.
        @param temp Temperature (K)
        @param scale Scale factor
        @return 0
        """

        try:
            temp = float(temp)
            scale = float(scale)
        except:
            logging.error(msg_param_fail)
            return 0

        self.sess.spec._bb_template(temp, scale)
        return 0
