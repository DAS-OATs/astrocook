import logging
import numpy as np
import sys
import wx

def main():

    # Must be extended path
    app = wx.App(False)
    from astrocook.gui import GUI
    try:
        ok
    except:
        pw = np.where([a[0]!='-' for a in sys.argv[1:]])
        fw = np.where([a[0]=='-' for a in sys.argv[1:]])
        paths = list(np.array(sys.argv[1:])[pw])
        flags = list(np.array(sys.argv[1:])[fw])
        gui = GUI(paths, flags)
        ok = gui._ok
    #except:
    #    logging.error("I found some problems loading this session.")
    #    ok = False

    if not ok:
        logging.warning("Re-starting with an empty session.")
        gui = GUI()
    app.MainLoop()


if __name__ == '__main__':
    main()
