import logging
import numpy as np
import sys
import traceback
#import wx

def main():

    # Must be extended path
    app = wx.App(False)
    from astrocook.gui import GUI
    try:
        pw = np.where([a[0]!='-' for a in sys.argv[1:]])
        fw = np.where([a[0]=='-' for a in sys.argv[1:]])
        paths = list(np.array(sys.argv[1:])[pw])
        flags = list(np.array(sys.argv[1:])[fw])
        gui = GUI(paths, flags)
        tb = None
    except:
        logging.error("I found some problems loading this session. Here's the "\
                      "traceback:\n")
        tb = traceback.format_exc()

    if tb is not None:
        print(tb)
        logging.warning("Re-starting with an empty session.")
        gui = GUI()
    app.MainLoop()


if __name__ == '__main__':
    main()
