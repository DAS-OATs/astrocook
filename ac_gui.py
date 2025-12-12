import logging
import numpy as np
import sys
import traceback
try:
    import wx
except ImportError:
    wx = None

def main():
    if wx is None:
        raise ImportError("You are trying to run a legacy feature! Please install wxPython first.")
    
    # Must be extended path
    app = wx.App(False)
    try:
        pw = np.where([a[0]!='-' for a in sys.argv[1:]])
        fw = np.where([a[0]=='-' for a in sys.argv[1:]])
        paths = list(np.array(sys.argv[1:])[pw])
        flags = list(np.array(sys.argv[1:])[fw])

        import astrocook.settings as settings 
        if '--legacy' in flags:
            settings.MODE = 'V1'
            logging.info("Astrocook starting in V1 LEGACY mode.")
        else:
            settings.MODE = 'V2'
            logging.info("Astrocook starting in V2 ARCHITECTURE mode.")

        from astrocook.legacy.gui import GUI
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
