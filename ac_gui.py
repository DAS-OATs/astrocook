import logging
import sys
import wx

def main():

    # Must be extended path
    app = wx.App(False)
    from astrocook.gui import GUI
    try:
        path = sys.argv[1]
        gui = GUI(path)
    except:
        logging.exception("I found some problems loading this session.")
        logging.warning("Re-starting with an empty session.")
        gui = GUI()
    app.MainLoop()

if __name__ == '__main__':
    main()
