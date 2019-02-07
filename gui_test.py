import wx

def main():

    # Must be extended path
    app = wx.App(False)
    from astrocook.gui import GUI
    gui = GUI()
    app.MainLoop()

if __name__ == '__main__':
    main()
