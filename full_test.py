import os

def launch(command):
    print("---> "+command)
    os.system(command)

launch('python frame_test.py -v')
launch('python spectrum_test.py -v')
