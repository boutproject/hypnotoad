To generate the hypnotoad_mainWindow.py file, need pyside2 installed as well as
Qt.py. Run:

    $ pyside2-uic hypnotoad_mainWindow.ui -o hypnotoad_mainWindow.py
    $ python -m Qt --convert hypnotoad_mainWindow.py

The second step converts from a pyside2-specific file to one using Qt.py which
can run with pyside, pyside2, PyQt4 or PyQt5.
