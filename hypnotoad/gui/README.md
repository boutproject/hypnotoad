The main functionality is implemented in `gui.py`. The windows are defined in the
auto-generated `hypnotoad_mainWindow.py` and `hypnotoad_preferences.py` files, and
hand-coded `matplotlib_widget.py`. The executable `hypnotoad-gui` is an 'entry-point'
defined in `setup.py`, implemented in the `__main__.py` file.

The `hypnotoad_mainWindow.py` file is auto-generated from the `hypnotoad_mainWindow.ui`
file. `hypnotoad_mainWindow.ui` is produced by, and can be edited by, the `qtcreator`
program, which should be available from your Linux distro's package manager.
```
$ qtcreator hypnotoad_mainWindow.ui
```

To generate the hypnotoad_mainWindow.py file, need pyside2 installed as well as
Qt.py. Run:

    $ pyside2-uic hypnotoad_mainWindow.ui -o hypnotoad_mainWindow.py
    $ python -m Qt --convert hypnotoad_mainWindow.py

The second step converts from a pyside2-specific file to one using Qt.py which
can run with pyside, pyside2, PyQt4 or PyQt5.

Similarly the `hypnotoad_preferences.py` is generated from the
`hypnotoad_preferences.ui` file.
