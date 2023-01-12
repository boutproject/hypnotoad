import setuptools
import versioneer

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hypnotoad",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="John Omotani, Ben Dudson and the BOUT++ team",
    author_email="john.omotani@ukaea.uk",
    description="Grid generator for BOUT++",
    license="OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/boutproject/hypnotoad",
    packages=setuptools.find_packages(include=("hypnotoad*",)),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "boututils~=0.1.7",
        "dill~=0.3,!=0.3.5,!=0.3.5.1",
        "func_timeout~=4.3",
        "matplotlib~=3.2",
        "netCDF4~=1.5",
        "numpy~=1.18",
        "optionsfactory~=1.0.11",
        "PyYAML>=5.1",
        "scipy~=1.6",
        "Qt.py~=1.2",
    ],
    extras_require={"gui-pyside2": ["pyside2~=5.13"], "gui-PyQt5": ["PyQt5~=5.12"]},
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "hypnotoad-circular = hypnotoad.scripts.hypnotoad_circular:main",
            "hypnotoad-geqdsk = hypnotoad.scripts.hypnotoad_geqdsk:main",
            "hypnotoad-torpex = hypnotoad.scripts.hypnotoad_torpex:main",
            "hypnotoad-plot-equilibrium = "
            "hypnotoad.scripts.hypnotoad_plot_equilibrium:main",
            "hypnotoad-plot-grid-cells = "
            "hypnotoad.scripts.hypnotoad_plot_grid_cells:main",
            "hypnotoad-recreate-inputs = "
            "hypnotoad.scripts.hypnotoad_recreate_inputs:main",
        ],
        "gui_scripts": ["hypnotoad-gui = hypnotoad.gui:main"],
    },
)
