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
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "boututils~=0.1.7",
        "func_timeout~=4.3",
        "matplotlib~=3.2",
        "netCDF4~=1.5",
        "numpy~=1.18",
        "optionsfactory~=1.0.1",
        "pyparsing~=2.4",
        "PyYAML~=5.1",
        "scipy~=1.4",
        "Qt.py~=1.2",
    ],
    extras_require={"gui-pyside2": ["pyside2~=5.13"], "gui-PyQt5": ["PyQt5~=5.12"]},
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "hypnotoad_geqdsk = hypnotoad.scripts.hypnotoad_geqdsk:main",
            "hypnotoad_torpex = hypnotoad.scripts.hypnotoad_torpex:main",
        ],
        "gui_scripts": ["hypnotoad-gui = hypnotoad.gui:main"],
    },
)
