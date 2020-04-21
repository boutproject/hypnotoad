import setuptools

import os.path

with open("README.md", "r") as fh:
    long_description = fh.read()

version_dict = {}
with open(os.path.join("hypnotoad2", "__version__.py")) as fh:
    exec(fh.read(), version_dict)
version = version_dict["__version__"]

setuptools.setup(
    name="hypnotoad2",
    version=version,
    author="John Omotani, Ben Dudson and the BOUT++ team",
    author_email="john.omotani@ukaea.uk",
    description="Grid generator for BOUT++",
    license="OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/boutproject/hypnotoad2",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "boututils~=0.1.4",
        "matplotlib~=3.2",
        "netCDF4~=1.5",
        "numpy~=1.18",
        "options~=1.4",
        "pyparsing~=2.4",
        "PyYAML~=5.1",
        "scipy~=1.4",
        "Qt.py~=1.2",
    ],
    extras_require={
        "gui-pyside2": ["pyside2~=5.14"],
        "gui-PyQt5": ["PyQt5~=5.12"],
    },
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "hypnotoad2_geqdsk = hypnotoad2.scripts.hypnotoad2_geqdsk:main",
            "hypnotoad2_torpex = hypnotoad2.scripts.hypnotoad2_torpex:main",
        ],
        "gui_scripts": [
            "hypnotoad2-gui = hypnotoad2.gui:main",
        ],
    },
)
