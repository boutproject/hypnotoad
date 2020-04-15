import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hypnotoad2",
    version="0.1",
    author="John Omotani, Ben Dudson and the BOUT++ team",
    author_email="john.omotani@ukaea.uk",
    description="Grid generator for BOUT++",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/boutproject/hypnotoad2",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL v3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "hypnotoad2_geqdsk = hypnotoad2.scripts.hypnotoad2_geqdsk:main",
            "hypnotoad2_torpex = hypnotoad2.scripts.hypnotoad2_torpex:main"
        ],
    }
)
