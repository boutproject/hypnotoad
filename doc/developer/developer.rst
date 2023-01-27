Developing hypnotoad
====================

When installing from the source code git repo for development, it is
recommended to do an editable install ::

    $ pip install --user -e .

This method means the installed modules and command line tools link to the
source code, so reflect changes you make, and avoids a possible error if the
repo has uncommitted changes (is 'dirty') where ``git diff`` fails to run on
the installed code. Note, if you use ``conda`` make sure all the dependencies
(see ``pyproject.toml``) are already installed to avoid them being
``pip``-installed and messing up dependency tracking.

.. toctree::

   equilibrium
   meshregion
   gui
   RELEASE_HOWTO
