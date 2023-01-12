Provenance tracking
===================

Describe outputs for provenance tracking.

.. list-table::
   :widths: 30 70

   * - ``hypnotoad_inputs``

     - String containing an plain text formatted table of the input parameters.

   * - ``hypnotoad_inputs_yaml``

     - String containing the inputs in YAML format. Can be used to recreate an
       input file that can be used to rebuild the grid, see
       :ref:`hypnotoad-recreate-inputs`.

   * - ``hypnotoad_input_geqdsk_file_contents``

     - String with a copy of the geqdsk file used to create the grid file. Can
       be extracted using :ref:`hypnotoad-recreate-inputs`.

   * - ``Python_version``

     - Python version that was used to run hypnotoad when the grid file was
       generated.

   * - ``module_versions``

     - String with a list of all the Python modules that were loaded at the
       time that hypnotoad was run, and their versions.

.. _hypnotoad-recreate-inputs:

``hypnotoad-recreate-inputs``
-----------------------------

.. argparse::
   :module: hypnotoad.scripts.hypnotoad_recreate_inputs
   :func: get_arg_parser
   :prog: hypnotoad-recreate-inputs
