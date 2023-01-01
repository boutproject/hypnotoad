Other configurations
====================

Command line utilities are provided for a couple of non-standard configurations.

Circular geometry
-----------------

``hypnotoad-circular`` provides an interface for creating a grid file for a
configuration with concentric, circular flux surfaces. Input is provided by a
YAML settings file, passed as a command line argument. The grid can be 'core
only' (only closed flux surfaces) when ``limiter = False`` or 'SOL only' (field
lines end on a limiter, no closed flux surfaces) when ``limiter = True``.

The full set of options are listed here: :ref:`_temp/circular-options:Circular
options`.
Nonorthogonal options can also be used, but are unlikely to be helpful as any
limiter present is orthogonal to the flux surfaces.

TORPEX X-point
--------------

``hypnotoad-torpex`` provides an interface for creating a grid file for a
configuration with an isolated X-point, where all four separatrix branches end
on a wall so there are no closed flux surfaces. A configuration of this type
was studied in the TORPEX device.

The full set of options are listed here: :ref:`_temp/torpex-options:TORPEX
options`. Note that at least ``psi_core`` and ``psi_sol`` are required to be
set (some values whose defaults say '*Required*' in the table will have default
values set from those two).

Nonorthogonal options can also be used:
:ref:`_temp/nonorthogonal-options:Nonorthogonal options`.
