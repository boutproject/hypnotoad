To generate a single null tokamak equilibrium:

$ ./tokamak_example.py single-null.yaml

A lower double null equilibrium:

$ ./tokamak_example.py lower-double-null.yaml

Those inputs can be modified to generate additional inputs.
The analytic poloidal flux functions can be chosen with
the geometry setting:

"lsn"   Lower Single Null
"usn"   Upper Single Null
"cdn"   Connected Double Null
"udn"   Upper Double Null
"ldn"   Lower Double Null
"udn2"  Upper Double Null, with a larger gap between separatrices.

