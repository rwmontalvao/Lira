# Lira
Lira: Rotational Invariant Shape and ElectrostaticDescriptors for Small Molecules and ProteinPockets based on Real Spherical Harmonics

### Prerequisites

1. Install [DUB package manager](https://github.com/dlang/dub/releases): `curl -fsS https://dlang.org/install.sh | bash -s ldc`


### Build

1. `git clone git@github.com:rwmontalvao/Lira.git`
1. `cd Lira/lira_coeffs`
1. `source ~/dlang/ldc-1.30.0/activate`
1. `dub --compiler=ldmd2 --build=release`
1. `deactivate`

### Test

1. `lira_coeffs --in_file ZINC67007472.tmesh --out_file /tmp/test.coeff --verbose`
1. `diff ZINC67007472.coeffs /tmp/test.coeff`
1. `rm /tmp/test.coeff`
