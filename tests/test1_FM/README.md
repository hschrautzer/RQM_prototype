# Test 1

This is the first simple test. The setup is:

- Square lattice (10 x 10 x 10) magnetic moments
- Interactions: Nearest neighbor exchange with $J=1$ eV.
- Magnetic configuration: Ferromagnetic aligned parallel to $\vec{e}_z$.
- Compute the lowest 4 modes using the RQM.

## Expected result

We expect two degenerate zero modes corresponding to two perpendicular global rotations of the magnetic configuration. Probably the next two modes are also degenerate.

**NOTE**: One should not cut degenerate subspaces. E.g. if we would try 3 lowest modes RQM has a problem by design. We are trying to approximate a *half broken subspace*. This should work in the end, but this is not a good first test.

## Spinaker

The result of Spinaker, used for comparison, are included in the folder `Spinaker` .

The ferromagnetic configuration is included also in this test directory and called `spin_i.dat`.

## RQM Prototype
