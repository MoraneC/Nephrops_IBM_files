# Nephrops_IBM_files
Repository containing a few files for parametrising the larval dispersal of Norway lobster

- **Opendrift_parametrisation_lobster_x_CMEMS.py** : Part of Opendrift framework, to prepare and start the simulation of particles representing the Nephrops norvegicus larvae. Call for hydrodynamic data through Marine Copernicus APIs (It needs an account). It can be modified.
  
- **IBM_NEP_Med.txt** : Text file containing the location for release of particles in Opendrift_parametrisation_lobster_x_CMEMS.py

- **Opendrift_module_lobster.py**: Part of Opendrift framework with biological models for simulating behaviour and growth of the larvae during the larval transport. After installing Opendrift, this script is added in the opendrift/models folder.

Outcome of the larval transport simulations are stored in Results_lobster_full.RData & Results_lobster_subset.RData available in Zenodo: 10.5281/zenodo.19652873


