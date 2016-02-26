Code to query the ESO web server and pull over grids of template spectra.
To use:
setup sims_photUtils, lsst_utils, and sims_skybrightness_data
Untar the all_spec.tgz file and run package_spec.py

That will overwrite the .npz files in sims_skybrightness_data with ones calculated with your current photUtils and bandpasses.

