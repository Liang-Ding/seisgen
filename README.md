# SEISGEN
![SEISGEN](https://github.com/Liang-Ding/seisgen/blob/main/doc/figs/seisgen.png)

SEISGEN is a python package to acquire and generate the Greens function and synthetic waveform from the stored receiver-side 3D database including the Strain Greens Tensor (SGT) database and the (displacement) Greens function (DGF) database. 
The Greens function and synthetic waveform could be utilized by inversion packages such as [pyCAPLunar](https://github.com/Liang-Ding/pyCAPLunar), MTUQ and gCAP-series packages to determine the parameters of single-point sources including the moment tensor and force and finite fault models. 

The 3D Greens function database is created by using the [SPECFEM3D_Cartesian software](https://geodynamics.org/resources/specfem3dcartesian) and [the python script](https://github.com/Liang-Ding/pyCAPLunar/blob/master/DSEM_Utils/merge_strainfield.py) in the project [pyCAPLunar](https://github.com/Liang-Ding/pyCAPLunar).

## Installation
For basic install:
```shell
git clone https://github.com/Liang-Ding/seisgen.git
cd seisgen
pip install -e .
```
or using pip 
```shell
pip install seisgen
```

## Package Structure
![SEISGEN](https://github.com/Liang-Ding/seisgen/blob/main/doc/figs/seisgen_structure.png)
 
The seisgen package requires the following database to work:
* The stored Greens function database, such as the SGT, DGF database
* The 3D background model utilized to create the SGT database 
* The HDF5 file storing the locations and interpolation parameters of the points within the 3D model. 

## Example of using the SGT database
* Initialize an instance of the data manager:
```shell
mgr = DSGTMgr(sgt_database_folder, model3D_folder, point_cloud_file)
```

* Acquire the Green's function in RTZ
```shell
greens = mgr.get_greens_function(station, origin)
```

* Acquire the F-K type greens function
```shell
greens = mgr.get_fk_greens_function(station, origin)
```

* Acquire the synthetic waveform in RTZ
```shell
stream = mgr.get_waveform(station, origin, mt_RTP)
```

## Example of using the DGF database
* Initialize an instance of the data manager:
```shell
mgr = DDGFMgr(sgt_database_folder, model3D_folder, point_cloud_file)
```

* Acquire the Green's function in ENZ
```shell
greens = mgr.get_greens_function(station, origin)
```

* Acquire the synthetic waveform in ENZ
```shell
stream = mgr.get_waveform(station, origin, force_enz)
# or in RTZ
stream = mgr.get_waveform(station, origin, force_enz, b_RTZ=True)
```
