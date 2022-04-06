# SEISGEN
![SEISGEN](https://github.com/Liang-Ding/seisgen/blob/main/doc/figs/seisgen.png)

The SEISGEN is a python package to acquire Green's function and generate synthetic seismic waveform from the stored receiver-side 3D Strain Green's Tensor (SGT) database.
The SGT database is created by using the [SPECFEM3D_Cartesian software](https://geodynamics.org/resources/specfem3dcartesian) and [the python script](https://github.com/Liang-Ding/pyCAPLunar/blob/master/DSEM_Utils/merge_strainfield.py) in the project [pyCAPLunar](https://github.com/Liang-Ding/pyCAPLunar).

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
* The strain Green's tensor (SGT) database
* The 3D background model utilized to create the SGT database 
* The HDF5 file storing the locations and interpolation parameters of the points within the 3D model. 

## Example of usage
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

