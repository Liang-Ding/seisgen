# SEISGEN
![SEISGEN](https://github.com/Liang-Ding/seisgen/blob/main/doc/figs/seisgen.png)

The SEISGEN is a python package to generate synthetic seismic waveform by utilizing the stored receiver-side Strain Green's Tensor (SGT) database that is created by using [the python script](https://github.com/Liang-Ding/pyCAPLunar/blob/master/DSEM_Utils/merge_strainfield.py) in the project [pyCAPLunar](https://github.com/Liang-Ding/pyCAPLunar) by merging the strain-filed data that is written out by the SPECFEM3D package upon waveform simulation in the 3D background model.

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

* Acquire the synthetic waveform in RTZ
```shell
stream = mgr.get_waveform(station, origin, mt_RTP)
```

