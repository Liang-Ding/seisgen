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
