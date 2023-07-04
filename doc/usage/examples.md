
## Examples of using the SGT database
* Initialize an instance of the data manager:
```shell
mgr = DSGTMgr(sgt_database_folder, model3D_folder, point_cloud_file)
```

* Acquire 3D Greens function in RTZ
```shell
greens = mgr.get_greens_function(station, origin)
```

* Acquire synthetic waveform in RTZ
```shell
stream = mgr.get_waveform(station, origin, mt_RTP)
```

## Examples of using the DGF database
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
