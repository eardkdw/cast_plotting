Routines for plotting ECMWF Grib data for the [CAST][1] Project.

The cast_plot.py file can be run on the command line as follows:
```
$ python cast_plot.py <GRIB-FILE>
```

The script will create a series of plots for each of the variables
defined inside the script as CAST_STANDARD_NAMES and CAST_LONG_NAMES.

## Requirements

This code is designed to run on a system using the [JASMIN Analysis Platform][2].  In particular it requires
 1. IRIS
 2. Cartopy
 3. ECMWF GribAPI
 4. All dependencies of the above

## TODO

 1. Allow configuration of which level from each variable to plot
 2. Better command-line interface
 3. Wind-vector plot of wind fields
 4. Better annotation of plots


[1]: http://www.faam.ac.uk/index.php/current-future-campaigns/384-cast-2014-co-ordinated-airborne-studies-in-the-tropics
[2]: http://proj.badc.rl.ac.uk/cedaservices/wiki/JASMIN/AnalysisPlatform
