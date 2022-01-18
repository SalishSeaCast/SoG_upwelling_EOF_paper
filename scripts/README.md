# Moore-Maley and Allen, 2022, Ocean Sci. scripts

This directory contains the Python scripts that generate the necessary datasets to perform the analysis presented in

>B. Moore-Maley and S. E. Allen: Wind-driven upwelling and surface nutrient delivery in a semi-enclosed coastal sea, Ocean Sci., 2022.

These scripts can be run by cloning this repository and setting up a conda environment, then executing each script at the command line with a results path provided. For example, if you choose to generate the aggregated results file rather than downloading from the data archive on the FRDR repository, the following command will aggregate the results files from the [SalishSeaCast erddap server](https://salishsea.eos.ubc.ca/erddap/)

```
$ python3 aggregate_results.py /path/to/files/
```

The following command will perform the PCA analysis.

```
$ python3 PCA.py /path/to/files/
```

## Licenses

The SalishSeaCast analysis and documentation are copyright 2013-2021 by the [Salish Sea MEOPAR Project Contributors](https://github.com/SalishSeaCast/docs/blob/master/CONTRIBUTORS.rst) and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.