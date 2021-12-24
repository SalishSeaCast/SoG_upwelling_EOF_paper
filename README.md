# Moore-Maley and Allen, 2021, Ocean Sci. analysis code

***

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

This repository contains a suite of interactive Jupyter notebooks and Python scripts that will allow the user to reproduce the results and figures presented in:

>B. Moore-Maley and S. E. Allen: Wind-driven upwelling and surface nutrient delivery in a semi-enclosed coastal sea, Ocean Sci., 2021.

[Conda](https://docs.conda.io/en/latest/miniconda.html) users can build the Python virtual environment necessary to run these scripts with the included `environment.yaml` file. From the Conda base environment, clone the repository to a working directory and build the environment:

```
$ git clone https://github.com/SalishSeaCast/SoG_upwelling_EOF_paper.git
$ cd SoG_upwelling_EOF_paper
$ conda update -n base conda
$ conda env create -f environment.yaml
$ source activate SoG_upwelling_EOF_paper
```

The SalishSeaCast and HRDPS results will need to be aggregated before running the analysis notebooks. The `scripts/aggregate_results.py` module is included for this task:

```
$ cd scripts
$ python3 aggregate_results.py /path/to/files
```

Finally, any of the notebooks can by run by starting a Jupyter session and navigating to the `notebooks` directory:

```
$ jupyter lab
```

## Licenses

***

The SalishSeaCast analysis and documentation are copyright 2013-2021 by the [Salish Sea MEOPAR Project Contributors](https://github.com/SalishSeaCast/docs/blob/master/CONTRIBUTORS.rst) and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.