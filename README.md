# Moore-Maley and Allen, 2022, Ocean Sci. analysis code

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

This repository contains a suite of interactive Jupyter notebooks and Python scripts that will allow the user to reproduce the results and figures presented in:

>B. Moore-Maley and S. E. Allen: Wind-driven upwelling and surface nutrient delivery in a semi-enclosed coastal sea, Ocean Sci., 2022.

A companion dataset has also been archived at the [Canadian Federated Research Data Repository (FRDR)](https://www.frdr-dfdr.ca/repo/):

>Moore-Maley, B., S. E. Allen. SalishSeaCast hourly surface along-axis wind velocity, temperature and nitrate summary 2015-2019, Federated Research Data Repository, https://doi.org/10.20383/102.0546, 2022.

This dataset needs to be obtained before running the notebooks in this code repository. Alternatively, the file can be generated using the provided scripts, however this process requires about a day to complete.

[Conda](https://docs.conda.io/en/latest/miniconda.html) users can build the Python virtual environment necessary to run these scripts with the included `environment.yaml` file. From the Conda base environment, clone the repository to a working directory and build the environment:

```
$ git clone https://github.com/SalishSeaCast/SoG_upwelling_EOF_paper.git
$ cd SoG_upwelling_EOF_paper
$ conda update -n base conda
$ conda env create -f environment.yaml
$ source activate SoG_upwelling_EOF_paper
```

Next, add the `scripts` directory to your `PYTHONPATH` environment variable to access the `tools` module:

```
$ export PYTHONPATH=$PYTHONPATH:/path/to/SoG_upwelling_EOF_paper/scripts
```

If you choose to generate the aggregated results file yourself, the `scripts/aggregate_results.py` module is included for this task:

```
$ cd scripts
$ python3 aggregate_results.py /path/to/files/
```

Either way, once you have obtained the aggregated results file, the PCA results can then be generated:

```
$ cd scripts
$ python3 PCA.py /path/to/files/
```

Finally, any of the notebooks can by run by starting a Jupyter session and navigating to the `notebooks` directory:

```
$ jupyter lab
```

## Licenses

The SalishSeaCast analysis and documentation are copyright 2013-2022 by the [Salish Sea MEOPAR Project Contributors](https://github.com/SalishSeaCast/docs/blob/master/CONTRIBUTORS.rst) and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.