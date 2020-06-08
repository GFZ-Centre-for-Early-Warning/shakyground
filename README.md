[![Build Status](https://travis-ci.com/gfzriesgos/shakyground.svg?branch=master)](https://travis-ci.com/gfzriesgos/shakyground)

# Shakyground

Openquake code is added here for convenience should be replaced in production
by an actual openquake installation on server

## Setup

We use the global vs30 grid created from the repository
[shakyground-grid-file](https://github.com/gfzriesgos/shakyground-grid-file).\
To be able to run shakyground locally put the resulting file as
`global_vs30.grd` within the shakyground directory (on same level as
`service.py`):

```shell
mv /path/to/valparaiso_updated_global_vs30.grd global_vs30.grd
```

You can look up some documentation about the used grid file for share wave
velocity on the usgs website:

[https://earthquake.usgs.gov/data/vs30/](https://earthquake.usgs.gov/data/vs30/)

Add to PYTHONPATH (had trouble using env.sh from openquake --> gmt was using
system gdal but openquake sqlite3 and complained):

```shell
export PYTHONPATH=$PYTHONPATH:/path/to/openquake/lib/python3.5/site-packages/
```

Required python libs:
* mock
* psutil

currently requires gmt --> will be removed

service can be tested with:

```
bash ./run_and_validate.sh
```

## Docker

To run shakyground within a Docker container you can do the following:

```shell
docker pull gfzriesgos/shakyground
docker run -it gfzriesgos/shakyground /bin/bash
# now you are in the container environment and can run shakyground by executing:
python3 service.py quakeml_test.xml
```
