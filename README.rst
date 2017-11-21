ssc
===

The `ssc` (Spin System Creator) package provides facilitates for grouping peaks that belong to the
same spin system into spin system clusters within a single peak list.

Requirements
~~~~~~~~~~~~

   * Python 3.4+
   * gcc 5.1+

Dependencies
~~~~~~~~~~~~

   * scipy_
   * pandas_
   * bokeh_

.. _scipy: https://www.scipy.org/
.. _pandas: http://pandas.pydata.org/
.. _bokeh: http://bokeh.pydata.org/en/latest/


To test `ssc` package
~~~~~~~~~~~~~~~~~~~~~

   * Extract from tarball and cd into directory:

   .. code:: bash

      $ tar -xvf ssc-0.1.0.tar.gx
      $ cd ssc-0.1.0

   * Install dependencies:

   .. code:: bash

      $ python3 -m pip install -r requirements.txt

Instructions to build docker `ssc` container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the root of this package is a `Dockerfile` containing instructions for docker.
Use it to build docker image!

Install docker
--------------

* Follow instructions to install docker for your system: https://docs.docker.com/engine/installation
* Enable docker service.


Clone or download ssc repository
--------------------------------

* Go to the following address and clone or download ssc repo: https://gitlab.cesb.uky.edu/smelandr/ProteinNMR.ssc


Build docker image using Dockerfile
-----------------------------------

* Go to project root and execute the following command:

.. code:: bash

   # docker build -t ssc .

This will create docker image with tag name `ssc` using `Dockerfile` instructions with
`ssc` code and necessary dependencies.


Running `ssc` container with peaklist files outside container
-------------------------------------------------------------

* To display `ssc` package help message run the following command:

.. code:: bash

   # docker run -t ssc --help

* To run `ssc` docker image with external (outside of container) peaklist file
  you need to mount it as an external volume so `ssc` package inside container can
  see it. Then you pass all other required arguments to perform spin system grouping.

.. code:: bash

   # docker run -v /absolute/path/to/jr19_HNcoCACB.txt:/ssc/jr19_HNcoCACB.txt -t ssc group --plpath=/ssc/jr19_HNcoCACB.txt --plformat=sparky --stype=HNcoCACB --dims=H,N,CA/CB --rdims=H,N --view

**Note**: make sure that you specify absolute path on the host machine and path where it will be
visible on the container side after `:`. Also make sure that path were you mount volume on the
container side matches `--plpath` command line argument.
