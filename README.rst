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