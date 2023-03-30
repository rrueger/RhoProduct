RhoProduct
----------

A sage module implementing a slightly more generalised EtaProduct with slower
performance and different features.

Install for usage with `python` with

    pip install https://github.com/rrueger/RhoProduct/raw/main/dist/sage_rhoproduct-0.9.0-py3-none-any.whl

Sage uses its own package hierarchy to match the version of `python` it is
shipped with. To install for usage with `sage` call

    sage -pip install https://github.com/rrueger/RhoProduct/raw/main/dist/sage_rhoproduct-0.9.0-py3-none-any.whl

*Note* Sage is not officially available as a python module on PyPi, therefore it
is not listed as a dependency in this module and will **not** be installed
alongside `RhoProduct`. Sage must be installed on the system independently,
before `RhoProduct` can be installed with `sage -pip` as above.
