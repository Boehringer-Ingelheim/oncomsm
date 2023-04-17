# oncomsm 0.1.4

* improved numerical stability (thanks to Andrew Johnson)


# oncomsm 0.1.3

* Added a `NEWS.md` file to track changes to the package.
* Made state name and censoring indicator configurable in `creat_srpmodel()`
* Reverted back to non-static C++ code
  (use rstantools to translate at build-time)
* Added `check_data()` function and using it to check data inputs where possible 
