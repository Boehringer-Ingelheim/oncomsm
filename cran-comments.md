## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* checking installed package size .. NOTE. The docs are 1.4 Mb due to some 
plots being included. The package uses rstantools to pre-compiled the stan
model used for inference this leads to a large size of the installed lib
folder.
* abort, exit, printf are found in the shared library. All of these are coming
from the rstantools dependency.
