**paRseMSF** is an R package to parse the information in MSF files resulting from a Thermo Proteome Discoverer analysis workflow. The MSF files use a SQLite, a lightweight, self-contained, relational database engine to store the results and the workflow parameters. The paRseMSF package for R uses the database interface drivers available in R ([RSQLite](http://cran.r-project.org/web/packages/RSQLite/index.html), a package that allows R to access SQLite databases) and [sqldf](https://code.google.com/p/sqldf/) package for R to access and process the tables stored in MSF files.

**To Install:**
  * Download the [zip file](http://parsemsf.googlecode.com/files/paRseMSF_1.0.zip)
  * In R, go to: Packages -> Install package(s) from local zip files.
  * Since paRseMSF is not yet up on CRAN, you may need to install the dependencies, i.e. the two R packages: RSQLite and sqldf.

---

  * Polpitiya et al, ASMS 2012, Vancouver, Canada, paRseMSF: R package for parsing Thermo MSF files. ([Poster](http://parsemsf.googlecode.com/files/paRseMSF_apolpitiya_2012_asms.pdf), [Abstract](http://parsemsf.googlecode.com/files/paRseMSF_apolpitiya_abstract.pdf))


---

![http://parsemsf.googlecode.com/files/screenshot2.png](http://parsemsf.googlecode.com/files/screenshot2.png)