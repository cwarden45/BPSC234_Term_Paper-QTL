`test_mpMap2-8parentRF-estimateMap.R` - learn about functions and data structure by working backwards to generate what is needed for `estimateMap()` function (using simulated data).

## Installing *R/mpMap2* 

Because *R/mpMap2*  is no longer available in *CRAN*, I found the lastest version [in the archive](https://cran.r-project.org/src/contrib/Archive/mpMap2/) (*v1.0.4*, from 9/11/2020) and I used the latest version of R that was available at that time (**R v4.0.2**).

I found that I needed to install the following dependencies using **Packages --> Install package(s) from local files ...**:

```
ggplot2_3.3.2.tar.gz
gtable_0.3.0.tar.gz
lifecycle_1.0.0.tar.gz #I needed to use this over lifecycle_0.2.0.tar.gz
MatrixModels_0.4-1.tar.gz
pbkrtest_0.4-8.6.tar.gz
pillar_1.6.2.tar.gz #I needed to use this over pillar_1.4.6.tar.gz
pryr_0.1.4.tar.gz
stringr_1.4.0.tar.gz
```

I could use `install.packages()` for other dependencies.

I then installed *R/mpMap2* using the **mpMap2_1.0.4.tar.gz** file using *Packages --> Install package(s) from local files ...*.

## Testing *R/mpMap2* `est_map()` Function (with simulated data)

Code executed using `test_mpMap2-simulated_8parentRF-est_map.R`, with output in *testPed5_mpMap2_estimateMap-kosambi-TEST.csv*.

Because a genetic map was required for simulation, the ground truth was also saved in *testPed5_qtl_sim.map-kosambi-TRUTH.csv*.

## QTL Analysis (*outside* of *R/mpMap2*?)

I noticed that there was a [as.mpInterval](https://rdrr.io/cran/mpMap2/man/as.mpInterval.html) function

Without testing functions in detail, I could test installation of [MPWGAIM](https://github.com/KlaraVerbyla/mpwgaim).  The last commit listed on GitHub for *mpwgaim* was from 3/29/2017.

First, I attempted the following:

```
library("devtools")
devtools::install_github("KlaraVerbyla/mpwgaim")
```

In order to successfully install `devtools` for the earlier version of R, I found the earlier version (**devtools_2.3.2.tar.gz**) and I installed the following dependencies manually:

```
gh_1.4.1.tar.gz
git2r_0.27.1.tar.gz
httr2_0.2.2.tar.gz
openssl_1.4.3.tar.gz
usethis_1.6.3.tar.gz
```

So, next, I found that I needed to find earlier versions of the following packages:

```
mpMap_1.14.tar.gz
vegan_2.5-6.tar.gz
```

**There is also a required `asreml` package.**  However, I believe that a license is needed for that dependency, such as [here](https://cran.r-project.org/web/packages/biometryassist/vignettes/installing-asreml-r.html).

If I had such a license, then the `devtools` installation may work and/or may at least help to reduce the number of dependencies to install manually?

Alternatively, the following may work, if you were able to install all of the dependencies:

 - Download .zip file for `mpwgaim`
 - Uncompress .zip file
 - Remove "master" from folder name.
 - Compress .zip file as .tar file
 - Further compress .tar file to be .tar.gz file (for manual local installation of the package)
