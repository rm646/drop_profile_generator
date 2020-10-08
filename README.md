## Introduction

This program creates an NxM pixels image of a pendant drop profile by solving
the Young-LaPlace equation numerically. The point is to make an idealised drop
profile, but mapped onto a pixel grid, so that one can generate fake
pendant drop profiles to test the analysis algorithms which extract surface
tension from those profiles. This program currently only deals with hanging drop
profiles, and does not allow rising drop profiles.

[Berry et al.'s  paper](https://doi.org/10.1016/j.jcis.2015.05.012) provides a
good background on pendant drop measurements and the equations solved here,
and also links to the [OpenDrop project](https://github.com/FrostadResearch/Pendant-drop-tensiometer-v2/tree/master/OpenDrop%20Software%20Package) for measuring surface tension from pendant drop
profiles.

## Usage

This is for Python 3:

```
python dpg_main.py [-h] [-nw NEEDLEWIDTH] [-ift INTERFACIALTENSION]
    [-dd DENSITYDIFF] [-R0 DROPDIM] [-N GRIDWIDTH] [-M GRIDHEIGHT] [-s SCALE]
```

The output is a png, with all parameters and the run date included in the
filename.

Also included is a simple bash script I used to scan a range of drop volumes.
