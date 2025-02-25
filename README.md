# orbitalthebox python package

## Short description
A package for using orbital parameters to carry out various irradiance calculations for Earth. Beta version. A release version is coming, which will also have a jupyter notebook tutorial.

## Install
You can install into your python environment using the `pip` terminal command (you may need to [install](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) `git` first):

`pip install git+https://github.com/bryanlougheed/orbitalthebox.git`

If the install has been successul then you should be able to import it in python in the usual way, e.g.:

`import orbitalthebox`

## License
The GPL license is permissive but please be reasonable.

## List of functions
Each function is has help documentation that you can access using the standard python `?` command. More functions will follow as I bug test them in python.

`getlaskar2004()`
`getlaskar2010()`
`solvekeplerE()`
`sollon2time()`
`time2sollon()`
`geographiclat()`
`dailymeanwm2()`
`intradaywm2()`
`thresholdjm2()`
`areaquad()`
