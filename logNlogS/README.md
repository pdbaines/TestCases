
## Analysis of astrophysical source populations

This code implements the model described in Wong et al (2013)
[see <http://arxiv.org/abs/1305.0979>]. The code to run
the analysis in `example.R` is contained in a package
`LNSSimple`, available from <https://github.com/pdbaines/LNSSimple>.
The package can either be downloaded and installed from source,
or else installed directly from GitHub using the `devtools` package
i.e., 

    require(devtools)
    install_github(username="pdbaines",repo="LNSSimple")

The code in `example.R` can be toggled between pure R (`useC=F`)
and C code for computationally intensive portions (`useC=T`). 

Note that the `LNSSimple` also contains multi-core versions
(and hence has lots of dependencies). These could potentially
be removed upon request.

Note that (depending on the settings) this code can take a *long* time to run! 

