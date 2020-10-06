This program aims to create an NxM text image of a pendant drop profile.
Necessary parameters will be passed in the command line.
Later, I hope to allow rising drop profiles to be generated.
Aim is to make perfect drop profiles to test analysis software.
E.g. shrink droplet and observe measured IFT.
To make the drop image, I'll first set up a numerical integration to get r(z).
Then, starting several R_0 above the apex, I'll fit the needle diameter to the profile: giving z_needle.
Together, I'll then be able to map the total z(r) function (which now includes a needle) onto the pixel grid.
