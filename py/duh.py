"""
This file is part of the Duh project.
Copyright 2014 David W. Hogg (NYU) and Dustin Lang (CMU).

## bugs
* not yet written
* does nothing
"""

import numpy as np
import astropy.io.fits as pf

def synthesize_image(image_to_synthesize, basis_image):
    """
    ## input:
    - `image_to_synthesize` - new image in which we are looking for anomalies
    - `basis_image` - old image that represents "truth" or whatever

    ## output:
    - `synthetic_image` - best-fit image created from `basis_image` to approximate `image_to_synthesize`

    ## comments:
    - Images must be input with the same size and shape and very close to perfectly aligned
      (this is really a "bug" not a "comment")!

    ## bugs:
    - Not yet written.
    - Ought to have user-settable order and options related to model complexity.
    """
    return None

if __name__ == "__main__":
    print "Hello World"
