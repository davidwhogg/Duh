"""
This file is part of the Duh project.
Copyright 2014 David W. Hogg (NYU) and Dustin Lang (CMU).

## bugs
* not yet written
* does nothing
"""

import numpy as np
import astropy.io.fits as pf

def populate_A_matrix(basis_image):
    """
    ## bugs:
    - Needs more information in this comment header.
    - Not yet written.
    - Ought to have user-settable order and options related to model complexity.
    - There must be a MUCH faster way to write this.
    """
    ny, nx = basis_image.shape
    npix = nx * ny
    A = np.zeros((10, npix))
    A[0] = np.ones(npix)
    A[1] = basis_image.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:, :-1] = basis_image[:, 1:]
    A[2] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:, 1:] = basis_image[:, :-1]
    A[3] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:-1] = basis_image[1:]
    A[4] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[1:] = basis_image[:-1]
    A[5] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:-1, :-1] = basis_image[1:, 1:]
    A[6] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:-1, 1:] = basis_image[1:, :-1]
    A[7] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[1:, :-1] = basis_image[:-1, 1:]
    A[8] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[1:, 1:] = basis_image[:-1, :-1]
    A[9] = foo.reshape(npix)
    return A

def synthesize_image(image_to_synthesize, weight_image, basis_image):
    """
    ## input:
    - `image_to_synthesize` - new image in which we are looking for anomalies
    - `weight_image` - image of inverse variance values or 1 (to use in fit) and 0 (to not use) values.
    - `basis_image` - old image that represents "truth" or whatever

    ## output:
    - `synthetic_image` - best-fit image created from `basis_image` to approximate `image_to_synthesize`

    ## comments:
    - Images must be input with the same size and shape and very close to perfectly aligned
      (this is really a "bug" not a "comment")!

    ## bugs:
    - Not yet written.
    - Ought to have user-settable order and options related to model complexity.
    - Can't take multiple `basis_image`s.  That would be cool!
    """
    assert weight_image.shape == image_to_synthesize.shape
    assert basis_image.shape == image_to_synthesize.shape
    ny, nx = image_to_synthesize.shape
    yg = np.arange(ny) + 0.5 - 0.5 * ny
    xg = np.arange(nx) + 0.5 - 0.5 * nx
    yg, xg = np.meshgrid(xg, yg)
    npix = ny * nx
    A = np.vstack((populate_A_matrix(basis_image),
                   populate_A_matrix(xg * basis_image),
                   populate_A_matrix(yg * basis_image)))
    W = weight_image.reshape(npix)
    ATA = np.dot(W[None, :] * A, A.T)
    ATb = np.dot(W[None, :] * A, image_to_synthesize.reshape(npix))
    pars = np.linalg.solve(ATA, ATb)
    return np.dot(A.T, pars).reshape((ny, nx))

if __name__ == "__main__":
    new = np.ones((3, 5))
    weight = np.ones((3, 5))
    old = np.random.uniform(size=(3, 5))
    print synthesize_image(new, weight, old)
