# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 10:32:58 2021

@author: aguli
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

#############################################################################
# Here I define some useful funtions to calculate the third grade of liberty
# for the Euler rotation given phi and theta angles, but chosen this third
# angle  gamma in such a way that
#  the direction of the x coordinate in the Magnetic Cloud frame fulfill
#  x_MC dot x_gse >0
#######################################################################


def calc_gamma(phi, theta):
    """
    Here we test if x_MC dot x_gse >0
    """
    epsilon = 10e-15
    gamma = np.arctan(-np.tan(phi) / np.sin(theta + epsilon))
    xmc_dot_xgse = np.cos(gamma) * np.sin(theta) * np.cos(phi) - np.sin(
        gamma
    ) * np.sin(phi)
    if xmc_dot_xgse < 0:
        gamma = np.arctan(-np.tan(phi) / np.sin(theta + epsilon)) + np.pi
        xmc_dot_xgse = np.cos(gamma) * np.sin(theta) * np.cos(phi) - np.sin(
            gamma
        ) * np.sin(phi)
        if xmc_dot_xgse < 0:
            raise ValueError("Error in gamma cuadrant")
    return gamma


######################################################################
# This funcion rotates any vector in gse coordinates according to theta,
# and phi rotation angles
#############################################################################


def rotation(x_gse, y_gse, z_gse, theta_deg, phi_deg):
    # from deg to radians
    theta = theta_deg * np.pi / 180.0
    phi = phi_deg * np.pi / 180.0
    gamma = calc_gamma(phi, theta)
    # building  the R matrix
    R = np.zeros((3, 3))
    R[0][0] = np.cos(gamma) * np.sin(theta) * np.cos(phi) - np.sin(
        gamma
    ) * np.sin(phi)
    R[0][1] = np.cos(gamma) * np.sin(theta) * np.sin(phi) + np.sin(
        gamma
    ) * np.cos(phi)
    R[0][2] = -np.cos(gamma) * np.cos(theta)
    R[1][0] = -np.sin(gamma) * np.sin(theta) * np.cos(phi) - np.cos(
        gamma
    ) * np.sin(phi)
    R[1][1] = -np.sin(gamma) * np.sin(theta) * np.sin(phi) + np.cos(
        gamma
    ) * np.cos(phi)
    R[1][2] = np.sin(gamma) * np.cos(theta)
    R[2][0] = np.cos(theta) * np.cos(phi)
    R[2][1] = np.cos(theta) * np.sin(phi)
    R[2][2] = np.sin(theta)
    # applying the rotation by hand
    x_mc = R[0][0] * x_gse + R[0][1] * y_gse + R[0][2] * z_gse
    y_mc = R[1][0] * x_gse + R[1][1] * y_gse + R[1][2] * z_gse
    z_mc = R[2][0] * x_gse + R[2][1] * y_gse + R[2][2] * z_gse

    return x_mc, y_mc, z_mc


##############################################################################
# Here I loaded the components of B in gse coordinates for a toy_cloud to test
# everything is working well
##############################################################################

bx = np.load("bx.npy")  # load bx in gse coordinates to test
by = np.load("by.npy")  # load by in gse coordinates to test
bz = np.load("bz.npy")  # load bz in gse coordinates to test


kind_mv = 1  # kind_mv=1 if its used to orient Magnetic Clouds

# Here I compute de Minimum Variance Matrix for the B components
# in GSE coordinates
####################################################
# MV matrix (M)
####################################################
# sigma^2_n=M*n*n
####################################################


def calculate_minimum_variance(bx, by, bz):
    minimum_variance_matrix = np.ones((3, 3))
    minimum_variance_matrix[0][0] = np.mean(bx ** 2) - np.mean(bx) * np.mean(
        bx
    )
    minimum_variance_matrix[1][1] = np.mean(by ** 2) - np.mean(by) * np.mean(
        by
    )
    minimum_variance_matrix[2][2] = np.mean(bz ** 2) - np.mean(bz) * np.mean(
        bz
    )
    minimum_variance_matrix[0][1] = np.mean(bx * by) - np.mean(bx) * np.mean(
        by
    )
    minimum_variance_matrix[0][2] = np.mean(bx * bz) - np.mean(bx) * np.mean(
        bz
    )
    minimum_variance_matrix[1][2] = np.mean(by * bz) - np.mean(by) * np.mean(
        bz
    )
    minimum_variance_matrix[1][0] = minimum_variance_matrix[0][1]
    minimum_variance_matrix[2][0] = minimum_variance_matrix[0][2]
    minimum_variance_matrix[2][1] = minimum_variance_matrix[1][2]

    return minimum_variance_matrix


def get_eigvals_and_eigvects(minimum_variance_matrix):
    eigvals, eigvecs = la.eig(minimum_variance_matrix)
    eigvals = eigvals.real
    sorted_index = np.argsort(eigvals)

    x_versor_nube = eigvecs[:, sorted_index[0]].reshape(3, 1)
    z_versor_nube = eigvecs[:, sorted_index[1]].reshape(3, 1)
    y_versor_nube = eigvecs[:, sorted_index[2]].reshape(3, 1)

    bx_nube = bx * x_versor_nube[0] + by * x_versor_nube[1] + bz * x_versor_nube[2]
    by_nube = bx * y_versor_nube[0] + by * y_versor_nube[1] + bz * y_versor_nube[2]
    bz_nube = bx * z_versor_nube[0] + by * z_versor_nube[1] + bz * z_versor_nube[2]

    NN = len(bz_nube) // 2
    if bz_nube[NN - 1] < 0 and bz_nube[NN] < 0 and bz_nube[NN + 1] < 0:
        z_versor_nube = -n_int
        # z_nube
        bx_nube = (
                bx * x_versor_nube[0] + by * x_versor_nube[1] + bz * x_versor_nube[2]
        )
        by_nube = (
                bx * y_versor_nube[0] + by * y_versor_nube[1] + bz * y_versor_nube[2]
        )
        bz_nube = (
                bx * z_versor_nube[0] + by * z_versor_nube[1] + bz * z_versor_nube[2]
        )

    return bx_nube, by_nube, bz_nube


def validate_cloud():
    pass


def get_rotation_angles():
    pass


minimum_variance_matrix = calculate_minimum_variance(bx, by, bz)


##############################################################
# Eigen vector and eigenvalues for the M Matrix are calculated here
# with these the problem is basically solved
#
##############################################################
eigvals, eigvecs = la.eig(minimum_variance_matrix)
eigvals = eigvals.real  # They have to be real ones

# I used these prints to test the results, so the prints can be deleted later
# print('-------------------------------------------')
# print(M @ eigvecs[:,0].reshape(3,1))
# print(eigvals[0]*eigvecs[:,0].reshape(3,1))
# print('-------------------------------------------')


#############################################################
# Here I order the eigen values from minimum to max value
##############################################################
data = eigvals
sorted_index = np.argsort(data)
# print(sorted_index) # sorted index
# print(data[sorted_index]) # values selected by the sorted index

######################################################################
# Here is the physics of the problem
# I asign the minimum eigen value and eigen vector
# to the minimum variance and the  direction where the variance is
# minimum
# the intermediate value to the intermediate variance direction
# and the maximum value to the maximum variance direction
# as expected for a Magnetic Cloud as in  [Bothmer & Schwenn, AnnGeophys
# 1998] paper
######################################################################

lambda1 = eigvals[sorted_index[0]]  # minimun
lambda2 = eigvals[sorted_index[1]]  # intermediate
lambda3 = eigvals[sorted_index[2]]  # maximum

v1 = eigvecs[:, sorted_index[0]].reshape(3, 1)
v2 = eigvecs[:, sorted_index[1]].reshape(3, 1)
v3 = eigvecs[:, sorted_index[2]].reshape(3, 1)

#####################################################################
# I renamed to the notation in physics where a direction is named n
######################################################################

n_min = v1
n_max = v3
n_int = v2

x_versor_nube = n_min
y_versor_nube = n_max
z_versor_nube = n_int

# Rotation of B (transformation to the cloud coordinates).
# coordinates  x_versor_nube, y_versor_nube, z_versor_nube.
bx_nube = bx * x_versor_nube[0] + by * x_versor_nube[1] + bz * x_versor_nube[2]
by_nube = bx * y_versor_nube[0] + by * y_versor_nube[1] + bz * y_versor_nube[2]
bz_nube = bx * z_versor_nube[0] + by * z_versor_nube[1] + bz * z_versor_nube[2]


#########################################################################
# Here I
# Check z_versor_nube sign: Bz(0)>0 otherwise I fixed that
# changing the sign z_versor_nube = -n_int
#########################################################################
NN = len(bz_nube) // 2
if bz_nube[NN - 1] < 0 and bz_nube[NN] < 0 and bz_nube[NN + 1] < 0:
    z_versor_nube = -n_int
    # z_nube
    bx_nube = (
        bx * x_versor_nube[0] + by * x_versor_nube[1] + bz * x_versor_nube[2]
    )
    by_nube = (
        bx * y_versor_nube[0] + by * y_versor_nube[1] + bz * y_versor_nube[2]
    )
    bz_nube = (
        bx * z_versor_nube[0] + by * z_versor_nube[1] + bz * z_versor_nube[2]
    )
    print("z_versor_nube es-n_int")

#######################################################
# Here I Check that  x_versor_nube sign: x_cloud=r(out_bound)
#######################################################
# The convention is x_cloud is in agreement
# with r_versor for the out-bound branch,
# This is checked computing alfa (angle between x_versor_nube and x_GSE)
# |alfa| needs to be lower than 90 deg.
# |alfa| = arc cos (x_versor_nube  x_GSE)
# If |alfa| < 90 is ok, if not sign in definition of x_versor_nube need to be
# changed
#########################################################
mod_alfa = np.arccos(x_versor_nube[0]) * 180 / np.pi
if mod_alfa < 85:
    print("Def. of x_versor_nube is Ok.")
if mod_alfa > 85:
    print("Def. of x_versor_nube is OPOSITE. INVERSE SIGN")
    x_versor_nube = -n_min
    # r
    bx_nube = (
        bx * x_versor_nube[0] + by * x_versor_nube[1] + bz * x_versor_nube[2]
    )
    by_nube = (
        bx * y_versor_nube[0] + by * y_versor_nube[1] + bz * y_versor_nube[2]
    )
    bz_nube = (
        bx * z_versor_nube[0] + by * z_versor_nube[1] + bz * z_versor_nube[2]
    )
    mod_alfa = np.acos(x_versor_nube[0]) * 180 / np.pi
    if mod_alfa >= 85:
        raise ValueError(
            "ERROR:change sgn x_versor didn´t do x_cloud=r(out_bound)"
        )
    else:
        raise ValueError("ERROR alfa GIVES near 90 degrees")
#################################################################
# Check if {x_versor_nube, y_versor_nube, z_versor_nube} is right-handed
#################################################################
MM_provi = np.ones((3, 3))
MM_provi[0, :] = x_versor_nube.reshape(1, 3)
MM_provi[1, :] = y_versor_nube.reshape(1, 3)
MM_provi[2, :] = z_versor_nube.reshape(1, 3)
MM_provi_det = np.linalg.det(MM_provi)
if abs(MM_provi_det - 1) < 1e-10:  # MM_provi_det == +1
    print("Def. d y_versor_nube Ok.")
if abs(MM_provi_det + 1) < 1e-10:  # MM_provi_det == -1
    print("Def.  y_versor_nube inverted. Sign changed.")
    y_versor_nube = -n_max
    # phi or y_mc
    bx_nube = (
        bx * x_versor_nube[0] + by * x_versor_nube[1] + bz * x_versor_nube[2]
    )
    by_nube = (
        bx * y_versor_nube[0] + by * y_versor_nube[1] + bz * y_versor_nube[2]
    )
    bz_nube = (
        bx * z_versor_nube[0] + by * z_versor_nube[1] + bz * z_versor_nube[2]
    )
    MM_provi[0, :] = x_versor_nube.reshape(1, 3)
    MM_provi[1, :] = y_versor_nube.reshape(1, 3)
    MM_provi[2, :] = z_versor_nube.reshape(1, 3)
    MM_provi_det = np.linalg.det(MM_provi)
    if abs(MM_provi_det - 1) > 1e-10:
        raise ValueError("ERROR: The change to righ hand-handed didnit work")


###############################################################33
# here I checked if it is used for Magnetic Clouds or to orient
# Shocks in the interplanetary medium
# for MC is kind_mv=1
#################################################################


if kind_mv == 1:
    main_versor = z_versor_nube
elif kind_mv == 2:
    main_versor = x_versor_nube
else:
    raise ValueError("kind_mv needs to be 1 or 2")


############################################################################
# Computation of theta, phi, and gamma from the computed rotation matrix
# theta: latitud, on ecliptic theta=0, range:[-90,+90]
# phi: angle between projection of z_n on ecpliptic, measured from x_gse
# (phi=0) to +y_gse (phi=90)
# and beyond, range:[0,+360].
#
# gamma: angle to rotate the perpendicular plane to the MC axis
#
#############################################################################
# theta:
# the scalar product between z_nube and z_gse is: n_z_nube(3)
# and thus z_versor_nube(3)=cos(theta_monio)
# theta_monio_min_var is inside [0,180]
theta_monio_min_var = np.arccos(main_versor[2]) * 180 / np.pi
# theta_monio: angle between z_gse and z_cloud, so:
theta_min_var = 90 - theta_monio_min_var
#############################################################################


mod = np.sqrt(main_versor[0] ** 2 + main_versor[1] ** 2)
p_zn_eclip_versor_x_gse = main_versor[0] / mod
p_zn_eclip_versor_y_gse = main_versor[1] / mod
cos_de_phi_min_var = p_zn_eclip_versor_x_gse
sin_de_phi_min_var = p_zn_eclip_versor_y_gse

print(
    '\phi (between "proyection Mc axis to the ecliptic plane" , "axis x_GSE"'
)
if sin_de_phi_min_var > 0 and cos_de_phi_min_var > 0:
    cuadrante = 1
    print("first cuadrant (0,90)")
    phi_min_var = np.arccos(cos_de_phi_min_var) * 180 / np.pi
elif sin_de_phi_min_var > 0 and cos_de_phi_min_var < 0:
    cuadrante = 2
    print("second Cuadrant (90,180)")
    phi_min_var = np.arccos(cos_de_phi_min_var) * 180 / np.pi
elif sin_de_phi_min_var < 0 and cos_de_phi_min_var < 0:
    cuadrante = 3
    print("Third Cuadrant (180,270)")
    phi_min_var = 360 - np.arccos(cos_de_phi_min_var) * 180 / np.pi
elif sin_de_phi_min_var < 0 and cos_de_phi_min_var > 0:
    cuadrante = 4
    print("Forth Cuadrant (270,360)")
    phi_min_var = 360 - np.arccos(cos_de_phi_min_var) * 180 / np.pi
else:
    raise ValueError("can not determine cuadrant Error")


##################################################################
# Rotate with the angles calculated earlier
##################################################################

[bx_n, by_n, bz_n] = rotation(bx, by, bz, theta_min_var, phi_min_var)

minimum_variance_matrix = MM_provi  # definitive variance Matrix

##################################################################
# here I get the Magnetic cloud components from the MV definitive
# matrix directly
# to check everything is working well
# I should be the same as the ones calculated from the rotation angles
#####################################################################
bx_mc = (
    minimum_variance_matrix[0][0] * bx
    + minimum_variance_matrix[0][1] * by
    + minimum_variance_matrix[0][2] * bz
)
by_mc = (
    minimum_variance_matrix[1][0] * bx
    + minimum_variance_matrix[1][1] * by
    + minimum_variance_matrix[1][2] * bz
)
bz_mc = (
    minimum_variance_matrix[2][0] * bx
    + minimum_variance_matrix[2][1] * by
    + minimum_variance_matrix[2][2] * bz
)

##################################################################
# Here I plot B in gse coordinates
##################################################################
plt.figure(1)
plt.plot(by)
plt.plot(bz)
plt.plot(bx)

###################################################################
# Here I plot B in Magnetic cloud frame from the MV definitive variance
# Matrix
#########################################################################
plt.figure(2)
plt.plot(by_mc)
plt.plot(bz_mc)
plt.plot(bx_mc)

#######################################################################
# Here I plot B in Magnetic Cloud frame from the angles calculated
# from the matrix eigen value and eigen problem
# I should be the same as figure (2)
##########################################################################
plt.figure(3)
plt.plot(by_nube)
plt.plot(bz_nube)
plt.plot(bx_nube)
