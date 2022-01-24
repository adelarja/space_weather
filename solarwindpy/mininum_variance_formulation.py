# -*- coding: utf-8 -*-
# This file is part of the
# 	Solarwindpy Project(https://github.com/adelarja/space_weather).

# Copyright (c) 2021, Adriana Gulisano, Adel Arja, Ricardo Pafundigit,
# Violeta Bazzano.
# All rights reserved.

# License: BSD 3-Clause License
# 	Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE


from enum import Enum

import numpy as np

import scipy.linalg as la

from solarwindpy.wind import MagneticField


class KindMv(Enum):
    X_VERSOR = 1
    Z_VERSOR = 2


KIND_MV = KindMv.Z_VERSOR


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


def correct_bz(wind_bz, wind_z_versor, bx, by, bz):
    middle_element_index = len(wind_bz) // 2
    if (
        wind_bz[middle_element_index - 1] < 0
        and wind_bz[middle_element_index] < 0
        and wind_bz[middle_element_index + 1] < 0
    ):
        wind_z_versor = -wind_z_versor
        wind_bz = (
            bx * wind_z_versor[0]
            + by * wind_z_versor[1]
            + bz * wind_z_versor[2]
        )
    return wind_bz, wind_z_versor


def correct_bx(wind_bx, wind_x_versor, bx, by, bz):
    mod_alfa = np.arccos(wind_x_versor[0]) * 180 / np.pi
    if mod_alfa > 85:
        wind_x_versor = -wind_x_versor
        wind_bx = (
            bx * wind_x_versor[0]
            + by * wind_x_versor[1]
            + bz * wind_x_versor[2]
        )
        mod_alfa = np.arccos(wind_x_versor[0]) * 180 / np.pi
        if mod_alfa >= 85:
            raise ValueError(
                "ERROR:change sgn x_versor didn´t do x_cloud=r(out_bound)"
            )
    return wind_bx, wind_x_versor


def correct_by(
    wind_by, wind_y_versor, wind_x_versor, wind_z_versor, bx, by, bz
):
    transposed_matrix_determinant = get_transposed_matrix_determinant(
        wind_x_versor, wind_y_versor, wind_z_versor
    )
    if abs(transposed_matrix_determinant + 1) < 1e-10:  # MM_provi_det == -1
        wind_y_versor = -wind_y_versor
        wind_by = (
            bx * wind_y_versor[0]
            + by * wind_y_versor[1]
            + bz * wind_y_versor[2]
        )
        transposed_matrix_determinant = get_transposed_matrix_determinant(
            wind_x_versor, wind_y_versor, wind_z_versor
        )
        if abs(transposed_matrix_determinant - 1) > 1e-10:
            raise ValueError(
                "ERROR: The change to righ hand-handed didnit work"
            )
    transposed_matrix = get_transposed_matrix(
        wind_x_versor, wind_y_versor, wind_z_versor
    )

    return wind_by, wind_y_versor, transposed_matrix


def get_transposed_matrix(wind_x_versor, wind_y_versor, wind_z_versor):
    transposed_matrix = np.ones((3, 3))
    transposed_matrix[0, :] = wind_x_versor.reshape(1, 3)
    transposed_matrix[1, :] = wind_y_versor.reshape(1, 3)
    transposed_matrix[2, :] = wind_z_versor.reshape(1, 3)
    return transposed_matrix


def get_transposed_matrix_determinant(
    wind_x_versor, wind_y_versor, wind_z_versor
):
    transposed_matrix = get_transposed_matrix(
        wind_x_versor, wind_y_versor, wind_z_versor
    )
    transposed_matrix_determinant = np.linalg.det(transposed_matrix)
    return transposed_matrix_determinant


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


def rotation(x_gse, y_gse, z_gse, theta_deg, phi_deg):
    # from deg to radians
    theta = theta_deg * np.pi / 180.0
    phi = phi_deg * np.pi / 180.0
    gamma = calc_gamma(phi, theta)
    # building the R matrix
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
    RT = np.transpose(R)
    # applying the rotation by hand
    x_mc = RT[0][0] * x_gse + RT[0][1] * y_gse + RT[0][2] * z_gse
    y_mc = RT[1][0] * x_gse + RT[1][1] * y_gse + RT[1][2] * z_gse
    z_mc = RT[2][0] * x_gse + RT[2][1] * y_gse + RT[2][2] * z_gse

    return x_mc, y_mc, z_mc


def get_main_versor(kind_mv, x_versor_nube, z_versor_nube):
    if kind_mv == KindMv.Z_VERSOR:
        main_versor = z_versor_nube
    elif kind_mv == KindMv.X_VERSOR:
        main_versor = x_versor_nube
    else:
        raise ValueError("kind_mv needs to be 1 or 2")
    return main_versor


class RotatedWind:
    def __init__(self, bgse0, bgse1, bgse2):
        self.bgse0 = bgse0
        self.bgse1 = bgse1
        self.bgse2 = bgse2

    @classmethod
    def get_rotated_wind(cls, wind: list[MagneticField]):
        bx = np.array([magnetic_field.bgse0 for magnetic_field in wind])
        by = np.array([magnetic_field.bgse1 for magnetic_field in wind])
        bz = np.array([magnetic_field.bgse2 for magnetic_field in wind])

        minimum_variance_matrix = calculate_minimum_variance(bx, by, bz)
        eigvals, eigvecs = la.eig(minimum_variance_matrix)
        eigvals = eigvals.real
        sorted_index = np.argsort(eigvals)

        wind_x_versor = eigvecs[:, sorted_index[0]].reshape(3, 1)
        wind_z_versor = eigvecs[:, sorted_index[1]].reshape(3, 1)
        wind_y_versor = eigvecs[:, sorted_index[2]].reshape(3, 1)

        wind_bx = (
            bx * wind_x_versor[0]
            + by * wind_x_versor[1]
            + bz * wind_x_versor[2]
        )
        wind_by = (
            bx * wind_y_versor[0]
            + by * wind_y_versor[1]
            + bz * wind_y_versor[2]
        )
        wind_bz = (
            bx * wind_z_versor[0]
            + by * wind_z_versor[1]
            + bz * wind_z_versor[2]
        )

        wind_bx, wind_x_versor = correct_bx(wind_bx, wind_x_versor, bx, by, bz)
        wind_bz, wind_z_versor = correct_bz(wind_bz, wind_z_versor, bx, by, bz)
        wind_by, wind_y_versor, _ = correct_by(
            wind_by, wind_y_versor, wind_x_versor, wind_z_versor, bx, by, bz
        )

        return cls(wind_bx, wind_by, wind_bz)

    @classmethod
    def get_rotated_wind_using_transposed_matrix(
        cls, wind: list[MagneticField]
    ):
        bx = np.array([magnetic_field.bgse0 for magnetic_field in wind])
        by = np.array([magnetic_field.bgse1 for magnetic_field in wind])
        bz = np.array([magnetic_field.bgse2 for magnetic_field in wind])

        minimum_variance_matrix = calculate_minimum_variance(bx, by, bz)
        eigvals, eigvecs = la.eig(minimum_variance_matrix)
        eigvals = eigvals.real
        sorted_index = np.argsort(eigvals)

        wind_x_versor = eigvecs[:, sorted_index[0]].reshape(3, 1)
        wind_z_versor = eigvecs[:, sorted_index[1]].reshape(3, 1)
        wind_y_versor = eigvecs[:, sorted_index[2]].reshape(3, 1)

        wind_bx = (
            bx * wind_x_versor[0]
            + by * wind_x_versor[1]
            + bz * wind_x_versor[2]
        )
        wind_by = (
            bx * wind_y_versor[0]
            + by * wind_y_versor[1]
            + bz * wind_y_versor[2]
        )
        wind_bz = (
            bx * wind_z_versor[0]
            + by * wind_z_versor[1]
            + bz * wind_z_versor[2]
        )

        wind_bx, wind_x_versor = correct_bx(wind_bx, wind_x_versor, bx, by, bz)
        wind_bz, wind_z_versor = correct_bz(wind_bz, wind_z_versor, bx, by, bz)
        wind_by, wind_y_versor, transposed_matrix = correct_by(
            wind_by, wind_y_versor, wind_x_versor, wind_z_versor, bx, by, bz
        )

        rotated_bx = (
            transposed_matrix[0][0] * bx
            + transposed_matrix[0][1] * by
            + transposed_matrix[0][2] * bz
        )
        rotated_by = (
            transposed_matrix[1][0] * bx
            + transposed_matrix[1][1] * by
            + transposed_matrix[1][2] * bz
        )
        rotated_bz = (
            transposed_matrix[2][0] * bx
            + transposed_matrix[2][1] * by
            + transposed_matrix[2][2] * bz
        )

        return cls(rotated_bx, rotated_by, rotated_bz)

    @classmethod
    def get_rotated_wind_using_angles(cls, wind: list[MagneticField]):
        bx = np.array([magnetic_field.bgse0 for magnetic_field in wind])
        by = np.array([magnetic_field.bgse1 for magnetic_field in wind])
        bz = np.array([magnetic_field.bgse2 for magnetic_field in wind])

        minimum_variance_matrix = calculate_minimum_variance(bx, by, bz)
        eigvals, eigvecs = la.eig(minimum_variance_matrix)
        eigvals = eigvals.real
        sorted_index = np.argsort(eigvals)

        wind_x_versor = eigvecs[:, sorted_index[0]].reshape(3, 1)
        wind_z_versor = eigvecs[:, sorted_index[1]].reshape(3, 1)
        wind_y_versor = eigvecs[:, sorted_index[2]].reshape(3, 1)

        wind_bx = (
            bx * wind_x_versor[0]
            + by * wind_x_versor[1]
            + bz * wind_x_versor[2]
        )
        wind_by = (
            bx * wind_y_versor[0]
            + by * wind_y_versor[1]
            + bz * wind_y_versor[2]
        )
        wind_bz = (
            bx * wind_z_versor[0]
            + by * wind_z_versor[1]
            + bz * wind_z_versor[2]
        )

        wind_bx, wind_x_versor = correct_bx(wind_bx, wind_x_versor, bx, by, bz)
        wind_bz, wind_z_versor = correct_bz(wind_bz, wind_z_versor, bx, by, bz)
        _, _, transposed_matrix = correct_by(
            wind_by, wind_y_versor, wind_x_versor, wind_z_versor, bx, by, bz
        )

        tita_deg = np.arcsin(transposed_matrix[2][2]) * 180 / np.pi
        tita = np.arcsin(transposed_matrix[2][2])
        fi = np.arctan2(
            transposed_matrix[1][2] / np.cos(tita),
            transposed_matrix[0][2] / np.cos(tita),
        )
        fi_deg = fi * 180 / np.pi
        if fi_deg < 0:
            fi_deg = 360 + fi_deg
        else:
            fi_deg = fi_deg

        bx_n, by_n, bz_n = rotation(bx, by, bz, tita_deg, fi_deg)

        return cls(bx_n, by_n, bz_n)
