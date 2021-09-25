# -*- coding: utf-8 -*-
# This file is part of the
#   Instituto de Astronomico y Fisica del Espacio (IAFE, CONICET-UBA)
# Copyright (c) 2021, Departamento de Fisica, FCEN (UBA).
# v1.2, by S. Dasso, Dec 20, 2007.
#   Full Text: https://github.com/adelarja/space_weather/blob/main/LICENSE.md

# =====================TEMPORARY (REMOVE AFTERWARDS)===========================

# =============================================================================


# =============================================================================
# IMPORTS
# =============================================================================

import sys
import numpy as np
import pandas as pd
import scipy.linalg as la
import cdflib

# ============================================================================
# CLASSES
# ============================================================================
class ElementosMatematicos:
    def __init__(self):
        self.M = pd.DataFrame({})
        self.epsilon = 1e-4
        pass

    def min_variance(self):
        """
        % MV matrix (M)
        % sigma^2_n=M*n*n
        M(1,1)=mean(bx.^2)-mean(bx)^2;
        M(2,2)=mean(by.^2)-mean(by)^2;
        M(3,3)=mean(bz.^2)-mean(bz)^2;
        M(1,2)=mean(bx.*by)-mean(bx)*mean(by);
        M(1,3)=mean(bx.*bz)-mean(bx)*mean(bz);
        M(2,3)=mean(by.*bz)-mean(by)*mean(bz);
        M(2,1)=M(1,2);
        M(3,1)=M(1,3);
        M(3,2)=M(2,3);
        """
        return 1

    def ordered_eigen(self, M):
        if M.shape[0] == M.shape[1]:
            eigvals, eigvecs = la.eig(M)
            # do ordering
        else:
            raise ValueError("Matrix should be square")


"""
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigenvalues and eigenvectors
[V,D] =eig(M);
% Eigenvalues nor ordered at this stage.
lambda(1)=D(1,1);
lambda(2)=D(2,2);
lambda(3)=D(3,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of three unit vectors to get orientation of the structure:
n(1,1)=V(1,1); n(1,2)=V(2,1); n(1,3)=V(3,1);
n(2,1)=V(1,2); n(2,2)=V(2,2); n(2,3)=V(3,2);
n(3,1)=V(1,3); n(3,2)=V(2,3); n(3,3)=V(3,3);
clear V
%% This convention satisfies: M*n(1,:)=lambda(1)*n(1,:)
%% (idem to 2 y 3).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here the eigenvalues and unit vector are ordered
[lambda_ordenado,indices_lambda]=sort(lambda);
lambda_ordenado;
n_min=n(indices_lambda(1),:);
n_int=n(indices_lambda(2),:);
n_max=n(indices_lambda(3),:);
eval1=lambda_ordenado(1);
eval2=lambda_ordenado(2);
eval3=lambda_ordenado(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ============================================================================
# FUNCTIONS
# ============================================================================
"""
