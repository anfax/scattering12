#!/usr/bin/csh
ifc -o integr_1D.x -check bounds -check overflow -check underflow -lcxml integr_1D.K_mat.f90 adda.f90 pleg.f90 besselnew.f90
