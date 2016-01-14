% Copyright (c) 2009 Sebastian Steiger, Integrated Systems Laboratory, ETH Zurich.
% Comments, suggestions, criticism or bug reports are welcome: steiger@purdue.edu. 
%
% This file is part of ANGEL, a simulator for LEDs based on the NEGF formalism.
% The software is distributed under the Lesser GNU General Public License (LGPL).
% ANGEL is free software: you can redistribute it and/or modify it under the terms 
% of the Lesser GNU General Public License v3 or later. ANGEL is distributed
% without any warranty; without even the implied warranty of merchantability or 
% fitness for a particular purpose. See also <http://www.gnu.org/licenses/>.

function mat = negf_read_real(file)

% read in an AQUA .xy-file
%file = '/home/steiger/src/release/amd64/aqua/results/qwr-kapon3_300K_tww3e13_kp8x8/qwr-kapon3_Source1.375V.xy';
fprintf(1,'Reading %s\n',file);

mat = load(file);
format short e;
fprintf('   matrix has %d columns and %d rows.\n',size(mat,2),size(mat,1));

