% RBR                                                                             
%                                                                                 
% Copyright (C) 2019  Haoyang Liu (liuhaoyang@pku.edu.cn)                         
%                     Zaiwen Wen  (wenzw@pku.edu.cn)                              
%                                                                                 
% This program is free software: you can redistribute it and/or modify            
% it under the terms of the GNU General Public License as published by            
% the Free Software Foundation, either version 3 of the License, or               
% (at your option) any later version.                                             
%                                                                                 
% This program is distributed in the hope that it will be useful,                 
% but WITHOUT ANY WARRANTY; without even the implied warranty of                  
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   
% GNU General Public License for more details.                                    
%                                                                                 
% You should have received a copy of the GNU General Public License               
% along with this program.  If not, see <http://www.gnu.org/licenses/>.  

%% A demo matlab file to invoke RBR MEX interface

% load matrix
load ../../examples/polblogs.mat
d = sum(A, 2);
opts = struct();
opts.verbose = 1;
opts.extract = 'rounding';
opts.full = 1;
k = 2;

% call rbr
id = rbr(A, d, k, opts);