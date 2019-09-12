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


function idx = rbr(A, d, k, opts)
    init;
    % check parameters
    if opts.p > k
        warning('p is set to k.');
        opts.p = k;
    end
    if opts.maxit < 5
        error('Max # of iterations is at least 5.');
    end
    
    if opts.verbose
        fprintf('Invoking rbr solver ... patient ...\n');
    end
    if strcmp(opts.extract, 'rounding')
        idx = mex_rbr(A, d, k, opts.p, opts.maxit, 0, opts.full);
        idx = idx + 1;
    elseif strcmp(opts.extract, 'kmeans')
        [~, U] = mex_rbr(A, d, k, opts.p, opts.maxit, 2, opts.full);
        % perform K-means
        if opts.verbose
            fprintf('Invoking K-means ... patient ...\n');
        end
        [idx, ~] = kmeans(U, k);
    else
        error(['unknown extract type: ', opts.extract]);
    end
    U = zeros(size(A, 1), k);
    for i = 1:size(U, 1)
        U(i, idx(i)) = 1;
    end
    if opts.verbose
        final_Q = cal_Q(A, d, 1 / sum(d), U);
        fprintf('Finished. Final Q: %f\n', final_Q);
    end
    
    function init
        if ~isfield(opts, 'p'); opts.p = 5; end
        if ~isfield(opts, 'maxit'); opts.maxit = 25; end
        if ~isfield(opts, 'verbose'); opts.verbose = 0; end
        if ~isfield(opts, 'extract'); opts.extract = 'rounding'; end
        if ~isfield(opts, 'full'); opts.full = 0; end
    end
end
