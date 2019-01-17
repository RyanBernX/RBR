function idx = rbr_kmeans(A, d, k, opts)
    init;
    
    if opts.verbose
        fprintf('Invoking rbr solver ... patient ...\n');
    end
    if strcmp(opts.extract, 'rounding')
        idx = mex_rbr(A, d, k, opts.p, opts.maxit, 0);
        idx = idx + 1;
    elseif strcmp(opts.extract, 'kmeans')
        [~, U] = mex_rbr(A, d, k, opts.p, opts.maxit, 2);
        % perform K-means
        if opts.verbose
            fprintf('Invoking K-means ... patient ...\n');
        end
        [idx, ~] = kmeans(U, k);

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
    end
end