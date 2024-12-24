function kr_product = khatri_rao(matrices)
    % Get the size of the first matrix
    [~, n_columns] = size(matrices{1});
    prod = 1;
    for k=1:length(matrices)
        prod = prod*size(matrices{k},1);
    end

    
    % Initialize the Khatri-Rao product matrix
    kr_product = zeros(prod, n_columns);
    
    % Loop over each column
    for i = 1:n_columns
        cum_prod = matrices{1}(:, i); % Initialize with the i-th column of the first matrix
        
        % Accumulate the Khatri-Rao product for the i-th column
        for j = 2:length(matrices)
            cum_prod = kron(cum_prod, matrices{j}(:, i)); % Kronecker product
        end
        
        % Assign the accumulated product to the i-th column of the result
        kr_product(:, i) = cum_prod;
    end
end