function cp_tensor = cp_normalize(cp_tensor)
    % Returns cp_tensor with factors normalized to unit length
    % Parameters
    % ----------
    % cp_tensor : CPTensor = (weights, factors)
    %     factors is a cell array of matrices, all with the same number of columns
    %     where each factor is of shape (s_i, R)
    
    % Validate the CP tensor and get rank
    [~, rank] = validate_cp_tensor(cp_tensor); % Ensure you have a validate_cp_tensor function implemented
    weights = cp_tensor.weights;
    factors = cp_tensor.factors;

    if isempty(weights)
        weights = ones(rank, 1); % Create a weight vector of ones
    end

    normalized_factors = cell(1, length(factors));
    for i = 1:length(factors)
        factor = factors{i};
        if i == 1
            factor = factor .* weights'; % Scale the first factor by weights
            weights = ones(rank, 1); % Reset weights to ones
        end

        % Compute the norms (scales) of the columns
        scales = vecnorm(factor, 2, 1); % Compute the Euclidean norm along columns
        scales_non_zero = scales;
        scales_non_zero(scales == 0) = 1; % Replace zeros with ones to avoid division by zero

        weights = weights .* scales'; % Update weights
        normalized_factors{i} = factor ./ scales_non_zero; % Normalize factors
    end

    cp_tensor.weights = weights; % Update weights in cp_tensor
    cp_tensor.factors = normalized_factors; % Update factors in cp_tensor
end
