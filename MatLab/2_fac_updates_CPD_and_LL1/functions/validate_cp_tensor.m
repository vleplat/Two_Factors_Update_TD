function [shape, rank] = validate_cp_tensor(cp_tensor)
    % Validates a cp_tensor in the form (weights, factors)
    % Returns the rank and shape of the validated tensor

    % Check if cp_tensor is a structure (CPTensor)
    if isfield(cp_tensor, 'weights') && isfield(cp_tensor, 'factors')
        % It's already been validated at creation
        shape = cp_tensor.shape;
        rank = cp_tensor.rank;
        return;
    elseif isnumeric(cp_tensor) && isscalar(cp_tensor) % 0-order tensor
        shape = 0;
        rank = 0;
        return;
    end

    weights = cp_tensor.weights;
    factors = cp_tensor.factors;

    % Determine rank based on the first factor
    if ndims(factors{1}) == 2
        rank = size(factors{1}, 2);
    elseif ndims(factors{1}) == 1
        rank = 1;
    else
        error('Got a factor with 3 dimensions but CP factors should be at most 2D, of shape (size, rank).');
    end

    shape = zeros(1, length(factors)); % Preallocate shape array
    for i = 1:length(factors)
        s = size(factors{i});
        if length(s) == 2
            current_mode_size = s(1);
            current_rank = s(2);
        else % The shape is just (size, ) if rank 1
            current_mode_size = s(1);
            current_rank = 1;
        end

        if current_rank ~= rank
            error('All the factors of a CP tensor should have the same number of columns. However, factors{1}.shape(2)=%d but factors{%d}.shape(2)=%d.', rank, i, size(factors{i}, 2));
        end

        shape(i) = current_mode_size; % Store the current mode size in shape
    end

    if ~isempty(weights) && ~isequal(size(weights), [rank, 1])
        error('Given factors for a rank-%d CP tensor but len(weights)=%d.', rank, size(weights, 1));
    end

    shape = num2cell(shape); % Convert shape to a cell array if needed
    shape = cell2mat(shape); % If you want it as a numeric array
end
