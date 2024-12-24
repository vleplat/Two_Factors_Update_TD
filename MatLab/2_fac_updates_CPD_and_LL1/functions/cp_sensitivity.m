function sensitivity = cp_sensitivity(CPtensor)
    N = length(CPtensor.factors);
    
    % Normalize the CP tensor
    CPtensor = cp_normalize(CPtensor); % Ensure you have a cp_normalize function implemented
    
    % Scale factors by the weights
    Factors = cell(1, N);
    for n = 1:N
        Factors{n} = CPtensor.factors{n} * diag(CPtensor.weights.^(1/N));
    end
    
    % Compute C matrices
    C = cell(1, N);
    for n = 1:N
        C{n} = Factors{n}' * Factors{n};
    end
    
    % Initialize sensitivity
    sensitivity = 0;
    
    for n = 1:N
        % Create an array of numbers from 1 to N
        array = 1:N;
        
        % Remove the element at index n
        array(n) = []; % MATLAB indexing is 1-based, so we adjust accordingly
        
        % Compute the product of matrices in C
        product = C{array(1)};
        for k = 2:length(array)
            product = product * C{array(k)};
        end
        
        % Update sensitivity
        sensitivity = sensitivity + trace(product);
    end
end
