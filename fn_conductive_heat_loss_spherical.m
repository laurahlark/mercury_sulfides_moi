function [r_all, q, T] = fn_conductive_heat_loss_spherical(T_top, T_bottom, k, H, r)
    % Gives heat flux and temperature profile through spherical shell with
    % fixed temperature boundary conditions. Shell is described in n
    % layers, each of which may have a different thermal conductivity and
    % rate of internal heating.
    %
    % Parameters:
    %   k    W/mK   array of n elements specifying thermal conductivity of 
    %               each layer.
    %   H    W/m^3  array of n elements specifying internal heating 
    %               generated in each layer.
    %   r    m      array of n+1 elements specifying boundaries of the n
    %               layers. r(1) should be inner radius of the total shell,
    %               r(n+1) should be outer radius.
    %
    % Returns:
    %   r_all  m      array giving radii at which other return values are
    %                 sampled. r_all(1) = r(1), r_all(end) = r(n+1).
    %   q      W/m^2  array same size as r_all giving heat flux profile.
    %   T      K      array same size as r_all giving temperature profile.
    %                 T(1) = T_bottom, T(end) = T_top.
    
    % Some naming and renaming.
    n = length(k);
    T_n = T_top;
    T_0 = T_bottom;
    r_0 = r(1);
    r = r(2:end);
    
    % Compute A constants.
    A = zeros(1, n);
    
    % First, compute A1 from known values.
    A(1) = (T_n - T_0 + H(n)*r(n)^2/(6*k(n)) - H(1)*r_0^2/(6*k(1)) ...
        - 1/(3*k(n)*r(n)) * sum(r(1:n-1).^3 .* (H(1:n-1) - H(2:n))) ...
        + sum(r(1:n-1).^2/6 .* (H(1:n-1)./k(1:n-1) - H(2:n)./k(2:n))) ...
        + sum(r(1:n-1).^2.*(H(1:n-1)-H(2:n))./(3*k(1:n-1))) ...
        + sum(1./(3.*r(1:n-1)).*(1./k(2:n) - 1./k(1:n-1)) .* ...
            cumsum(r(1:n-1).^3 .* (H(1:n-1) - H(2:n))))) ...
        / (-1/(k(n)*r(n)) + 1/(k(1)*r_0) ...
        - sum((1./r(1:n-1)) .* (1./k(1:n-1) - 1./k(2:n))));
    
    % Use A1 to compute A2-An.
    A(2:n) = A(1) - cumsum(r(1:n-1).^3/3 .* (H(1:n-1) - H(2:n)));
    
    % Compute B constants.
    B = zeros(1, n);
    B(1) = T_0 + H(1)*r_0^2/(6*k(1)) + A(1)/(k(1)*r_0);
    B(2:n) = B(1) - cumsum(...
        r(1:n-1).^2/6 .* (H(1:n-1)./k(1:n-1) - H(2:n)./k(2:n)) + ...
        1./r(1:n-1) .* (A(1:n-1)./k(1:n-1) - A(2:n)./k(2:n)));
    
    % Checks:
    check = true;
    if check
        eps = .0000001;
        % Try building A1 up slowly.
        A1_numerator = T_n - T_0 + H(n)*r(n)^2/(6*k(n)) - H(1)*r_0^2/(6*k(1));
        for i=1:n-1
            A1_numerator = A1_numerator - 1/(3*k(n)*r(n)) * r(i)^3 * (H(i) - H(i+1));
            A1_numerator = A1_numerator + r(i)^2/6 *(H(i)/k(i) - H(i+1)/k(i+1));
            A1_numerator = A1_numerator + r(i)^2*(H(i)-H(i+1))/(3*k(i));
            alt_inner_sum = 0;
            for j=1:i
                alt_inner_sum = alt_inner_sum + r(j)^3*(H(j)-H(j+1));
            end
            A1_numerator = A1_numerator + 1/(3*r(i))*(1/k(i+1)-1/k(i))*alt_inner_sum;
        end
        A1_denominator = -1/(k(n)*r(n)) + 1/(k(1)*r_0);
        for i=1:n-1
            A1_denominator = A1_denominator - 1/r(i) * (1/k(i) - 1/k(i+1));
        end
        if abs((A(1) - A1_numerator/A1_denominator)/A(1)) > eps
            fprintf("A1 value other than expected: %.9g, %.9g\n", A(1), A1_numerator/A1_denominator);
        end
        Aip1 = A(1:n-1) - r(1:n-1).^3/3.*(H(1:n-1) - H(2:n));
        if sum(abs(A(2:n) - Aip1)) / sum(Aip1) > eps
            fprintf("A values other than expected:\n");
            A(2:n)
            Aip1
        end
        Bip1 = B(1:n-1) - r(1:n-1).^2/6.*(H(1:n-1)./k(1:n-1) - H(2:n)./k(2:n)) - 1./r(1:n-1).*(A(1:n-1)./k(1:n-1) - A(2:n)./k(2:n));
        if sum(abs(B(2:n) - Bip1)) > eps
            fprintf("B values other than expected:\n");
            B
            Bip1
        end
        Bn = T_n + H(n)/(6*k(n))*r(n)^2 + A(n)/(k(n)*r(n));
        if abs(B(n) - Bn)/B(n) > eps
            fprintf("B_n values other than expected: %.15g, %.15g\n", B(n), Bn);
        end

        if n == 1
            Acheck = (T_n - T_0 - H(1)/(6*k(1)) * (r_0^2 - r(1)^2)) / (1/(k(1)*r_0) - 1/(k(1)*r(1)));
            if abs((A(1) - Acheck)/A(1)) > eps
                fprintf("A does not match expected for single layer: %g, %g\n", ...
                    A(1), Acheck);
            end
            Bcheck = T_0 + H(1)/(6*k(1))*r_0^2 + (T_n - T_0 - H(1)/(6*k(1)) * (r_0^2-r(1)^2)) / (1 - r_0/r(1));
            if abs((B(1) - Bcheck)/B(1)) > eps
                fprintf("B does not match expected for single layer: %g, %g\n", ...
                    B(1), Bcheck);
            end
        end
    end
    
    % Compute r_all, q, T.
    r_all = sort(unique([r, r_0:1e3:r(n)]));  % Every km, but don't miss values in r.
    q = zeros(size(r_all));
    T = zeros(size(r_all));
    for i=1:n
        if i == 1; r_lo = r_0; else; r_lo = r(i-1); end
        q_match = q(r_all == r_lo);
        T_match = T(r_all == r_lo);
        all_i = find(r_all >= r_lo & r_all <= r(i));
        q(all_i) = H(i) * r_all(all_i) / 3 - A(i) ./ r_all(all_i).^2;
        T(all_i) = -H(i)/(6*k(i)) .* r_all(all_i).^2 - A(i)./(k(i)*r_all(all_i)) + B(i);
        
        % Check matching at layer borders.
        if check
            if i > 1
                if abs(q_match - q(r_all == r_lo)) > eps
                    fprintf("Heat flux mismatch at %g: %g vs %g\n", ...
                        r(i), q_match, q(r_all == r_lo));
                end
                if abs(T_match - T(r_all == r_lo)) > eps
                    fprintf("Temperature mismatch at r=%g km: %g vs %g\n", ...
                        r(i)/1e3, T_match, T(r_all == r_lo));
                end
            end
        end
    end
    
    if check
        if abs(T(1) - T_0) > eps
            fprintf("Bottom temperature does not match: %g, %g\n", T_0, T(1));
        end
        if abs(T(end) - T_n) > eps
            fprintf("Top temperature does not match: %g, %g\n", T_n, T(end));
        end
    end
end


function out=sum_so_far(vector)
    out = cumsum(vector);
    % value at i gives sum of vector(1:i).
%     out = transpose(sum(tril(ones(length(vector), 1) .* vector), 2));
end