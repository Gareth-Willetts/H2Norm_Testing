function H2ns = numerical_H2_norm(a, c)
% SOLVE_FRACTION_FREE_H2_NORM 
%
% @brief:
%   Adjusted to implement Algorithm 1: Numerical computation of the H2 norm.
%   This calculates the square of the H2 norm.
%
% Input:
%   a - Row vector of denominator coefficients (descending powers).
%   c - Row vector of numerator coefficients (descending powers).
%
% Output:
%   H2ns    - The calculated square of the H2 Norm (Result).

    % n <- deg(a)
    n = length(a) - 1;
    [p1, p2] = get_p1_p2(n, a);
    
    % Store initial leading coefficients for the final H2 formula
    % coeff(p1, deg(p1)) and coeff(p2, deg(p2))
    coeff_p1 = p1(1);
    coeff_p2 = p2(1);

    % --- Construct polynomials ce and co from c(s) ---
    % The pseudocode creates ce(s) from even-indexed coefficients of c(s)
    % and co(s) from odd-indexed coefficients.
    % ce(s) = c_0*s^0 + c_2*s^1 + c_4*s^2 + ...
    % co(s) = c_1*s^0 + c_3*s^1 + c_5*s^2 + ...
    c_rev = fliplr(c);
    ce_coeffs_ascending = c_rev(1:2:end);
    co_coeffs_ascending = c_rev(2:2:end);
    
    ce = fliplr(ce_coeffs_ascending);
    co = fliplr(co_coeffs_ascending);
    
    % z0(s) <- (ce(s))^2 - s * (co(s))^2
    ce_squared = conv(ce, ce);
    co_squared = conv(co, co);
    s_times_co_squared = [co_squared, 0];
    
    % Pad with zeros to subtract polynomials of different lengths
    len1 = length(ce_squared);
    len2 = length(s_times_co_squared);
    max_len = max(len1, len2);
    
    ce_squared_padded = [zeros(1, max_len - len1), ce_squared];
    s_co_sq_padded = [zeros(1, max_len - len2), s_times_co_squared];
    
    z0 = ce_squared_padded - s_co_sq_padded;
    % ----------------------------------------------------
    % d1 <- deg(p1), d2 <- deg(p2)
    d1 = length(p1) - 1;
    d2 = length(p2) - 1;
    
    % Stability check setup
    s_flag = 1;
    if isempty(p2) || p2(1) == 0 || p1(1) * p2(1) <= 0
        s_flag = 0;
    end
    
    % Call the updated algorithm implementation
    % We pass the initial coefficients needed for the final formula
    [s_flag, pn_plus_1_val, zn_minus_1_val, H2_result] = ...
        get_algorithm1_terms(p1, p2, z0, d1, d2, s_flag, coeff_p1, coeff_p2);
    
    if s_flag
        H2ns = H2_result; % Return the computed H2 norm square
    else
        H2ns = [];
    end
end

function [p1, p2] = get_p1_p2(n, a)
% GET_P1_P2
    ae_coeffs = a(end:-2:1);
    ao_prime_coeffs = a(end-1:-2:1);
    
    ae = fliplr(ae_coeffs);
    if isempty(ao_prime_coeffs)
        ao = 0;
    else
        ao = [fliplr(ao_prime_coeffs), 0];
    end
    
    if mod(n, 2) == 1 % n is odd
        p1 = ao;
        p2 = ae;
    else % n is even
        p1 = ae;
        p2 = ao;
    end
end

function [s, pb_out, z_out, H2] = get_algorithm1_terms(p1, p2, z0, da, db, s, cp1, cp2)
% GET_ALGORITHM1_TERMS
% Implements the loop and logic from Algorithm 1.
    
    % Initialisation according to Algorithm 1
    pa = p1;
    pb = p2;
    z  = z0;
    n_sum = da + db; % n in the algorithm
    q = mod(n_sum, 2);
    
    % --- Initial steps ---
    % Here we choose F = LC() (leading coefficient)
    % vz <- F1(z)
    vz = z(1); 
    if abs(vz) < 1e-12
        s = 0; 
    end % Stability check against divide by zero
    
    % z <- z / vz
    if s
        z = z / vz; 
    end
    
    % h <- vz (Combined variable for mu and nu)
    h = vz;
    
    % vp <- F2(pb)
    vp = pb(1);
    if abs(vp) < 1e-12
        s = 0;
    end
    
    % pb <- pb / vp
    if s
        pb = pb / vp;
    end
    
    % h <- h / vp
    if s
        h = h / vp;
    end
    
    k = 2;

    % --- Main Loop ---
    while (k <= n_sum) && s
        % Update z
        % s^(floor((n-k)/2)) aligns degrees. 
        % In MATLAB, we append zeros to pa to align it with z.
        shift_z = floor((n_sum - k)/2);
        
        % The coefficient ratio: coeff(z, n+1-k) / coeff(pa, da)
        % This is z(1) / pa(1)
        % if abs(pa(1)) < 1e-12, s = 0; break; end
        ratio_z = z(1) / pa(1);
        
        pa_aligned = [pa, zeros(1, shift_z)];
        
        % z <- z - ratio * pa
        z = poly_subtract(z, ratio_z * pa_aligned);
        
        % Find LC again
        vz = z(1);
        % if abs(vz) < 1e-12
        %      s = 0; break; 
        % end
        z = z / vz;
        h = h * vz;
        
        % Update pc (which replaces pb)
        % coeff(pa, da) / coeff(pb, db)
        % if abs(pb(1)) < 1e-12, s = 0; break; end
        ratio_p = pa(1) / pb(1);
        
        pb_aligned = [pb, zeros(1, q)];
        
        % pc <- pa - s^q * ratio * pb
        pc = poly_subtract(pa, ratio_p * pb_aligned);
        
        % Find LC again
        vp = pc(1);
        % if abs(vp) < 1e-12, s = 0; break; end
        pc = pc / vp;
        h = h / vp;
        
        % Shift and update indices
        pa = pb;
        pb = pc;
        
        k = k + 1;
        q = ~q; % Toggle q
        da = db;
        if q
            db = db - 1;
        end
        
        % (Optional: explicit stability check, will implement once this is correct)
    end
    
    if s
        pb_out = pb;
        z_out  = z;
        
        % --- Calculate Result: H2 norm squared ---
        % Formula: (h * coeff(p1) * z) / (2 * coeff(pa) * coeff(p2) * pb)

        numerator = h * cp1 * z; 
        denominator = 2 * pa(1) * pb;
        
        H2 = numerator / denominator;
    else
        pb_out = [];
        z_out = [];
        H2 = [];
    end
end

function p_out = poly_subtract(poly1, poly2)
    % Subtracts poly2 from poly1.
    % Padding is done in the main while loop.
    len1 = length(poly1);
    len2 = length(poly2);
    if len1 == len2
        p_out = poly1 - poly2;
    else
        s = 0;
        error('Mismatched polynomial lengths.');
    end
    
    % Remove leading zeros from the result to keep polynomials tidy.
    first_nonzero = find(p_out ~= 0, 1, 'first');
    if isempty(first_nonzero)
        p_out = 0; % The polynomial is the scalar zero.
    else
        p_out = p_out(first_nonzero:end);
    end
end