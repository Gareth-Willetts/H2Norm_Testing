function [p1, p2] = get_p1_p2(n, a)
% GET_P1_P2
%
% @brief:
%   This function implements the logic from Function 2. It takes
%   a polynomial a(s), represented as a vector of coefficients in descending
%   order of power, and splits it into two new polynomials, p1(s) and p2(s).
%
%   The algorithm constructs one polynomial, a^e(s), from the coefficients of
%   the even-powered terms of a(s), and another, a^o(s), from the
%   odd-powered terms. The final assignment to p1(s) and p2(s) depends on the
%   degree of the input polynomial a(s).
%
% Inputs:
%   a - A row vector representing the coefficients of the input polynomial
%       in descending order of power (e.g., [3 2 1] for 3x^2 + 2x + 1).
%
% Outputs:
%   p1 - The first output polynomial (coefficient vector).
%   p2 - The second output polynomial (coefficient vector).
%
% Example 1: n is even (n=4)
%   % Let a(x) = 5x^4 + 4x^3 + 3x^2 + 2x + 1
%   a = [5 4 3 2 1];
%   n = 4;
%   [p1, p2] = get_p1_p2(n, a);
%   % Expected output:
%   % p1 = [5 3 1]   (from a^e)
%   % p2 = [4 2 0]   (from a^o)
%
% Example 2: n is odd (n=5)
%   % Let a(x) = 6x^5 + 5x^4 + 4x^3 + 3x^2 + 2x + 1
%   a = [6 5 4 3 2 1];
%   n = 5;
%   [p1, p2] = get_p1_p2(n, a);
%   % Expected output:
%   % p1 = [6 4 2 0] (from a^o)
%   % p2 = [5 3 1]   (from a^e)

    % The pseudocode builds new polynomials ae(s) and ao(s).
    % ae(s) = a_0*s^0 + a_2*s^1 + a_4*s^2 + ...
    % ao(s) = a_1*s^1 + a_3*s^2 + a_5*s^3 + ...
    
    % To implement this efficiently in MATLAB, we extract the coefficients
    % of a(x) corresponding to even and odd powers.
    
    % If a = [a_n, a_{n-1}, ..., a_1, a_0], then the even-powered
    % coefficients (a_0, a_2, ...) can be found by starting from the end
    % of the vector and taking every second element.
    ae_coeffs = a(end:-2:1);
    
    % Similarly, the odd-powered coefficients (a_1, a_3, ...) are found
    % by starting from the second-to-last element.
    ao_prime_coeffs = a(end-1:-2:1);
    
    % For a^e(s), the power of s is k, corresponding to the original power
    % 2k. This is a direct mapping, so we just reverse the extracted
    % coefficients to get the standard descending power representation.
    ae = fliplr(ae_coeffs);

    % For a^o(s), the power of s is k+1 for the original coefficient a_{2k+1}.
    % This means a^o(s) is s * (a_1*s^0 + a_3*s^1 + ...).
    % Multiplying by 's' in polynomial terms is equivalent to appending a 
    % zero to the coefficient vector.
    if isempty(ao_prime_coeffs)
        % Handle cases where 'a' has no odd-powered terms.
        ao = 0;
    else
        ao = [fliplr(ao_prime_coeffs), 0];
    end

    % Final assignment based on the parity of the degree n.
    if mod(n, 2) == 1 % n is odd
        p1 = ao;
        p2 = ae;
    else % n is even
        p1 = ae;
        p2 = ao;
    end

end