function [HilbD] = HilbDim(NSites,NPart)

% combinatorics for choosing how to split up NPart Bosons on NSites
% also sometimes thought of as the balls and sticks sorting, where the
% balls are the bosons and the sticks are the way they separated into
% lattice sites. So for NSites you have NSites-1 dividers, hence the
% choosing by NSites-1 = numer of sticks

HilbD = round(factorial(NPart+NSites-1)/factorial(NPart)/factorial(NSites-1));

end
