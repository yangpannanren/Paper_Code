function [ numbers ] = laprnd( num )
%LAPRND Summary of this function goes here
%   Generates random numbers according to the Laplacian distribution.
%   Standard deviation is set to 1.

    temp = exprnd(1/sqrt(2),num,2);
    numbers = temp(:,1)-temp(:,2);

end
