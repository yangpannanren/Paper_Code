function [X,Y]=circle(center,radius,NOP)

if (nargin <3),
 error('Please see help for INPUT DATA.');
end;
THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X(:)+center(1);
Y=Y(:)+center(2);