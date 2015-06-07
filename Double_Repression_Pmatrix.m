% ------------------------------------------------------------------------------------- %
% Copyright (c) 2015 Lucks Lab,
% School of Chemical and Biomolecular Engineering,
% Cornell University, Ithaca NY 14853 USA.
%
% Details about this software can be found in 
% Hu et al. “Generating effective models and parameters for RNA genetic circuits” 
% ACS SynBio, 2015, http://dx.doi.org/10.1021/acssynbio.5b00077 .
%
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions: 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%
% Double_Repression_Pmatrix.m
% A function that solves the P_matrix analytically with given ODE and parameters 
% ------------------------------------------------------------------------------------- %


function [ P_mtrx ] = Double_Repression_Pmatrix(P,C,x,nstep)

P_mtrx=zeros(8,15,(nstep+1)); % initialize the Jacobian matrix 

for i=2:(nstep+1)
 
    %taking x values starting from the 2nd row   
    P_mtrx(1,1,i)=1;
    P_mtrx(3,2,i)=(1-x(i,2)./(P(4)+x(i,2)))^C(4).*(1-P(15)).^C(4);
    P_mtrx(5,3,i)=(1+(-1).*P(15)).^C(3).*(1+(-1).*x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-1)+(-1).* ...
    x(i,4).*(P(5)+x(i,4)).^(-1)+x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-1).*x(i,4).*(P(5)+x(i,4)).^(-1)) ...
    .^C(3);
    P_mtrx(3,4,i)=((C(4).*x(i,2)).*P(2).*(1-P(15)).^C(4).*(P(4)./(P(4)+x(i,2))).^(C(4)-1))./(P(4)+x(i,2)).^2;
    P_mtrx(5,5,i)=C(3).*(1+(-1).*P(15)).^C(3).*P(3).*(P(14).^(-1).*x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^( ...
    -2)+x(i,4).*(P(5)+x(i,4)).^(-2)+(-1).*x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-1).*x(i,4).*(P(5)+ ...
    x(i,4)).^(-2)+(-1).*P(14).^(-1).*x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-2).*x(i,4).*(P(5)+ ...
    x(i,4)).^(-1)).*(1+(-1).*x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-1)+(-1).*x(i,4).*(P(5)+ ...
    x(i,4)).^(-1)+x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-1).*x(i,4).*(P(5)+x(i,4)).^(-1)).^((-1)+ ...
    C(3));
    P_mtrx(7,6,i)=x(i,6);
    P_mtrx(6,6,i)=-x(i,6);
    P_mtrx(1,7,i)=-x(i,1);
    P_mtrx(2,7,i)=-x(i,2);
    P_mtrx(3,8,i)=-x(i,3);
    P_mtrx(4,8,i)=-x(i,4);
    P_mtrx(5,9,i)=-x(i,5);
    P_mtrx(6,9,i)=-x(i,6);
    P_mtrx(7,10,i)=-x(i,7);
    P_mtrx(8,10,i)=x(i,7);
    P_mtrx(5,11,i)=-x(i,5);
    P_mtrx(6,11,i)=x(i,5);
    P_mtrx(1,12,i)=-x(i,1);
    P_mtrx(2,12,i)=x(i,1);
    P_mtrx(3,13,i)=-x(i,3);
    P_mtrx(4,13,i)=x(i,3);
    P_mtrx(5,14,i)=C(3).*(1+(-1).*P(15)).^C(3).*P(3).*((-1).*P(14).^(-2).*P(5).*x(i,2).*(P(14).^(-1).* ...
    P(5)+x(i,2)).^(-2)+P(14).^(-2).*P(5).*x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-2).*x(i,4).*(P(5)+ ...
    x(i,4)).^(-1)).*(1+(-1).*x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-1)+(-1).*x(i,4).*(P(5)+ ...
    x(i,4)).^(-1)+x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-1).*x(i,4).*(P(5)+x(i,4)).^(-1)).^((-1)+ ...
    C(3));
    P_mtrx(5,15,i)=(-1).*C(3).*(1+(-1).*P(15)).^((-1)+C(3)).*P(3).*(1+(-1).*x(i,2).*(P(14).^(-1).* ...
    P(5)+x(i,2)).^(-1)+(-1).*x(i,4).*(P(5)+x(i,4)).^(-1)+x(i,2).*(P(14).^(-1).*P(5)+x(i,2)).^(-1) ...
    .*x(i,4).*(P(5)+x(i,4)).^(-1)).^C(3);
    P_mtrx(3,15,i)=P(2)*C(4)*(-(1-P(15))^(C(4)-1))*(P(4)/(P(4)+x(i,2)))^(C(4));
end

end

