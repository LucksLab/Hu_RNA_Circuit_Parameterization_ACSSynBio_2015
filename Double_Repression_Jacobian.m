% ------------------------------------------------------------------------------------- %
% Copyright (c) 2015 Lucks Lab,
% School of Chemical and Biomolecular Engineering,
% Cornell University, Ithaca NY 14853 USA.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions: % The above copyright notice and this permission notice shall be included in
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
% Double_Repression_Jacobian.m
% A function that solves the Jacobian matrix analytically with given ODE and parameters 
% ------------------------------------------------------------------------------------- %


function [ J_mtrx ] = Double_Repression_Jacobian(P,C,x,nstep)

J_mtrx=zeros(8,8,(nstep+1));     % initialize the Jacobian matrix 

J_mtrx(1,1,:)=-P(12)-P(7);       % These elements are constants.
J_mtrx(2,1,:)=P(12);
J_mtrx(2,2,:)=-P(7);
J_mtrx(3,3,:)=-P(13)-P(8);
J_mtrx(4,3,:)=P(13);
J_mtrx(4,4,:)=-P(8);
J_mtrx(5,5,:)=-P(11)-P(9);
J_mtrx(6,5,:)=P(11);
J_mtrx(6,6,:)=-P(9)-P(6);
J_mtrx(7,6,:)=P(6);
J_mtrx(7,7,:)=-P(10);
J_mtrx(8,7,:)=P(10);


for i=2:(nstep+1)
    J_mtrx(3,2,i)=-(C(4).*P(2).*(1-P(15)).^C(4).*(P(4)./(P(4)+x(i,2))).^C(4))./(P(4)+x(i,2));

    J_mtrx(5,2,i)=-(C(3).*P(3).*(1-P(15)).^C(3).*((P(5)/P(14)).*P(5)./((x(i,2)+(P(5)/P(14))).*(x(i,4)+P(5)))).^C(3))./((P(5)/P(14))+x(i,2));

    J_mtrx(5,4,i)=-(C(3).*P(3).*(1-P(15)).^C(3).*((P(5)/P(14)).*P(5)./((x(i,2)+(P(5)/P(14))).*(x(i,4)+P(5)))).^C(3))./(P(5)+x(i,4));

end

end

