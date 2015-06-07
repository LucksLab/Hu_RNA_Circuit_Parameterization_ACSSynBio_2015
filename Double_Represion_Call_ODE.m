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
% Double_Represion_Call_ODE.m
% Contains ODEs for the system 
% ------------------------------------------------------------------------------------- %


function [ t, x ] = Double_Represion_Call_ODE(DF,tspan,P,CSTR_LV)
%This function solves the system 

x0=DF.Initial_Conditions;
C=DF.Constants; 
    function dxdt=MassbalanceEqns(t,x)
            
        if CSTR_LV==3
            
            dxdt(1,1) = P(1)-P(12).*x(1)-P(7).*x(1);
            dxdt(2,1) = P(12).*x(1)-P(7).*x(2);
            dxdt(3,1) = P(2).*((1-x(2)./(P(4)+x(2))).^C(4)).*(1-0).^C(4)-P(13).*x(3)-P(8).*x(3); %No autotermination in this exp. setup
            dxdt(4,1) = P(13).*x(3)-P(8).*x(4);
            dxdt(5,1) = P(3).*(1-((x(4)./(P(5)+x(4)))+(x(2)./((P(5)./P(14))+x(2)))-((x(4)./(P(5)+x(4)))).*(x(2)./((P(5)./P(14))+x(2))))).^C(3).*(1-P(15)).^C(3)-P(11).*x(5)-P(9).*x(5)+P(6).*x(6);
            dxdt(6,1) = P(11).*x(5)-P(6).*x(6);
            dxdt(7,1) = P(6).*x(6)-P(10).*x(7);
            dxdt(8,1) = P(10).*x(7);
         
        else
            dxdt(1,1) = P(1)-P(12).*x(1)-P(7).*x(1);
            dxdt(2,1) = P(12).*x(1)-P(7).*x(2);
            dxdt(3,1) = P(2).*((1-x(2)./(P(4)+x(2))).^C(4)).*(1-P(15)).^C(4)-P(13).*x(3)-P(8).*x(3);
            dxdt(4,1) = P(13).*x(3)-P(8).*x(4);
            dxdt(5,1) = P(3).*(1-((x(4)./(P(5)+x(4)))+(x(2)./((P(5)./P(14))+x(2)))-((x(4)./(P(5)+x(4)))).*(x(2)./((P(5)./P(14))+x(2))))).^C(3).*(1-P(15)).^C(3)-P(11).*x(5)-P(9).*x(5)+P(6).*x(6);
            dxdt(6,1) = P(11).*x(5)-P(6).*x(6);
            dxdt(7,1) = P(6).*x(6)-P(10).*x(7);
            dxdt(8,1) = P(10).*x(7);
        end
            

    end
    
    
[t,x]=ode15s(@MassbalanceEqns,tspan,x0);

end