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
% fmincon_constrains.m
% This function contains the constraints that are used in fmincon
% ------------------------------------------------------------------------------------- % 

function [ A,b ] = fmincon_constrains(pset,P)
%This function contains fmincon constraints

Est_num=length(pset);
A=zeros(2*(Est_num),(Est_num));
b=zeros(2*(Est_num),1);

for i=1:Est_num
    x=pset(i);
    if x==4 %P4 & P5 are in the same order of magnitude
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=3.33*P(5);
        b(2*i)=-0.33*P(5);
    elseif x==5%P4 & P5 are in the same order of magnitude
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=3.33*P(4);
        b(2*i)=-0.33*P(4);       
    elseif x==6 %elongation rate <= initiation rate 
        A(2*i-1,i)=1;
        b(2*i-1)=P(11);
    elseif x==7 %P7 & P8 are 50% around each other
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=2*P(8);
        b(2*i)=-0.5*P(8);
    elseif x==8 %P7 & P8 are 50% around each other
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=2*P(7);
        b(2*i)=-0.5*P(7);        
    elseif x==9 %dagradation rate of mRNA<degradation rate of antisense RNA
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=P(7);
        b(2*i)=-0.0001;       
    elseif x==11 % elongation rate <= initiation rate <= 10*elongation rate  
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=10*P(6);
        b(2*i)=-P(6);
    elseif x==12 %P12 & P13 in the same order of magnitude
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=3.33*P(13);
        b(2*i)=-0.33*P(13);
    elseif x==13 %P12 & P13 in the same order of magnitude
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=3.33*P(12);
        b(2*i)=-0.33*P(12); 
    elseif x==14 %0.1%~20% crosstalk
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=0.2;
        b(2*i)=-0.001;         
    elseif x==15 % 0.1%~70% autotermination
        A(2*i-1,i)=1;
        A(2*i,i)=-1;
        b(2*i-1)=0.7;
        b(2*i)=-0.001;  
    end
end
end