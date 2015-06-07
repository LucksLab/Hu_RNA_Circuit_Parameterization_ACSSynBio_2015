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
% Datafile.m
% Contains a struct of data that's used in Main scripts.
% ------------------------------------------------------------------------------------- %


function [ DF ] = Datafile(CSTR_LV )

iter_num=1000;

%initial guess of all parameters
file=fopen('guess.txt');
P_i= fscanf(file, '%f',[1,15])'; 


P_i(1)=28*P_i(3);
P_i(2)=8*P_i(3);
pset=[3;6;9;10;11]; % identifiable parameters at LV1
if CSTR_LV==1
    P_i(1)=0;
    P_i(2)=0;
elseif CSTR_LV==2
    P_i(1)=0;
    
elseif CSTR_LV==3    % parameter shifting 
    P_i(1)=0;
    P_i(5)=P_i(4);
    P_i(8)=P_i(7);
    P_i(13)=P_i(12);

end


IC=[0           %A'2_0 
    0           %A2_0
    0           %A'1_0 
    0           %A1_0
    0           %M_0 
    0           %Mi_0
    0           %G_0 
    0];         %GM_0 

C=[ 0       
    0            
    1               % Level 1 order
    1] ;            % Level 2 order


P_i_lb=P_i.*0.85;    %artificial range, within 15% of the best guess
P_i_ub=P_i.*1.15;       
NParameters=length(P_i);
Nstates=length(IC);
 
a=P_i_lb;            %generate a set of P from 15% range about primary guess, first round, all random
 b=P_i_ub ;
 for j=1:(iter_num)
     for i=1:(NParameters)
         P(j,i)=a(i)+(b(i)-a(i)).*rand(1,1);
     end
 end
 
DF.Initial_Parameters=P_i;
DF.Initial_Conditions=IC;
DF.Constants=C;
DF.Num_States=Nstates;
DF.Num_Parameters=NParameters;
DF.Construct=CSTR_LV;
DF.P_i_LB=P_i_lb;
DF.P_i_UB=P_i_ub;
DF.Parameter_library=P;
DF.iter_num=iter_num;
DF.pset=pset;
return;
