% ------------------------------------------------------------------------------------- %
% Copyright (c) 2015 Lucks Lab,
% School of Chemical and Biomolecular Engineering,
% Cornell University, Ithaca NY 14853 USA.
%
% Details about this software can be found in 
% Hu et al. “Generating effective models and parameters for RNA genetic circuits” 
% ACS SynBio, 2015, http://dx.doi.org/10.1021/acssynbio.5b00077 .
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
% Analysis_func.m
% Reads estimated parameters and passes them to Analysis_Main.m to recreate
% parameterization trajectries.
% ------------------------------------------------------------------------------------- %


function [DF] = Analysis_func(CSTR_LV)

iter_num=100;
filename=['P_solution.txt'];
file=fopen(filename);
P_est= fscanf(file, '%f',[iter_num,15]); 

if CSTR_LV==1    
 P_est(:,1)=0;
 P_est(:,2)=0;
elseif CSTR_LV==3
 
 P_est(:,8)=P_est(:,7);
 P_est(:,5)=P_est(:,4);
 P_est(:,13)=P_est(:,12);
 P_est(:,1)=0;
 
elseif CSTR_LV==2
 P_est(:,1)=0;

end


IC=[0           %A'2_0 
    0           %A2_0
    0           %A'1_0 
    0           %A1_0
    0           %M_0 
    0           %Mi_0
    0           %G_0 
    0];         %GM_0 




C=[ 0            % k_c,the crosstalk repression coefficient of A_2 on att-1
    0            % p_t,the probability of a attentuator auto-determinates itself
    1               % Level 1 order
    1] ;            % Level 2 order


DF.Num_Parameters=length(P_est(1,:));
DF.Initial_Conditions=IC;
DF.Constants=C;
DF.Construct=CSTR_LV;
DF.Parameter_library=P_est;
DF.iter_num=iter_num;

end


