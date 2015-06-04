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
% Data_Histogram.m
% Makes Histograms for all estimated parameters
% ------------------------------------------------------------------------------------- %
close all
CSTR_LV=4;

[DF] = Analysis_func(CSTR_LV); %Load iteration number from struct
iter_num=DF.iter_num;

P_est=zeros((iter_num),15);
filename=['P_solution.txt'];
file=fopen(filename);
P_est= fscanf(file, '%f',[(iter_num),15]); 



IC=[0           %A'2_0 
    0           %A2_0
    0           %A'1_0 
    0           %A1_0
    0           %M_0 
    0           %Mi_0
    0           %G_0 
    0];         %GM_0 




C=[0          % k_c,the crosstalk repression coefficient of A_2 on att-1
   0            % p_t,the probability of a attentuator auto-determinates itself
   1               % Level 1 order
   1] ;            % Level 2 order

DF.Num_Parameters=length(P_est(1,:));
DF.Initial_Conditions=IC;
DF.Constants=C;
DF.Construct=CSTR_LV;
DF.Parameter_library=P_est;
DF.iter_num=iter_num;


for i=1:15
    x=P_est(:,i);
    figure (i)
    hist(x,50);
    xlabel('value');
    ylabel('frequency');
    title_name=['Paremeter ',num2str(i)];
    title(title_name);
    hold all
end

