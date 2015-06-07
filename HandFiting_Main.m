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
% HandFitting_Main.m
% This program plots the trajectories with given initial guessed parameters
% stored in the txt file "guess.txt". By comparing trajectories, you can generate a
% set of guesses manually.
% ------------------------------------------------------------------------------------- % 

clear all
close all

CSTR_LV=1;   %CSTR_LV==construct level
t_i=0;
t_f= 6000;   % run time 100min
t_inc = 300; % interval 5min
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;


while CSTR_LV<=4  
    if CSTR_LV==4
        [ DF ] = Datafile(CSTR_LV );
        Exp_Data_DF= Exp_Data( CSTR_LV );
    elseif CSTR_LV==3 
        DF = Datafile(CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV );
    elseif CSTR_LV==2
        DF =  Datafile(CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV );
    elseif CSTR_LV==1
        DF = Datafile( CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV );
        Data=Exp_Data_DF.Data;
        Exp_max=max(max(Data));
        Exp_min=min(min(Data));
    end
    P=DF.Initial_Parameters;
    [t,x] = Double_Represion_Call_ODE(DF,tspan,P,CSTR_LV);
    GM_opt=transpose(x(:,8));
    for j=1:(nstep+1)
        GM_opt_Norm(j)=(GM_opt(j)-GM_opt(1))./(Exp_max-Exp_min); 
    end
    
    if CSTR_LV==1
        color='g'; 
    elseif CSTR_LV==2
        color='b';
    elseif CSTR_LV==3
        color='m';
    elseif CSTR_LV==4
        color='r';
    end
    figure(CSTR_LV)   %generate a profile 
    plot(t,GM_opt_Norm ,color,'LineWidth',2);
    xlabel('time'), ylabel('concentration')
    hold all
    GM_exp=Exp_Data_DF.Data;
    GM_exp_norm=GM_exp./(Exp_max);
%     axis([0 6000 0 2])
    plot(t,GM_exp_norm,color,'LineStyle','--'); 
    hold on
    CSTR_LV=CSTR_LV+1;
end
%------------------------------------------------
   