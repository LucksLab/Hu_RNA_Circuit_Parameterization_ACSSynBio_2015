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
% Find_Primary_Guess_Main.m
% This program uses fmincon function to estimate parameters based on
% provided experimental trajectories. With each trial, it finds the best
% estimated parameter set and updates the best guesses.
% ----CAUTION----
% Running this script would change the primary paramters set.
% To avoid this, comment out the very last line.
% ----CAUTION----
% ------------------------------------------------------------------------------------- % 

clear all
close all
CSTR_LV=1;

t_i=0;
t_f= 6000;   
t_inc = 300;  
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;

[DF] = Analysis_func(CSTR_LV); 
loop_num=20;
Sum_error=zeros(loop_num,1); 

while CSTR_LV<=4
    if CSTR_LV==4   
        [DF] = Analysis_func(CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV );
    elseif CSTR_LV==3
        [DF] = Analysis_func(CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV );
    elseif CSTR_LV==2
        [DF] = Analysis_func(CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV );
    elseif CSTR_LV==1
        [DF] = Analysis_func(CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV );
        Data=Exp_Data_DF.Data; %all trajactories are normalized with the max level of LV1
        Exp_max=max(max(Data));
        Exp_min=min(min(Data));
    end
%load the Data file
    
    Norm_Exp_Avg=Exp_Data_DF.avg./Exp_max;
%--------------------------------------------------
    
    iter_num=DF.iter_num;
    loopsz=iter_num;
    GM_opt_Norm_sum=zeros(1,(nstep+1));
    GM_error=zeros(loop_num,1);
    
    for i =1:loop_num
        P_lib=DF.Parameter_library;
        NParameters=DF.Num_Parameters;
        P=P_lib(i,:);
        P_saved(i,:)=P; 
        [t,x] = Double_Represion_Call_ODE(DF,tspan,P,CSTR_LV);
        GM_opt=transpose(x(:,8));
        for j=1:(nstep+1)
            GM_opt_Norm(j)=(GM_opt(j)-GM_opt(1))./(Exp_max-Exp_min);  
            GM_error(i)= GM_error(i)+abs(GM_opt(j)- Exp_Data_DF.avg(j));
        end
    
        GM_Mtx(i,:)=GM_opt_Norm;  
        GM_opt_Norm_sum=GM_opt_Norm_sum+GM_opt_Norm;      
    end

    for j=1:(nstep+1)
        GM_std(j)=std(GM_Mtx(:,j));  
    end 


    GM_MEAN=GM_opt_Norm_sum./loop_num; 
    GM_UB=GM_MEAN+GM_std;
    GM_LB=GM_MEAN-GM_std;

    figure(CSTR_LV)   %generate a profile 
    plot(t,GM_Mtx ,'g');
    hold on
    plot(t,GM_MEAN,'k','LineWidth',3);
    plot(t,GM_UB,'k--',t,GM_LB,'k--');
    xlabel('time'), ylabel('concentration')
    hold all

    GM_exp=Exp_Data_DF.Data;
    GM_exp_norm=GM_exp./(Exp_max);

    plot(t,GM_exp_norm,'b--'); 
    hold on
    
    Sum_error=Sum_error+GM_error; 

    CSTR_LV=CSTR_LV+1;
end

[C,Index]=min(Sum_error);
Best_fit=P_saved(Index,:);
% dlmwrite('guess.txt',Best_fit','delimiter','\t');

