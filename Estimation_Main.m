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
% Estimation_Main.m
% The main script that estimates parameters
% ------------------------------------------------------------------------------------- % 
clear all
close all

CSTR_LV=1;

t_i=0;
t_f= 6000;   
t_inc = 300;  
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;

randn('state',sum(100*clock));
rand('state',sum(100*clock));

while CSTR_LV<=4
    
    if CSTR_LV==4   
        [DF] = Opt_datafile(DF,pset,P_Estimated_mtr,CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV );
    elseif CSTR_LV==3
        [DF] = Opt_datafile(DF,pset,P_Estimated_mtr,CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV );
    elseif CSTR_LV==2
        [DF] = Opt_datafile(DF,pset,P_Estimated_mtr,CSTR_LV);
        Exp_Data_DF= Exp_Data( CSTR_LV);
    elseif CSTR_LV==1
        [DF] = Datafile (CSTR_LV);
        Exp_Data_DF= Exp_Data(CSTR_LV);
        Data=Exp_Data_DF.Data;
        Exp_max=max(max(Data));
        Exp_min=min(min(Data));
        
    else disp ('Error,Construct level ranges from 1-4');
    end
    clearvars P_Estimated_mtr
   
    P=DF.Initial_Parameters;  

    [t,x] = Double_Represion_Call_ODE(DF,tspan,P,CSTR_LV); % The ODE function doesnt change with construct
    GM_t= transpose(x(:,8));

    for i=1:(nstep+1)
        GM_t_Norm(i)=(GM_t(i)-GM_t(1))./(Exp_max-Exp_min); %normalize Gm_t by experimental data 
    end

    figure(CSTR_LV)   %generate a profile with initial guess, pink line
    plot(t,GM_t_Norm,'m-','LineWidth',2); 
    xlabel('time'), ylabel('concentration')
    hold on

    GM_exp=Exp_Data_DF.Data;
    GM_exp_norm=GM_exp./(Exp_max);

    plot(t,GM_exp_norm,'b--'); %Experimental lines
    hold on

pset=DF.pset;
%------------------------------------------------------
%This part of the program Estimates the identifiable parameters using the experimental
%data 

    loopsz=DF.iter_num;
    GM_opt_Norm_sum=zeros(1,(nstep+1));
    for i =1:(loopsz)
        [P_Estimated,P_all,GM_opt,fval,exitflag]=Estimation_fmincon(Exp_Data_DF,DF,tspan,pset,i);
        P_Estimated_mtr(i,:)=P_Estimated;
        if CSTR_LV==4
            P_solution(i,:)=P_all; %Gathering final solution
        end
        
        for j=1:(nstep+1) %normalize GM_opt
            GM_opt_Norm(j)=(GM_opt(j)-GM_opt(1))./(Exp_max-Exp_min);       
        end
        GM_Mtx(i,:)=GM_opt_Norm;  
        GM_opt_Norm_sum=GM_opt_Norm_sum+GM_opt_Norm;
    end

    for j=1:(nstep+1)
        GM_std(j)=std(GM_Mtx(:,j));  % find standard deviation of the plot
    end 

    GM_MEAN=GM_opt_Norm_sum./(loopsz); 
    GM_UB=GM_MEAN+1.96*GM_std; %95% confidence interval
    GM_LB=GM_MEAN-1.96*GM_std;


    plot(t,GM_Mtx ,'g');%individual modeled lines
    hold on
    plot(t,GM_MEAN,'k','LineWidth',2);% mean modeled line
    hold on
    plot(t,GM_UB,'k--',t,GM_LB,'k--');% 95% confidence
    hold all
% --------------------------------------------------------
% part end

    if CSTR_LV==4
        dlmwrite('P_solution.txt',P_solution','delimiter','\t'); 
    end

    CSTR_LV = CSTR_LV+1; % move on to the next construct


 end


