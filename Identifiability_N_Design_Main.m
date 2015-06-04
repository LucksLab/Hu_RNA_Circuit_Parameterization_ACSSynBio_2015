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
% Identifiability_N_Design_Main.m
% This program uses the initial guess parameter set to compute sensitivity
% matrices for each proposed construct. Then use "Identifiability.m" function
% to compute the identifiable paramters of each proposed construct
% ------------------------------------------------------------------------------------- % 

clear all
close all

CSTR_LV=input('Input CSTR_LV:');
Prev_ID=input('previously identified :');%(Include ALL previsouly identified parameter indices (single integers), input the array in a [], separated by ;)
Prev_ID=sort(Prev_ID);


t_i=0;
t_f= 6000;   
t_inc = 30; %Using finer time step for sensitivity matrix calculation, though not experimentally sampling at this fine of an interval  
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;

if CSTR_LV==4     
    DF = Datafile(CSTR_LV);
    Exp_Data_DF= Exp_Data( CSTR_LV );
elseif CSTR_LV==3
    DF = Datafile(CSTR_LV);
    Exp_Data_DF= Exp_Data( CSTR_LV );
elseif CSTR_LV==2
    DF = Datafile(CSTR_LV);
    Exp_Data_DF= Exp_Data( CSTR_LV );
elseif CSTR_LV==1
    DF = Datafile(CSTR_LV);
    Exp_Data_DF= Exp_Data( CSTR_LV );
        
else disp ('Error,Construct level ranges from 1-4');
end

P=DF.Initial_Parameters;  %Obtain the initial guess vector from struct

[t,x] = Double_Represion_Call_ODE(DF,tspan,P,CSTR_LV); % The ODE function doesnt change along with construct
GM_t= transpose(x(:,8));

[ SSM ] = Double_Repression_Analytical_SSM( DF,x,nstep,t_inc,Prev_ID);  %this function uses the analytical method to find SSM



[List,pset]=Identifiability (SSM);

if CSTR_LV==3 % return real psets for CSTR_LV3
    for i=1:size(pset)
        if pset(i)==5 
            pset(i)=4;
            SSM(:,[4,5])=SSM(:,[5,4]);
        elseif pset(i)==8
            pset(i)=7;
            SSM(:,[7,8])=SSM(:,[8,7]);
        elseif pset(i)==13 
            pset(i)=12;
            SSM(:,[12,13])=SSM(:,[13,12]);
        end
    end 
     
end

%--------SSM heat map
A=abs(SSM');
y=1:15;
x=0:5:100;
figure
imagesc(x,y,A)
caxis([0 0.5]) 
xlabel('Parameters'), ylabel('Time(min)')
colorbar
fig=figure(1);
filename=['SSM',num2str(CSTR_LV)];
print (fig,filename,'-dpng');
%--------

disp('Identifiable Parameters are:');
disp(pset);
