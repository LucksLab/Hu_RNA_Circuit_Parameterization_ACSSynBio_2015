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
% Estimation_fmincon.m
% A function that estimates identifiable parameters by minimizing errors to experimental data 
% ------------------------------------------------------------------------------------- %

function [P_Estimated,P_all,GM_opt,fval,exitflag]=Estimation_fmincon(Exp_Data_DF,DF,tspan,pset,i)


    CSTR_LV=DF.Construct;
    P_lib=DF.Parameter_library;
    NParameters=DF.Num_Parameters;
    iter_num=DF.iter_num;
    P=P_lib(i,:); %Columns are the different parameters, rows are different sets
        
    function [Error] = fminconfunc( y )
             
            timesize=Exp_Data_DF.timestep; % get experimental data 
        
            for k=1:length(pset)  %update P with identifiable parameters
                P(pset(k))=y(k);
            end
            
            [ t, x ] = Double_Represion_Call_ODE(DF,tspan,P,CSTR_LV);
            GM_est=transpose(x(:,8)); %Basing estimate error on GFP prediction
            Error=0; 
           
            for j=1:(timesize)
               Error=Error+(GM_est(j)-Exp_Data_DF.avg(j)).^2; %Sum of Error 
            end
            
            
    end

    
    y0=P(pset(1));
    
    [ A,b ] = fmincon_constrains(pset,P);
    for m=2:length(pset)
        y0=[y0,P(pset(m))];
    end

    lb=0.33*y0; 
    ub=3.33*y0;

    
    [y,fval,exitflag] = fmincon(@fminconfunc,y0,A,b,[],[],lb,ub);
    P_Estimated=y;

    P_all=P;
    for i=m:(length(pset))
        P_all((pset(m)))=P_Estimated(m);
    end
    
    P=P_all;
    [ t, x ] = Double_Represion_Call_ODE(DF,tspan,P,CSTR_LV);
    GM_opt=transpose(x(:,8));    
    
end


