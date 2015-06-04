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
% Opt_datafile.m
% Only the first construct would read data_file.m, all the following
% constructs read this file and is updated from the previous construct.
% ------------------------------------------------------------------------------------- % 


function [DF] = Opt_datafile(DF,pset,P_Estimated_mtr,CSTR_LV)

P_i=DF.Initial_Parameters; %get initial guesses from last round
P_lib=DF.Parameter_library;


for i=1:(length(pset))
    P_lib(:,pset(i))= P_Estimated_mtr(:,i);   %update the paramter library
end
    
if CSTR_LV==3

    P_lib(:,2)=8*P_lib(:,3); 
    P_i(2)=8*P_i(3);
    
    DF.P5_lib_saver=P_lib(:,5);
    DF.P5_i_saver=P_i(5);
    P_lib(:,5)=P_lib(:,4);
    P_i(5)=P_i(4);
    
    DF.P8_lib_saver=P_lib(:,8);
    DF.P8_i_saver=P_i(8);
    P_lib(:,8)=P_lib(:,7);
    P_i(8)=P_i(7);
    
    DF.P13_lib_saver=P_lib(:,13);
    DF.P13_i_saver=P_i(13);
    P_lib(:,13)=P_lib(:,12);
    P_i(13)=P_i(12);

    pset=[4;7;12]; %Identifiable param at LV2+LV3
elseif CSTR_LV==2;
    
    P_lib(:,2)=8*P_lib(:,3); 
    P_i(2)=8*P_i(3);
    pset=[5;8;13;15]; %Identifiable param at LV1+LV2
 

elseif CSTR_LV==4 
    
    P_lib(:,5)=DF.P5_lib_saver;
    P_i(5)=DF.P5_i_saver;
    
    P_lib(:,8)=DF.P8_lib_saver;
    P_i(8)=DF.P8_i_saver;
    
    P_lib(:,13)=DF.P13_lib_saver;
    P_i(13)=DF.P13_i_saver;
    
    P_lib(:,1)=28*P_lib(:,3);
    P_i(1)=28*P_i(3);
    pset=[14]; %Identifiable param at LV1+LV2+LV3
end

DF.Construct=CSTR_LV;   
DF.Initial_Parameters=P_i;
DF.Parameter_library=P_lib;
DF.pset=pset;

end

