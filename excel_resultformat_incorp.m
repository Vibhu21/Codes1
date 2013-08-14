% © Vibhu Prasad, University of Zurich 
% To extract and analyse the the raw data from HTM image processing analysis. 

%Last edited: 4 May 2013 1917hrs

function [excel_resultformat_incorp] = excel_resultformat_incorp(num1,DataFolder)

% fluoindex and No of select nuclei from the HTM script

% extracting numerical values from the cell array 
idx = cellfun(@ischar,num1);
num2= zeros(size(num1));
num2(~idx) = cell2mat(num1(~idx));
num2(idx) = str2double(regexprep(num1(idx),'\D',''));

   
[num,txt]=xlsread([DataFolder,'platedesign.xls']);
x=unique(txt);

%Following two for loops determins the mean average nuclear intensity of
%the non-infected smaples from the excel sheet. It is important to put the
%determiner for non-infected samples as 'NI'. For exepriments where only
%single non-infected are present can also be used.
D=[];
C=[];
for y=1:1
    [row,col] = find(ismember(txt,'NI')==1);
        for j=1:length(col)
            A=col(j,:);
            B=row(j,:);
            NI_value(j,1)=num2(B+1,A+1);
            
        end
        
    %Check whether instances for non-infected sampels appeared once or more
    %times. If once, then there is no need to mean it. Similarly, standard
    %deviation is only calculated if the instance occured more than once.
            if length(col)==1
                NI_meanvalue(y,:)=(NI_value);            
            else
                NI_meanvalue(y,:)=mean(NI_value);
            end
end


% Following 2 for loops run through the whole length of unique cases (
% or unique siRNA targets and determine their Mean Nuclear Intensity
% substracted with mean nuclear intensity of non-infecteds. It is improtant
% to put the same detemriners for duplicate siRNAs and not siRNA1 & 2 and
% so on.
for y=1:length(x)
    C=0;
    z=x(y,:);
    [row,col] = find(ismember(txt,z)==1);
        for j=1:length(col)
            A=col(j,:);
            B=row(j,:);
            C(j,1)=num2(B+1,A+1)-NI_meanvalue;
            C(j,2)=num2(B+10,A+1);
            C(j,3)=num2(B+19,A+1);
            A=0;
            B=0;
        end
    %Check whether instances for particular siRNA appeared once or more
    %times. If once, then there is no need to mean it. Similarly, standard
    %deviation is only calculated if the instance occured more than once.
    
    if length(col)==1
        D(y,:)=C;
        E(y,:)=std(C,0,1);
    else
        D(y,:)=mean(C);
        E(y,:)=std(C,0,1);
    end
end


%normalize the data here
colnames_nucnorm={};
colnames_stdnorm={};
colname_final={};
X_norm=[];
D_norm=[];
E_norm=[];
num_nrmlze=input('Enter the number of normalizations for the experiment : ');

if num_nrmlze>0
    for v=1:num_nrmlze
        nrmlze=cell(10,1);
        nrmlze{v}=input('Enter the normalization identifier :','s');
        [row,col] = find(ismember(x,nrmlze{v})==1); %This line finds where in the
        %column of unique siRNA names does the normalzining siRNA comes
            if isempty(row)==1
            error('Wrong normalizing indentifier entered, Please recheck (note that indentifiers are case-sensitive)');
            end
    
        nrmlzng_AvNuc(v)=D(row,1); %extracts the mean nuc intensity of normalizing siRNA

        %This loops normalizes every mean nuc intensity of individual targets
        %with the normalizing Nuc Intensity
    
            for y=1:length(x)
                D_norm(y,0+v)=(100/nrmlzng_AvNuc(v))*D(y,1);
                
                E_norm(y,0+v)=(E(y,1)/D(y,1))*D_norm(y,v);
            end
            
            colnames_nucnorm{1,v}=cellstr(['Norm_AvNuc ' num2str(0+v)]);
            colnames_stdnorm{1,v}=cellstr(['Norm_Stdev ' num2str(0+v)]);            
    end
    colnames={'Identifiers','Mean Nuc Int','Mean Fluo Index','Mean Cell Number', 'stdev nuc int' 'stdev fluo index','stdev cell number'};
    colname_final=horzcat(colnames{:},(horzcat(colnames_nucnorm{:},colnames_stdnorm{:})));
    X_norm=horzcat(horzcat(num2cell(D_norm)),num2cell(E_norm));
end
            

X=horzcat(horzcat(x,num2cell(D)),num2cell(E));
X_final=horzcat(X,X_norm);
excel_resultformat_incorp=vertcat(colname_final,X_final);
disp('Analysis finished successfully');
