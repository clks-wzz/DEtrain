function [ err ] = EnergyFunction( varargin)
para=varargin;
%AUCpara=w1; Fpara=w2; MAEpara=w3
global iteration_num;
global besterr;

str_out='';
str_out=strcat(str_out,['itr  ' num2str(iteration_num) ':']);
for i=1:length(varargin);
    str_out=strcat(str_out,['  ' num2str(varargin{i})]);
end
disp(str_out);
%disp(['itr  ' num2str(iteration_num) ':  ' num2str(w1) '  ' num2str(w2) '  ' num2str(w3)  '  ' num2str(w4)]);
result=LossFunction(para);

iteration_num=iteration_num+1;
err=result;
disp(['Accuracy is ' num2str(err)]);
if(err<besterr)
    besterr=err;
    save ./bestpara.mat para err;
end

end



