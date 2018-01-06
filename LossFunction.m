function [score]=LossFunction(para_)
%% LossFunction(): aim at minimal value
    para=zeros(length(para_));
    for i=1:length(para)
        para(i)=para_{i};
    end
%% The loss function
    dot_=(para(1)-1)^2+(para(2)-2)^2+(para(3)-3)^2+(para(4)-4)^2;
    score=sum(dot_);
end