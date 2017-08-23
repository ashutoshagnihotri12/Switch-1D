function [ f ] = newKorenLimiter(theta)
A = 2*theta;
B = (ones(1,length(theta)) + 2*theta)/3;
C = 2*ones(1,length(theta));
D = min([A;B;C]);
f = max([zeros(1,length(theta));D]);
end

