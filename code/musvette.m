function [ output_args ] = musvette( input_args )
%UNTÝTLED2 Summary of this function goes here
%   Detailed explanation goes here

% load('Omer-Config0-Path1.mat','SB');
% L=SB;
% L{2}
% load('Omer-Config0-Path2.mat','SB');
% O=SB;
% O{2}

% a=1:5;
% plot(1.5:0.5:4.5,[1 1 1 1 1 1 1],1.5:0.5:4.5,[1 1 1 1 1 1 1]*2)
% axis([1 5 -2 5])
% legend('mal','Location','North','FontSize',5)


poly_degrees = [4,8];
all_errors=[ 1 3 1 3 1 3 ];
a=cell(1,1)
a{1}{1}=[11245 5635]
a{2}{1}=[1145 5635]
a{1}{2}=[245 5635]
a{2}{2}=[45 5635]
figure
plot([0:2:6],all_errors([1:2:6,6]),'-+')
ac=gca;
ac.XTick = 0:2:6-1;
save('ada.mat','poly_degrees')
save('all_errors.mat','all_errors')
save('aaaa.mat','a')

end

