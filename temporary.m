
clear all;
clc;
Delta=[1 0 0 0 0;1 1 0 0 0;1 1 1 0 0;1 1 1 1 0;1 1 1 1 1];



Sbin7=[1 0 0 0 0;
       1 1 0 0 0;
       0 1 0 0 0;
       0 1 0 0 0;
       0 1 0 0 1];

Sbin7=[1 0 0 0 0;
       1 1 0 0 0;
       0 1 0 0 0;
       0 1 0 0 0;
       0 1 0 0 1];



Sbin=Sbin7;
QI=test_QI(Sbin,Delta)

if(QI==1)
    disp('STOOOOOOOOOOOOOOOOOOOOOOOOOP')
    disp('STOOOOOOOOOOOOOOOOOOOOOOOOOP')
    disp('STOOOOOOOOOOOOOOOOOOOOOOOOOP')
    disp('STOOOOOOOOOOOOOOOOOOOOOOOOOP')
    disp('STOOOOOOOOOOOOOOOOOOOOOOOOOP')
    disp('STOOOOOOOOOOOOOOOOOOOOOOOOOP')
    disp('STOOOOOOOOOOOOOOOOOOOOOOOOOP')
end