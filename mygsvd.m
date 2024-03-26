function [e1,e2] = mygsvd(A,B)
Q1 = myqr(A,1e-4,100);
Q2 = myqr(B,1e-4,100);
r1 = size(Q1,2);
QA=Q1'*A;
QB=Q2'*B;
[L0,~]=qr([QA;QB],0);
L1=L0(1:r1,:);
e0=eig(L1*L1');
e1=sort(sqrt(e0),'ascend');
e2=sort(sqrt(1-e0),'descend');
end