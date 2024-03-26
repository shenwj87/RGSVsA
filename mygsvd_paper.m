% data=unnamed; clearvars unnamed
clc;clear;
for i=1:1%size(data,1)

%% 构造矩阵 A低秩B低秩[A;B]满秩 rank(A)+rank(B)>n
m = 1500; p = 1500; n = 1500; ra = 800; rb = 800;
aa=sort(rand(n-ra,1)*1e-15,'ascend');
bb=sort(rand(n-rb,1)*1e-15,'descend');
b=[sort(1-aa,'descend');sort(rand(ra+rb-n,1),'descend');bb];
a=diag(sort(sqrt(1-b.^2),'ascend'));
a(1:n-ra,1:n-ra)=diag(aa);
a(rb+1:end,rb+1:end)=diag(sqrt(1-bb));
sigma1=diag(a);
sigma2=b;
Sigma1=a;
Sigma2=diag(b);

if m<n && p<n
    Sigma1=Sigma1(n-m+1:end,:);
    Sigma2=Sigma2(1:p,:);
    U=orth(randn(m)+randn(m)*1i);
    V=orth(randn(p)+randn(p)*1i);
elseif m>n && p<n
    Sigma2=Sigma2(1:p,:);
    U=orth(randn(m,n)+randn(m,n)*1i);
    V=orth(randn(p)+randn(p)*1i);
elseif m<n && p>n
    Sigma1=Sigma1(n-m+1:end,:);
    U=orth(randn(m)+randn(m)*11);
    V=orth(randn(p,n)+randn(p,n)*1i);
else
    U=orth(randn(m,n)+randn(m,n)*1i);
    V=orth(randn(p,n)+randn(p,n)*1i);
end

X=randn(n)+randn(n)*1i;
A=U*Sigma1*X;
B=V*Sigma2*X;
ra=rank(A);rb=rank(B);rab=rank([A;B]);

tol = 1e-3;
s1=sigma1;
s1(s1<tol)=[];
s2=sigma2;
s2(s2<tol)=[];

%% 投影方法
tic;
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
t(i,1) = toc;
qn(i,1)=norm(A-Q1*Q1'*A,'fro');
qn(i,2)=norm(B-Q2*Q2'*B,'fro');
e1(e1<tol)=[];
e2(e2<tol)=[];
no(i,1)=norm(e1-s1(1:size(e1,1)),'fro');
no(i,2)=norm(e2-s2(size(s2,1)-size(e2,1)+1:end),'fro');
clearvars Q1 Q2 QA QB r1 r2 L0 L1 e0 e1 e2
%% gsvd
tic;
[~,~,~,S,X]=gsvd(A,B);
t(i,2)=toc;
S=S(S~=0);%s=S(size(S):-1:1);%S降序
X=X(X~=0);%x=X(size(X):-1:1);%X升序
S(S<tol)=[];
X(X<tol)=[];
no(i,3)=norm(S-s1(1:size(S,1)),'fro');
no(i,4)=norm(X-s2(size(s2,1)-size(X,1)+1:end),'fro');
clearvars S X s x
%% gsvd0
tic;
[~,~,~,S0,X0]=gsvd(A,B,0);
t(i,3)=toc;
S0=S0(S0~=0);
X0=X0(X0~=0);
S0(S0<tol)=[];
X0(X0<tol)=[];
no(i,5)=norm(S0-s1(1:size(S0,1)),'fro');
no(i,6)=norm(X0-s2(size(s2,1)-size(X0,1)+1:end),'fro');
clearvars S0 X0 s0 x0
%% gsvd in SIAMX
tic;
[phi,psi,rsiamx]=siamx(A,B);
t(i,4)=toc;
phi=sort(diag((phi)),'ascend');
psi=sort(diag((psi)),'descend');
phi(phi<tol)=[];
psi(psi<tol)=[];
no(i,7)=norm(phi-s1(1:size(phi,1)),'fro');
no(i,8)=norm(psi-s2(size(s2,1)-size(psi,1)+1:end),'fro');
clearvars phi psi rsiamx
end
value=[m,p,n,t,qn,no];
% value=[data(1:i,:),rA',rB',t,qn,no,rAB'];

