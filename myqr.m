function Q1 = myqr(A,tol,b)
    n = size(A,2);
    Y = A*randn(n,b);
    [Q,~] = qr(Y,0);
    Q1 = [];
    ny1 = norm(Y,"fro");
    ny = ny1;
    while ny>tol
        Y = A*randn(n,b);
        Q1 = [Q1,Q];
        Y = Y-Q1*(Q1'*Y);
        ny2= norm(Y,"fro");
        if ny2>=ny1
            break
        elseif size(Q1,2)==n
                break
        else
            [Q,~] = qr(Y,0);
            ny = ny2;
        end
        
    end
end