function M=Gauss(M)
    size_M=size(M);
    n=size_M(1);
    r=rank(M);
    cout=r;
    n_z=0;
    for i1=1:n
        if abs(M(i1,i1))<=eps && i1>cout
            n_z=n_z+1;
            record(n_z)=i1;
            continue
        elseif abs(M(i1,i1))<=eps && i1<=cout
            flag=0;
            for position=i1+1:n
                if abs(M(position,i1))>eps
                    temp=M(position,1:n);
                    M(position,1:n)=M(i1,1:n);
                    M(i1,1:n)=temp;
                    flag=1;
                    break
                end
            end
            if flag==0
                cout=cout+1;
                n_z=n_z+1;
                record(n_z)=i1;
                continue
            end
        end
        
        for i2=1:n
            if abs(M(i2,i1))<=eps || i2==i1
                continue
            end
            co=M(i2,i1)/M(i1,i1);
            M(i2,1:n)=M(i2,1:n)-co*M(i1,1:n);
        end
        
    end
end