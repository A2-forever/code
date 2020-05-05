function[M]=Normalize(M)
    n=size(M);
    n=n(1);
    for i=1:n
        s=M(i,1:n).*M(i,1:n);
        temp_Normalization=sum(s);
        temp_Normalization=sqrt(temp_Normalization);
        M(i,1:n)=M(i,1:n)/temp_Normalization;
        if M(i,1)<-eps
            M(i,1:n)=-M(i,1:n);
        end
    end
end