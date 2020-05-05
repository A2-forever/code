eps = 1e-6;
fa = "¦Õ";
syms x;

data=importdata('Input_conjugate.txt');
n=data(1,1);
m=data(1,2);
ne=data(2,1);

l=size(data);
con_M=zeros(n);

for i=3:l(1)
    con_M(data(i,1), data(i,2)) = 1;
    con_M(data(i,2), data(i,1)) = 1;
end

E = eye(n);
for i=1:(m-1)
    con_M(i, i + 1) = 1;
    con_M(i + 1, i) = 1;
end

p = con_M + E * x;
d = det(p);
solution = solve(d, x);
solution = double(solution);
n_solution = size(solution);

k = 1;
store_c = zeros(n);
energy = -1*ones(1,n);
for i=1:n_solution
    
    if i > 1
        if abs(solution(i) - solution(i - 1)) < eps
            continue
        end
    end

    M = con_M + E * solution(i);
    r = rank(M);
    n_zero = n - r;
    k=k+n_zero;
    energy(k-n_zero:k-1) = solution(i);
    
    record=-1*ones(n);
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

    flag = 1;
    for I=k-n_zero:k-1
        for J=1:n
            if J == record(flag)
                store_c(I,J) = 1;
            elseif ismember(J,record)~=1
                store_c(I,J) = -M(J, record(flag)) / M(J,J);
            end
        end
        flag=flag+1;
    end
    
end
 
    
for i =1:n - 1
    if abs(energy(i) - energy(i + 1)) < eps
            store_c(i,1:n) = store_c(i,1:n) - store_c(i+1,1:n);
            store_c(i+1, 1:n) = store_c(i, 1:n) + 2*store_c(i+1, 1:n);
    end
end

store_c=Normalize(store_c);

t1 = ne / 2;
if abs(energy(t1) - energy(t1 + 1)) > eps
    t2 = t1;
else
    t2 = t1 + 1;
end

e_density = zeros(1,n);
e_bond = zeros(n);

for i=1:n
    s=store_c(1:t1-1,i) .* store_c(1:t1-1,i);
    e_density(i) = 2 * sum(s(:));
    e_density(i) =e_density(i)+ store_c(t1, i) * store_c(t1, i) + store_c(t2, i) * store_c(t2, i);
end

for i=1:n
    for j=i:n
        if con_M(i, j) > eps
            s=store_c(1:t1-1, i) .* store_c(1:t1-1, j);
            e_bond(i, j)=2 *sum(s(:));
            e_bond(i, j) =e_bond(i, j)+ store_c(t1, i) * store_c(t1, j) + store_c(t2, i) * store_c(t2, j);
        end
    end
end

e_free = ones(1,n) * sqrt(3);
for i=1:n
    for j=1:n
        if con_M(i,j) > eps
            e_free(i) = e_free(i) - e_bond(i,j)- e_bond(j,i);
        end
    end
end

K=[-energy',store_c];
xlswrite('Output_conjugate.xls',K);
fid = fopen('Output_conjugate.txt','w');
fprintf(fid,"Molecular Energy Level and Orbital: ");

for i=1:n
    fprintf(fid,"\n");
    k=-energy(i);
    k=round(k,4);
    k=num2str(k);
    if -energy(i) < -eps
        fprintf(fid,"Energy: ¦Á"+k+"¦Â £º");
    else
        fprintf(fid,"Energy: ¦Á+"+k+"¦Â £º");
    end
    
    for j=1:n
        k=store_c(i,j);
        k=round(k,4);
        k=num2str(k);
        if store_c(i,j) > eps &&j ~= 1 
            fprintf(fid,"+"+k+fa+num2str(j)+" ");
        else
            fprintf(fid,k+fa+num2str(j)+" ");
        end
    end
end

fprintf(fid,"\n");
fprintf(fid,"\n");

fprintf(fid,"Charge Density: ");
fprintf(fid,"\n");
for i=1:n
    k=e_density(i);
    k=round(k,4);
    k=num2str(k);
    fprintf(fid,"Carbon"+num2str(i)+" : "+k);
    fprintf(fid,"\n");
end
fprintf(fid,"\n");

fprintf(fid,"Bond Level: ");
fprintf(fid,"\n");
for i=1:n
    for j=i+1:n
        if con_M(i,j) > eps
            k=e_bond(i,j);
            k=round(k,4);
            k=num2str(k);
            fprintf(fid,"Carbon " +num2str(i)+ " and Carbon "+ num2str(j)+" is: "+k+"    ");
        end
    end
    fprintf(fid,"\n");
end

fprintf(fid,"Free Valence: ");
fprintf(fid,"\n");
for i=1:n
    k=e_free(i);
    k=round(k,4);
    k=num2str(k);
    fprintf(fid,"Carbon "+num2str(i)+" : "+k);
    fprintf(fid,"\n");
end
fprintf(fid,"\n");

fclose(fid);