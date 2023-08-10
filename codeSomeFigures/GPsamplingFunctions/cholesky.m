function L=cholesky (a)
% Descomposición Cholesky
% a: Matriz cuadrada
[n]=size(a,1);
[L]=zeros(n,n);
for k = 1 : n
   for i = 1 : k-1
     soma = 0;
     for j = 1 : i-1
       soma = soma + a(i,j) * a(k,j);
     end
     a(k, i) = (a(k, i) - soma)/a(i, i);
   end
   soma = 0;
   for j = 1 : k -1
     soma = soma + a(k,j)^2;
   end
a(k,k) = (a(k,k) - soma)^.5;
end
L=tril(a);

