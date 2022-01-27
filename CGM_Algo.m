% Conjugate Gradient Method
p = [1;1;1]; % Guess parameters
N = size(p, 1);
I = 200;
Y = fModelCorrect;          % measurement vector
delp = 10^(-5);
J = zeros(I, N);            % sensitivity matrix 
d = 0;
k = 0;
while true
    k = k + 1
    T = fModelGuess(p);
    S = sum((Y-T).^2);
    e = 5;
 if abs(S) < e
     break;
 else
 for i = 1:I                % Constructing sensitivity matrix
     for j = 1:N
         pInc = p;
         pInc(j) = pInc(j) + delp*pInc(j);
         Tdelpj = fModelGuess(pInc);
         Tpj = fModelGuess(p);
         J(i,j) = (Tdelpj(i) - Tpj(i))/(delp*p(j));
     end
 end
 
 grad_S = -2*(J'*(Y-T)); 
     
 if k == 1
     gamma = 0;
     d = grad_S;
 else
 gamma = norm(grad_S)/norm(grad_prev);
 d = grad_S + gamma*d;
 end
 beta = ((J*d)'*(T-Y))/((J*d)'*(J*d));
 p = p - beta*d;
 grad_prev = grad_S ;
     
end   
end
