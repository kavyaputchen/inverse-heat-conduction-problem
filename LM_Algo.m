% Levenberg-Marquardt Method 

p = [1;1;1]; % Guess parameters
N = size(p, 1);
I = 200;
Y = fModelCorrect;          % measurement vector
delp = 10^(-5);
J = zeros(I, N);            % sensitivity matrix 
e1 = 3;
e2 = 1.5;
e3 = 3;
mu = 0.0001;
plot(1:200, Y);

count = 0;
while true
 count = count + 1
 T = fModelGuess(p);
 S = sum((Y-T).^2);         % Squared Error 

 
 for i = 1:I                % Constructing sensitivity matrix
     for j = 1:N
         pInc = p;
         pInc(j) = pInc(j) + delp*pInc(j);
         Tdelpj = fModelGuess(pInc);
         Tpj = fModelGuess(p);
         J(i,j) = (Tdelpj(i) - Tpj(i))/(delp*p(j));
     end
 end

 omega = diag(diag(J'*J));   

 pupdt = p;

 while true
     
     pupdt = pupdt + (pinv(J'*J + mu*omega))*(J')*(Y - fModelGuess(pupdt));   
     T = fModelGuess(pupdt);         % updated T
     Supdt = sum((Y-T).^2);
     if Supdt >= S         
       mu = 10*mu;
     else
       mu = 0.1*mu;  
       break;
     end
 end

 if (Supdt<e1 || norm(J'*(Y - fModelGuess(p)))<e2 || norm(pupdt-p)<e3)
     pEstimate = pupdt
     break;
 else
     p = pupdt;
 end
 
end