% Matlab code ising.m
%
% Written by William J. DeMeo on 12/15/97
% last modified 2013.10.19
% Inputs:
%         d = number of nodes of the ising lattice (e.g. d=100)
%         n = number of iterations (e.g. n=1000)
%         beta = Annealing schedule (e.g. beta = 2,
%            when beta -> 0 will always go to nee state
%            when beta -> infty will never go to state of higher energy)
%
X = zeros(d,1);
Energy = zeros(n,1);
Ave=0;

% Initialize using hot start
for i=1:d,
  if rand < .5
    X(i) = -1;
  else
    X(i) = 1;
  end
end

% Initialize at state of high energy (not used)
if 0
  X(1:2:d-1)=1;
  X(2:2:d)=-1;
end

% Initialize energy
H=0;
for i=1:(d-1)
  H = H - X(i)*X(i+1);
end

% Initialize observable
phi = zeros(n+1,1);

for i=1:d
  if X(i) == 1
    phi(1)= phi(1)+1;
  end
end

for j=1:n
  site=0;
  while site==0
    site = round(rand*d);
  end
  % Compute new energy (with new X(site) = -X(site)):
  if site == 1
    Hnew = H + 2*X(1)*X(2);
  elseif site == d
    Hnew = H + 2*X(d)*X(d-1);
  else
    Hnew = H + 2*X(site)*(X(site-1) + X(site+1));
  end

  % Change to new state in two ways:
  % with probability 1/2, flip the switch
  if rand < .5
    if (Hnew <= H) | (rand < exp(-beta*(Hnew-H)))
      X(site) = -X(site);
      H = Hnew;
      % Compute new value of observable
      if X(site) == 1
        phi(j+1) = phi(j) + 1;
      else
        phi(j+1) = phi(j) -1;
      end
    else phi(j+1)=phi(j);
    end
  else phi(j+1)=phi(j);
  end
  Energy(j) = H;
end
plot(Energy)

% A second observable
phi2 = zeros(n+1,1);
for i=1:n+1
  phi2(i) = phi(i)*phi(i);
end

% Write first observable to phi.dat:
fid = fopen('phi.dat','w');
fprintf(fid,'%f\n',phi);
fclose(fid);

% Write second observable to phi2.dat:
fid = fopen('phi2.dat','w');
fprintf(fid,'%f\n',phi2);
fc1ose(fid);

% Compute a few Lanczos coefficients by new method
START=1001;
V1=zeros(2);
C1=zeros(2);
C2=zeros(2);
C3=zeros(2);

V1 = cov(phi(START:n+1),phi(START:n+1));
V1 = V1(1,1);
C1 = cov(phi(START:n),phi(START+1:n+1));
C1 = C1(1,2);
C2 = cov(phi(START:n-1),phi(START+2:n+1));
C2 = C2(1,2);
C3 = cov(phi(START:n-2),phi(START+3:n+1));
C3 = C3(1,2);

alpha1 = C1/V1;
beta1 = sqrt(C2/V1 -(C1/V1)^2);
alpha2 = (1/(beta1^2))*(C3/V1 - 2*alpha1*(C2/V1) + alpha1^3);

T1 = [alphal beta1; beta1 alpha2]
V1 = cov(phi2(START:n+1),phi2(START:n+1));
V1 = V1(1,1);
C1 = cov(phi2(START:n),phi2(START+1:n+1));
C1 = C1(1,2);
C2 = cov(phi2(START:n-1),phi2(START+2:n+1));
C2 = C2(1,2);
C3 = cov(phi2(START:n-2),phi2(START+3:n+1));
C3 = C3(1,2);

alpha1 = C1/V1;
beta1 = sqrt(C2/V1 - (C1/V1)^2);
alpha2 = (1/(beta1^2))#(C3/V1 - 2*alpha1*(C2/V1) + alpha1^3);

T2 = [alpha1 beta1; beta1 alpha2]

fid = fopen('Tmats.dat','v');
fprintf(fid,'%f\n',T1);
fprintf(fid,'%f\n',T2);
fclose(fid);

%%% end ising.m
