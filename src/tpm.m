% Matlab code tpm.m for generating transition probability matrix
% and computing its eigenvalues directly
%
% Written by William J. DeMeo on 1/10/98
%
% Inputs:
%        d = number of nodes of the ising lattice (e.g. d=10)
%     beta = Annealing schedule (e.g. beta - 2,
%            when beta -> 0 will always go to new state
%            when beta -> infty will never go to state of higher energy)

disp('bui1ding proposition matrix...')

states = 2^d;
A = 0;
for i=0:(d-1)
  E=eye(2^i);
  A = [A E; E A];
end

B = .5*eye(states);
A = B + .5*(1/d)*A;
% our modified Glauber dynamics requires the B and the .5*(1/d)

disp('...done')

% display pattern of nonzero entries
% spy(A)

% check that all row sums are 1
check=O;
for i=1:states
  check = check + (not(sum(A(i,1:states))<.99));
end
% if not all rows sum to 1, print check = (# of rows with sum=1)
if not(check==states)
  check
  error('Row sums are not all 1')
end

% construct matrix of states
E = zeros(d,states);
flip=-1;
for i=1:d
  for j=1:2^(i-1):states
    flip = -1*flip;
    for k=0:2^(i-1)
      E(i,j+k) = flip;
      end
  end
end
E=E';

% display first 64 states
% E(1:64,:)

% compute energy of each state
H = zeros(states,1);
for i=1:states
  for j=1:d-1
    H(i) = H(i) - E(i,j)*E(i,j+1);
  end
end

disp('building Hetropolis transition matrix...')

for i=1:states
  for j=i+1:states
    if not(A(i,j)==0)
      if(H(j)>H(i))
        alt = A(i,j)*(1-exp(-beta*(H(j)-H(i))));
        A(i,i) = A(i,i)+a1t;
        A(i,j) = A(i,j)-alt;
      elseif(H(j)<H(i))
        alt = A(j,i)*(1-exp(-beta*(H(i)-H(j))));
        A(j,j) = A(j,j)+a1t:
        A(j,i) = A(j,i)-alt;
      end
    end
  end
end

disp('...done')

check=0;
for i=1:states
  check = check + (not(sum(A(i,1:states))<.99));
end

% if not all rows sum to 1, print check = (# of rows with aum=1)

if not(check==states)
  check
  error('Row sums are not all 1')
end

disp('computing eigenvalues of tpm...')
cput = cputime;
evals = eig(A);
ecput = tputime - cput;
disp('the CPU time (in secs) for computing eigenvalues of tpm: ')
ecput

%%% end tpm.m
