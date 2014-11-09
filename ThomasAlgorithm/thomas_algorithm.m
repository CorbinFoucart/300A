 % Thomas Algorithm
 % Corbin Foucart
 %
 %
 % Given: A,d (A is tridiagonal)
 % We look to solve the system Ax = d. 
 % The system has form:
 %
 % [ b1 c1 0   ...            0  ] [ x1 ]    [ d1 ]
 % [ a2 b2 c2 0   ...            ] [ x2 ]    [ d2 ]
 % [ 0  a3 b3 c3 0  ...       .  ] [ .  ]    [ d3 ]
 % [ .        .               .  ] [ .  ] =  [  . ]
 % [ .           .               ] [ .  ]    [  . ]
 % [ .              .         0  ] [    ]    [  . ]
 % [                       c_n-1 ] [    ]    [    ]
 % [ 0    ...        0  an   bn  ] [ xn ]    [ dn ]
 %
 % We will refer to this nomenclature in the comments.
 %
 % To this end, we will construct a c' vector which will contain the upper
 % diagonal entries, as well as a d' vector, which is the modified b
 % vector
 %
 % We then solve the system by back-substitution and return the solution.
 %
 % Here A,d are a storage device that we will refer to for values but not
 % explicitly change over the course of the algorithm. In constructing the
 % solution, we know that the lower diagonal entry will vanish and the main
 % diagonal entry will become one. No point wasting computations editing
 % the values of these entries in the matrix A.
 %
 % The recursive formulas:
 %          [ c_i/b_i                                 for i = 1
 %   c'_i = [      
 %          [ c_i / (b_i - a_i * c'_i-1)              for i = 2:n-1
 %
 %          [ d_i/b_i                                 for i = 1
 %   d'_i = [
 %          [(d_i - a_i*d'_i-1) / (b_i - a_i * c'_i)  for i = 2:n
 %
 %          [ d'_i                                    for i = n
 %   x_i =  [
 %          [ d'_i - c'_i*x_i-1                       for i = n-1 ... 1

 function rx = thomas_algorithm(A, d)
 
 % get the size of A
 n = length(A(:,1));
 
 % initialize the vectors to hold the solutions
 % c' contains n - 1 entries, d' contains n
 c_prime = zeros(n - 1, 1);
 d_prime = zeros(n , 1);
 
 x = zeros(n, 1);
 
 % n step process to form upper triangular matrix
 for i = 1:n
     % if the first step, need only to divide through by b1
     if (i == 1)
         %c'_1 = c1/b1
         c_prime(i) = A(1, 2)/A(1, 1);
         d_prime = d(1)/A(1,1);
     % steps 2 -> n    
     else
         % select the correct values out of the stored
         b_i = A(i, i);
         a_i = A(i, i - 1);
         d_i = d(i);
         
         % compute the next value of c', d' according to the 
         % recursive formula defined above, which actually
         % amounts to zeroing out the a value of the row
         % and dividing to make the b value 1.
         % note that c only runs through n - 1
         if (i ~= n)
             c_i = A(i, i + 1);
             c_prime(i) = c_i/(b_i - a_i * c_prime(i - 1));
         end
         
         d_prime(i) = (d_i - a_i*d_prime(i - 1))...
             /(b_i - a_i*c_prime(i - 1));
     end
 end
 
 % back substitute for the solution
 for i = 0 : n - 1
     if (i == 0)
         % the last entry is already solved
         x(n) = d_prime(n);
     else
         % back substitute working backwards; each successive 
         % iteration requires the value in the solution vector
         % that was computed immediately before it.
         x(n - i) = d_prime(n - i)- c_prime(n - i)*x(n - i + 1);
     end
 
  % return our solution
  rx = x;
 
 end
 