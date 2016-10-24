function [A, b] = convert_to_cspace(A, b)
  A = [A, zeros(size(A, 1), 1);
       zeros(2, size(A, 2) + 1)];
  A(end-1,3) = 1;
  A(end,3) = -1;
  b = [b; 0; 0];
end