function mat = crossSkewSymMat(v)
% return the matrix mat such that cross(v,u) = mat*u for any u
if(any(size(v) ~= [3,1]))
  error('v should be 3 x 1');
end
mat = [0 -v(3) v(2);...
       v(3) 0 -v(1);...
       -v(2) v(1) 0];
end