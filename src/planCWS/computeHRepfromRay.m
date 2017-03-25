function [A,b] = computeHRepfromRay(ray,plane)
% compute the halfspace representation A*w<=b from the rays 
% @param plane  We use the plane plane(1:end-1)*w=plane(end) to determine a
% intersecting hyperplane with the cone convcone(ray)
dim = size(ray,1);
num_ray = size(ray,2);
if(any(size(plane) ~= [1,dim+1]))
  error('plane should be of size 1 x %d,dim+1');
end
ray_len = plane(end)./(plane(1:end-1)*ray);
if(any(isnan(ray_len)) || any(isinf(ray_len)))
  error('The plane does not intersect with the ray');
end
ray_vert = bsxfun(@times,ones(dim,1),ray_len).*ray;
plane_null = null(plane(1:end-1));
plane_pt = plane(1:end-1)\plane(end);
ray_plane_pt = plane_null\(ray_vert-bsxfun(@times,ones(1,num_ray),plane_pt));
K = convhulln(ray_plane_pt');
num_facet = size(K,1);
A = zeros(num_facet,dim);
b = zeros(num_facet,1);
valid_row = true(num_facet,1);
for i = 1:num_facet
  null_ray = null(ray(:,K(i,:))')';
  if(size(null_ray,1)~=1)
    valid_row(i) = false;
  else
    A(i,:) = null_ray;
    if(sum(A(i,:)*ray)>0)
      A(i,:) = -A(i,:);
    end
  end
end
A = A(valid_row,:);
b = b(valid_row,:);
M = A*A';
M_triu = M.*(~tril(ones(size(M))));
valid_row = ~any(M_triu>1-1e-5,1);
A = A(valid_row,:);
b = b(valid_row,:);
end