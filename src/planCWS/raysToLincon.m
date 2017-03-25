function [Ain,bin] = raysToLincon(rays)
dim = size(rays,1);
vert = [zeros(dim,1) rays];
K = convhulln(vert',{'QJ'});
K = K(any(K==1,2),:);
Ain = [];
for i = 1:size(K,1)
  vert_i = vert(:,K(i,:));
  null_vec = null(vert_i')';
  null_vec = null_vec.*bsxfun(@times,-sign(sum(null_vec*vert,2)),ones(1,dim));
  null_vec = null_vec(max(null_vec*rays,[],2)<1e-5,:);
  Ain = [Ain;null_vec];
end
bin = zeros(size(Ain,1),1);
end