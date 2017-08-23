function f = PoiSolv1D(src,LBC,RBC,h)
	N = size(src);
	S = N(1)*N(2);
	src(1) = src(1) + (2*LBC)/(h^2);             % corr. in source at LB
	src(end) = src(end) + (2*RBC)/(h^2);         % corr. in source at RB
	A = spdiags(ones(S,1)*[1 -2 1],-1:1,S,S);
	A(1,1) = -3;    A(end,end) = -3;             % Dirichlet b.c.
	A1D = A/(h^2);
	f = -A1D\(src');
	f = f';
end