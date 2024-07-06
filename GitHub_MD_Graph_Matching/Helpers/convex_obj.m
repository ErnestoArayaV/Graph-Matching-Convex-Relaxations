function Val=convex_obj(A,B,X)
Val=norm(A*X-X*B,'fro')^2;
