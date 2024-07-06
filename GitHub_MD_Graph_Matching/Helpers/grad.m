function val_grad = grad(X)
A = evalin('base','A');
B = evalin('base','B');
val_grad = (A^2)*X-X*(B^2);
end 

