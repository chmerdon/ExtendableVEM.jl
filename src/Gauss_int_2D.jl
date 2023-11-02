### Gauss integration for 2D element  
### f: function, k: accuracy, T: element nodes coordinate, center: element center coordinate
function Gauss_int_2D_tri(f::Function,k,T,center)
	el_int = 0
	el_n = size(T,2)
	Tb = [T;T[:,2:end] T[:,1]]
	for n in 1:el_n
		int = 0
		x1 = Tb[1,n]
		x2= Tb[3,n]
		x3 = center[1]
		y1 = Tb[2,n]
		y2 = Tb[4,n]
		y3 = center[2]
		J=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
	
		if k==3
			x=[1/2 0;1/2 1/2;0 1/2]
			A=[1/6 1/6 1/6]
		elseif k==4
			x=[(1/sqrt(3)+1)/2 (1-1/sqrt(3))*(1+1/sqrt(3))/4;(1/sqrt(3)+1)/2 (1-1/sqrt(3))*(1-1/sqrt(3))/4;(-1/sqrt(3)+1)/2 (1+1/sqrt(3))*(1+1/sqrt(3))/4;(-1/sqrt(3)+1)/2 (1+1/sqrt(3))*(1-1/sqrt(3))/4]
			A=[(1-1/sqrt(3))/8 (1-1/sqrt(3))/8 (1+1/sqrt(3))/8 (1+1/sqrt(3))/8]
		elseif k==9
			x=[(1+0)/2 (1-0)*(1+0)/4;(1+sqrt(3/5))/2 (1-sqrt(3/5))*(1+sqrt(3/5))/4;(1+sqrt(3/5))/2 (1-sqrt(3/5))*(1-sqrt(3/5))/4;(1-sqrt(3/5))/2 (1+sqrt(3/5))*(1+sqrt(3/5))/4;(1-sqrt(3/5))/2 (1+sqrt(3/5))*(1-sqrt(3/5))/4;(1+0)/2 (1-0)*(1+sqrt(3/5))/4;(1+0)/2 (1-0)*(1-sqrt(3/5))/4;(1+sqrt(3/5))/2 (1-sqrt(3/5))*(1+0)/4;(1-sqrt(3/5))/2 (1+sqrt(3/5))*(1+0)/4]
			A=[64/81*(1-0)/8 100/324*(1-sqrt(3/5))/8 100/324*(1-sqrt(3/5))/8 100/324*(1+sqrt(3/5))/8 100/324*(1+sqrt(3/5))/8 40/81*(1-0)/8 40/81*(1-0)/8 40/81*(1-sqrt(3/5))/8 40/81*(1+sqrt(3/5))/8]
		end
		xx = x1.+(x2-x1).*x[:,1].+(x3-x1).*x[:,2]
		yy = y1.+(y2-y1).*x[:,1].+(y3-y1).*x[:,2]

		for j in 1:k
			int += f(xx[j],yy[j]).*A[j]
		end
		int = J*int
		el_int += int
	end
	return el_int
end