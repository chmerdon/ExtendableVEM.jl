### error function for scalar  
### eu: exact solution function, a: numerical solution coefficient, p:0 L2, p:1 H1,
### center: element center coordinate, diameter: element diameter, k: error order,
function er_f(x, eu, a, p, center, diamter, k)  ##  k max 2
	# k == 2 npolys = (k+2)*(k+1) ÷ 2  
	if k==1
		e = a[1] .* poly_base(x[1],x[2],1,p,center,diamter) + a[2] .* poly_base(x[1],x[2],2,p,center,diamter) + a[3] .* poly_base(x[1],x[2],3,p,center,diamter)
	elseif k==2
		e = a[1] .* poly_base(x[1],x[2],1,p,center,diamter) + a[2] .* poly_base(x[1],x[2],2,p,center,diamter) + a[3] .* poly_base(x[1],x[2],3,p,center,diamter) + a[4] .* poly_base(x[1],x[2],4,p,center,diamter) + a[5] .* poly_base(x[1],x[2],5,p,center,diamter) + a[6] .* poly_base(x[1],x[2],6,p,center,diamter)
	end
		# for i in 1:npolys
	# 	e .-= a[i] .* poly_base(x[1],x[2],i,p,center,diamter)
	# end
	er = eu(x) - e
	er = er .^ 2
	return er
end

### error function for vector  
### eu: exact solution function, a: numerical solution coefficient, p:0 L2, p:1 H1,
### center: element center coordinate, diameter: element diameter, k: error order,
function er_f_2(x, eu, a, p, center, diamter, k)  ## k=2
	# k==2 npolys = (k+2)*(K+1) ÷ 2
	e = a[1] .* poly_base(x[1],x[2],1,p,center,diamter) + a[2] .* poly_base(x[1],x[2],2,p,center,diamter) + a[3] .* poly_base(x[1],x[2],3,p,center,diamter) + a[4] .* poly_base(x[1],x[2],4,p,center,diamter) + a[5] .* poly_base(x[1],x[2],5,p,center,diamter) + a[6] .* poly_base(x[1],x[2],6,p,center,diamter)

	er1 = eu(x,1) - e

	# for i in 1:npolys
	# 	er1 .-= a[i] .* poly_base(x[1],x[2],i,p,center,diamter)
	# end
	er1 = er1 .^ 2

	er2 = eu(x,2) - e
	# for i in 1:npolys
	# 	er2 .-= a[i] .* poly_base(x[1],x[2],i,p,center,diamter)
	# end
	er2 = er2 .^ 2

	er = sqrt.(er1.^2 + er2.^2)
	return er
end

### error function for scalar  
### eu: exact solution function, a: numerical solution coefficient, p:0 L2, p:1 H1,
### center: element center coordinate, diameter: element diameter, k: error order,
function er_fp(x, eu, a, p, center, diamter, k)  ## k=1
	# npolys = (k+2)*(K+1)/2
	e = a[1] .* poly_base(x[1],x[2],1,p,center,diamter) + a[2] .* poly_base(x[1],x[2],2,p,center,diamter) + a[3] .* poly_base(x[1],x[2],3,p,center,diamter)
	er = eu(x,3) - e
	# for i in 1:npolys
	# 	er .-= a[i] .* poly_base(x[1],x[2],i,p,center,diamter)
	# end
	er = er .^ 2
	return er
end