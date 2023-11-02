
### basic 2D function for P_k in each element  
### 2D:[x,y], poly_u: α-p_α, poly_k:0 p, poly_k:[1 2] ∇p, 
### center: element center coordinate, diameter: element diameter
function poly_base(x,y,poly_u,poly_k,center,diameter)
 
	m1 = 1
	m2 = (x-center[1]) ./ diameter
	m3 = (y-center[2]) ./ diameter
	dm = 1 / diameter

	if poly_k == 0 
		if poly_u == 1
			m = m1
		elseif poly_u == 2
			m = m2
		elseif poly_u == 3
			m = m3
		elseif poly_u == 4
			m = m2.^2
		elseif poly_u == 5
			m = m2.*m3
	    elseif poly_u == 6
			m = m3.^2
		end
	elseif poly_k == 1
		if poly_u == 1
			m = 0;
		elseif poly_u == 2
			m = dm;
		elseif poly_u ==3
			m = 0;
		elseif poly_u ==4
			m = 2*m2.*dm;
		elseif poly_u ==5
			m = dm.*m3;
		elseif poly_u ==6
			m = 0;	
		end
	elseif poly_k == 2
		if poly_u == 1
			m = 0;
		elseif poly_u == 2
			m = 0;
		elseif poly_u ==3
			m = dm;
		elseif poly_u ==4
			m = 0;
		elseif poly_u ==5
			m = m2.*dm;
		elseif poly_u ==6
			m = 2*m3.*dm;	
		end
	end

	return m
end