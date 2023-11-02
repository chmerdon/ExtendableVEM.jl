### find the boundary edge number
function findm(a,b)
	an = size(a,2)
	bn = size(b,2)
	c = zeros(Int32, an)
	for i in 1:an
		aa = a[:,i]
		for j in 1:bn
			if aa[1] == b[1,j]
				if aa[2] == b[2,j]
					c[i] = j
				end
			end
		end
	end
	return c
end
