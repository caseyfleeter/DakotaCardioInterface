import numpy as np

# Python Runge-Kutta 4 method
#    dydt is a (system of) differential eqn(s) of the form dydt = f(t,y,params)
#    t0 is the starting time, h is the time step, n is the number of time steps
#    y0 is the initial condition (vector)
#    params are other parameters needed for the dydt function
def rk4(dydt,t0,h,n,y0,params=['None']):
	hh = h/2.
	# solution (vector) at inital time
	w = [y0]
	# time steps
	t = [t0]
	for i in xrange(0,n):
		# if np.mod(i,100) == 0:
		# 	print 'CurrentTime %d/%d ' % (i,n)

		# keep track of time steps
		t.append(t0 + (i+1)*h)
		#print('rk4 '+ str(t0 + (i+1)*h))

		# Implement RK4 algorithm
		# if dydt function does not require additional inputs
		if params == ['None']:
			k1 = [h*dy for dy in dydt(t[i],w[i])]
			
			wtemp = [ww + kk1/2. for ww,kk1 in zip(w[i],k1)]
			k2 = [h*dy for dy in dydt(t[i]+hh,wtemp)]

			wtemp = [ww + kk2/2. for ww,kk2 in zip(w[i],k2)]
			k3 = [h*dy for dy in dydt(t[i]+hh,wtemp)]

			wtemp = [ww + kk3 for ww,kk3 in zip(w[i],k3)]
			k4 = [h*dy for dy in dydt(t[i]+h,wtemp)]
		# if dydt does requires additional inputs
		else:
			k1 = [h*dy for dy in dydt(t[i+1],w[i],*params)]
			
			wtemp = [ww + kk1/2. for ww,kk1 in zip(w[i],k1)]
			
			k2 = [h*dy for dy in dydt(t[i+1]+hh,wtemp,*params)]

			wtemp = [ww + kk2/2. for ww,kk2 in zip(w[i],k2)]
			
			k3 = [h*dy for dy in dydt(t[i+1]+hh,wtemp,*params)]

			wtemp = [ww + kk3 for ww,kk3 in zip(w[i],k3)]
			
			k4 = [h*dy for dy in dydt(t[i+1]+h,wtemp,*params)]
		# print([ww + (1/6.)*(kk1+2.*kk2+2.*kk3+kk4) for ww,kk1,kk2,kk3,kk4 in zip(w[i],k1,k2,k3,k4)])
		w.append([ww + (1/6.)*(kk1+2.*kk2+2.*kk3+kk4) for ww,kk1,kk2,kk3,kk4 in zip(w[i],k1,k2,k3,k4)])
		
	return w,t
