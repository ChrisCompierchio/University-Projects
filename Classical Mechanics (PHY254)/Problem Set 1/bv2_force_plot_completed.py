from pylab import *
import numpy as np
# Solution to 2020 PHY254 Homework 1 question 4 part (b).  
# Plotting analytic solutions of
# a(t), v(t) and x(t) for F = - b v^2 = ma 

#constant for force, just to pick a number
b = 0.25

# mass of the particle.  Things only depend on b/m.
m = 1.0

#initial conditions
x0= 0.0
v0= 3.0  # just to choose a number

#time step
dt =  0.01 # you can vary this to see the effects.  0.5 is too long.  0.01 is good.
t = arange(0.0,10.0,dt)      # an array of times to plot up to 10 seconds

npts = len(t)  # number of time steps

#  analytic expressions: you found these in part (a) of the question
#  full arrays will be generated automatically, with entries for each t.

#velocity
v = v0/(((b*t)/m)+1) # insert your expression here

#position
x = (m/b)*log(((v0*b*t)/m)+1) # insert your expression here

#acceleration
a = -(v0**2)*(b/m)*(1/((((v0*b*t)/m)+1)**2)) # insert your expression here

# numerical solution via time stepping

an = zeros(npts)  # makes arrays full of zeros
vn = zeros(npts)
xn = zeros(npts)

vn[0] = v0  # set initial condition as first member of array

for i in range(1,npts):    # range function omits last index at npts, starts at 1 not zero
    vn[i] = vn[i-1]+a[i-1]*dt # insert your expression here

xn[0] = x0  # set initial condition as first member of array

for i in range(1,npts): 
    xn[i] = xn[i-1] + v[i]*dt # insert your expression here

for i in range(0,npts):     # go over the whole range 0 .. npts-1
    an[i] = -(b/m)*(v[i]**2)# insert your expression here


#set up first plotting window
subplot(3,1,1)
#plot a(time), label it, turn on the grid.

plot(t,a, t,an) # comparison plot for analytic and numerical results
ylabel('a(t)')
grid('on')

#set up next plotting window, plot v(time)
subplot(3,1,2)
#plot(t, v)
plot(t,v, t,vn) # comparison plot for analytic and numerical results
ylabel('v(t)')
grid('on')

#and x(time); label y and x axis for bottom plot.
subplot(3,1,3)
#plot(t, x)
plot(t,x, t, xn)  # comparison plot for analytic and numerical results
legend((r'exact',r'numerical approx'),loc="lower right")
ylabel('x(t)')
xlabel('t')
grid('on')

# super title (goes at top of window)
# formatting trick puts dt in the title
suptitle('Chris Compierchios bv squared force solutions for time step ' + '%.3f' % dt )
# save a pdf of the figure
# formatting trick puts dt in the title
savefig('Chris_Compierchio_bv2_force_dt_'+'%.3f' % dt + '.pdf')
# show the plot on the screen
show()
