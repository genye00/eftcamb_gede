import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root

def df(t,vec):
  omr = vec[0]
  x = vec[1]
  z = vec[2]
  lmd = vec[3]
  y = 256*xiv0*x/(z*lmd**4)
  om = 1 - omr - x**2/6. - z/3. - x**2*y
  hp = -((3*om*(1 + 6*y) + 4*omr*(1 + 6*y) + x*(x + 12*x*y + 27*x*y**2 + lmd*y*z))/(2 + 12*y + 3*x**2*y**2))
  # d = dOmega_r, dx, dy, dz, dp
  d = np.empty(len(vec))
  d[0] = -4*omr - 2*hp*omr
  d[1] = -((-((x + 9*x*y)*(3*om + 4*omr + x**2*(1 + 3*y))) + (2 - x**2*y)*(x*(3 + 9*y) + lmd*z))/(2 + 12*y + 3*x**2*y**2))
  d[2] = -2*hp*z + lmd*z*x
  d[3] = -0.25*lmd**2*x
  return d

xiv0 = -1
phi_i = 0.7

omr_i = 1-1e-9
z_i = 1e-7
l_i = 4/phi_i
x_i = (-1+np.sqrt(1-9.6*xiv0*phi_i**3))/6/xiv0/256*z_i*l_i**4

ti=0
tf=3.5
ts = np.linspace(ti,tf,num=2000)
ax = plt.figure().add_subplot(projection='3d')
num = 30
cm = plt.get_cmap('cool', num)
for i in range(num):
  # x_i = np.random.uniform(1,1)*(-0.136/xiv0*z_i)
  z_i = 10**np.random.uniform(-5,-8)
  x_c = (-1+np.sqrt(1-9.6*xiv0*phi_i**3))/6/xiv0/256*z_i*l_i**4
  x_i = x_c * np.random.normal(1, 0.3)
  # p_i = np.random.normal(0.6, 0.1) 
  sol = solve_ivp(df, t_span=[ti,tf], y0=[omr_i,x_i,z_i,l_i], t_eval=ts)
  vec = np.ones(len(sol.t))
  ax.plot(np.log10(np.abs(sol.y[1])), sol.y[2]/3, 256*xiv0*sol.y[1]/sol.y[2]/sol.y[3]**4, c=cm(i))

  ax.plot(np.log10(np.abs(sol.y[1])), 0.14*vec, 256*xiv0*sol.y[1]/sol.y[2]/sol.y[3]**4, lw=0.5, alpha=0.8, c=cm(i))
  ax.plot(-10*vec, sol.y[2]/3, 256*xiv0*sol.y[1]/sol.y[2]/sol.y[3]**4, lw=0.5, alpha=0.8, c=cm(i))
  ax.plot(np.log10(np.abs(sol.y[1])), sol.y[2]/3, 0*vec, lw=0.5, alpha=0.8, c=cm(i))
  # n0 = 500
  # n1 = 510
  # pos = lambda n: np.array([np.log10(np.abs(sol.y[1][n])),sol.y[2][n]/3, 256*xiv0*sol.y[1][n]/sol.y[2][n]/sol.y[3][n]**4])
  # ax.quiver(pos(n0)[0], pos(n0)[1], pos(n0)[2], (pos(n1)-pos(n0))[0], (pos(n1)-pos(n0))[1], (pos(n1)-pos(n0))[2], lw=5, length=1)

def f(xvec, omr):
    x = xvec[0]
    z = xvec[1]
    y = np.empty(2)
    y[0] = 18*x**14*xiv0**2 + 15*x**9*xiv0*z + 2*x**7*xiv0*(9 + 9*omr - 13*z)*z + x**4*z**2 + 16*z**3 - 2*x**2*z**2*(3 - omr + z)
    y[1] = -6*x**12*xiv0**2 - 12*x**7*xiv0*z - x**2*z**2 + 2*z**2*(1 - omr + z) + 4*x**5*xiv0*z*(3 - 3*omr + 5*z)
    return y
num = 2
xz = np.empty([num,2])
omrs = np.linspace(0,1,num=num, endpoint=True)
for i, val in enumerate(omrs):
    r = root(f, x0=[0.5, 0.2], args=(val,)) 
    xz[i,:] = r.x
# ax.plot(np.log10(xz[:,0]), xz[:,1]/3, -xiv0*xz[:,0]**5/xz[:,1],marker='.',ms=5, c='black')

# ax.plot(-10*np.ones(2), xz[:,1]/3, -xiv0*xz[:,0]**5/xz[:,1], marker='x',ms=3, c='black', alpha=0.8)
# ax.plot(np.log10(xz[:,0]), 0.14*np.ones(2), -xiv0*xz[:,0]**5/xz[:,1], marker='x',ms=3, c='black', alpha=0.8)
# ax.plot(np.log10(xz[:,0]), xz[:,1]/3, 0*np.ones(2), marker='x',ms=3, c='black', alpha=0.8)
ax.set_xlabel(r'$|\dot{\phi}/HM_p|$',fontsize=15)
ax.set_ylabel(r'$V/3M_p^2H^2$',fontsize=15)
ax.set_zlabel(r'$\xi H\dot{\phi}$',fontsize=15)
ax.set_zlim(0,0.2)
ax.set_ylim(0,0.14)
ax.set_xlim(-10,0)
ax.set_zticks([0, 0.05, 0.1, 0.15, 0.2])
ax.set_xticks([-10,-6, -4, -2, 0])
ax.set_yticks([0,0.04,0.08,0.12])
ax.set_xticklabels([r'$10^{-10}$',r'$10^{-6}$',r'$10^{-4}$',r'$10^{-2}$', '1'])

plt.show()