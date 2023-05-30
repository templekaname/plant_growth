import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
pi = np.pi

class ParamPlant:
  def __init__(self, E, rho, r_0, dr, dl, rw_stress):
    # Check E
    if (type(E) == int) or (type(E) == float):
      self.E = float(E)
    else:
      raise Exception("Datatype of E not allowed")
    # Check rho
    if (type(rho) == int) or (type(rho) == float):
      self.rho = float(rho)
    else:
      raise Exception("Datatype of rho not allowed")
    # Check r_0
    if (type(r_0) == int) or (type(r_0) == float):
      self.r_0 = float(r_0)
    else:
      raise Exception("Datatype of r_0 not allowed")
    # Check dr
    if (type(dr) == int) or (type(dr) == float):
      self.dr = float(dr)
    else:
      raise Exception("Datatype of dr not allowed")
    # Check dl
    if (type(dl) == int) or (type(dl) == float):
      self.dl = float(dl)
    else:
      raise Exception("Datatype of dl not allowed")
    # Check rw_stress
    if (type(rw_stress) == int) or (type(rw_stress) == float):
      self.rw_stress = float(rw_stress)
    else:
      raise Exception("Datatype of rw_stress not allowed")
  
  def summary(self):
    print("Summary of the plant parameters")
    print("         E = {}".format(self.E))
    print("       rho = {}".format(self.rho))
    print("       r_0 = {}".format(self.r_0))
    print("        dr = {}".format(self.dr))
    print("        dl = {}".format(self.dl))
    print(" rw_stress = {}".format(self.rw_stress))

class ParamFoliage:
  def __init__(self, dh, h_end):
    # Check dh
    if (type(dh) == int) or (type(dh) == float):
      self.dh = float(dh)
    else:
      raise Exception("Datatype of dh not allowed")
    # Check h_end
    if (type(h_end) == int) or (type(h_end) == float):
      self.h_end = float(h_end)
    else:
      raise Exception("Datatype of h_end not allowed")

  def summary(self):
    print("Summary of foliage load")
    print("        dh = {}".format(self.dh))
    print("     h_end = {}".format(self.h_end))

class ParamAngle:
  def __init__(self, theta_0, theta_C, theta_P):
    # Check theta_0
    if (type(theta_0) == int) or (type(theta_0) == float):
      self.theta_0 = float(theta_0)
    else:
      raise Exception("Datatype of theta_0 not allowed")
    # Check theta_C
    if (type(theta_C) == int) or (type(theta_C) == float):
      self.theta_C = float(theta_C)
    elif (theta_C == "straight"):
      self.theta_C = "straight"
    else:
      raise Exception("Datatype of theta_C not allowed")
    # Check theta_P
    if (type(theta_P) == int) or (type(theta_P) == float):
      self.theta_P = float(theta_P)
    else:
      raise Exception("Datatype of theta_P not allowed")
  
  def summary(self):
    print("Summary of preferred angles")
    print("   theta_0 = {}".format(self.theta_0))
    print("   theta_C = {}".format(self.theta_C))
    print("   theta_P = {}".format(self.theta_P))



class GrowthModel:
  def __init__(self, N:int, pm:ParamPlant, pf:ParamFoliage, pa:ParamAngle):
    self.N = N
    self.pm = pm
    self.pf = pf
    self.pa = pa
    self.l = np.zeros((N+1,N+1))
    self.theta = np.zeros((N+1,N+1))
    self.x = np.zeros((N+1,N+2))
    self.y = np.zeros((N+1,N+2))
    # Set the first coordinates
    l = self.l
    x = self.x
    y = self.y
    theta = self.theta
    dl = self.pm.dl 
    t0 = self.pa.theta_0

    x[0,0] = 0
    x[0,1] = dl * np.cos(t0)
    y[0,0] = 0
    y[0,1] = -dl * np.sin(t0)
    l[0,0] = ((x[0,1]-x[0,0])**2 + (y[0,1]-y[0,0])**2)**0.5
    theta[0,0] = t0
    

  def summary(self):
    self.pm.summary()
    self.pf.summary()
    self.pa.summary()

  def grow(self,k):
    print("Starting growth step {}".format(k))
    E = self.pm.E
    rho = self.pm.rho
    r_0 = self.pm.r_0
    dr = self.pm.dr
    dl = self.pm.dl
    dh = self.pf.dh
    h_end = self.pf.h_end
    theta_C = self.pa.theta_C
    theta_P = self.pa.theta_P
    rw_stress = self.pm.rw_stress

    self.l[k,:k] = deepcopy(self.l[k-1,:k])
    self.l[k,k] = dl
    self.theta[k,:k] = deepcopy(self.theta[k-1,:k])
    if theta_C == "straight":
      self.theta[k,k] = self.theta[k,k-1]
    else:
      self.theta[k,k] = theta_C
    self.x[k,:k+1] = deepcopy(self.x[k-1,:k+1])
    self.x[k,k+1] = self.x[k,k] + self.l[k,k]*np.cos(self.theta[k,k])
    self.y[k,:k+1] = deepcopy(self.y[k-1,:k+1])
    self.y[k,k+1] = self.y[k,k] + self.l[k,k]*np.sin(self.theta[k,k])

    l = deepcopy(self.l[k,:k+1])
    theta = deepcopy(self.theta[k,:k+1])
    x = deepcopy(self.x[k,:k+2])
    y = deepcopy(self.y[k,:k+2])

    dp = np.zeros(k+2)
    domega = np.zeros(k+1)
    dMo = np.zeros(k+2)
    dMs = np.zeros(k+1)
    I = np.zeros(k+1)
    alpha = np.zeros(k+1)
    C = np.zeros(k+1)

    newx = np.zeros(k+2)
    newy = np.zeros(k+2)
    newl = np.zeros(k+1)
    newtheta = np.zeros(k+1)

    for i in range(k):
      S1 = np.sum(l[i:k]*(2*r_0 + (2*k - 2*i - 1)*dr))
      S3 = dh*(k-i+1)
      dp[i] = rho*pi*dr*S1 + rho*pi*r_0**2*dl + S3 + h_end
    dp[k] = rho*pi*r_0**2*dl + dh + h_end
    dp[k+1] = h_end

    for i in range(k):
      domega[i] = rho*pi*dr*(2*r_0 + (2*k - 2*i - 1)*dr)*np.cos(theta[i])
    if theta_C == "straight":
      domega[k] = rho*pi*r_0**2*np.cos(theta[k-1])
    else:
      domega[k] = rho*pi*r_0**2*np.cos(theta_C)

    for i in range(k+1):
      r_bar = r_0 + (k - i - 1)*dr
      dMs[i] = -0.5*r_bar**2*pi*dr*rw_stress*np.sin(theta[i] - theta_P)
    dMs[k] = 0

    for i in range(k+1):
      dMo[i] = -np.sum(0.5*domega[i:]*l[i:]**2 + dp[i+1:]*l[i:]*np.cos(theta[i:]) - dMs[i:])
      dMo[k+1] = 0

    for i in range(k+1):
      I[i] = pi/4*(r_0 + (k-i)*dr)**4

    for i in range(k+1):
      trm1 = domega[:i+1]*l[:i+1]**3/(6*E*I[:i+1])
      trm2 = dp[1:i+2]*l[:i+1]**2/(2*E*I[:i+1])
      trm3 = -dMo[1:i+2]*l[:i+1]/(E*I[:i+1])
      trm4 = -dMs[:i+1]*l[:i+1]/(E*I[:i+1])
      alpha[i] = np.sum(trm1+trm2+trm3+trm4)

    C = -domega*l**4/(24*E*I) - dp[1:]*l**3*np.cos(theta)/(6*E*I) + dMo[1:]*l**2/(2*E*I) + dMs*l**2/(2*E*I) + alpha*l

    for i in range(0,k+1):
      rot = np.array([[np.cos(theta[i]),np.sin(theta[i])],[-np.sin(theta[i]),np.cos(theta[i])]])
      v = np.array([l[i],C[i]])
      prod = np.dot(rot,v) + np.array([newx[i],newy[i]])
      newx[i+1] = prod[0]
      newy[i+1] = prod[1]

    self.x[k,:k+2] = newx
    self.y[k,:k+2] = newy
    self.l[k,:k+1] = ((newx[1:]-newx[:-1])**2 + (newy[1:]-newy[:-1])**2)**0.5
    for i in range(k+1):
      if newx[i+1] - newx[i] < 0 and newy[i+1] - newy[i] < 0 :
        self.theta[k,i] = np.arctan(-(newy[i+1]-newy[i])/(newx[i+1]-newx[i])) - pi
      elif newx[i+1] - newx[i] < 0 and newy[i+1] - newy[i] > 0:
        self.theta[k,i] = np.arctan(-(newy[i+1]-newy[i])/(newx[i+1]-newx[i])) + pi
      else:
        self.theta[k,i] = np.arctan(-(newy[i+1]-newy[i])/(newx[i+1]-newx[i]))

    print("Finished growth step {}".format(k))

  def simulate(self):
    for k in range(1,self.N+1):
      self.grow(k)
    print("Simulation finished")

  def reset(self):
    N = self.N
    self.l = np.zeros((N+1,N+1))
    self.theta = np.zeros((N+1,N+1))
    self.x = np.zeros((N+1,N+2))
    self.y = np.zeros((N+1,N+2))
    # Set the first coordinates
    l = self.l
    x = self.x
    y = self.y
    theta = self.theta
    dl = self.pm.dl 
    t0 = self.pa.theta_0

    x[0,0] = 0
    x[0,1] = dl * np.cos(t0)
    y[0,0] = 0
    y[0,1] = dl * np.sin(t0)
    l[0,0] = ((x[0,1]-x[0,0])**2 + (y[0,1]-y[0,0])**2)**0.5
    theta[0,0] = t0
    print("Reset data to default")

  def track_growth(self,frames):
    plt.gca().invert_yaxis()
    for k in frames:
      plt.plot(self.x[k,:k+2],self.y[k,:k+2],alpha=0.5,c='k')
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.draw()
    plt.grid()
