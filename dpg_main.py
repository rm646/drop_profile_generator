#---modules
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import sys
import argparse
import math
from scipy.integrate import odeint
import scipy.integrate as integrate
import scipy.special as special
import datetime
import imageio

#---argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('-nw','--needlewidth',help="Diameter that drop will attach to, in mm.", type=float)
parser.add_argument('-ift','--interfacialtension', help="Interfacial tension of drop, in mN/m.", type=float)
parser.add_argument('-dd','--densitydiff',help="Density_drop - Density_surr, in g/mL. Must be positive.", type=float)
parser.add_argument('-R0','--dropdim',help="Drop dimension", type=float)
parser.add_argument('-N','--gridwidth', help="Number of pixels in X direction of output image.", type=int)
parser.add_argument('-M','--gridheight', help="Number of pixels in Y direction of output image.", type=int)
parser.add_argument('-s','--scale', help="Image scale in px/mm.", type=float)
args=parser.parse_args()

#---set defaults
defaultN=780
defaultM=580
defaultScale=85 #px/mm
defaultnw=1.83 #mm
defaultdd=0.2509 #g/mL
defaultR0=2.4 #mm
defaultift=53.5 #mN/m
tolerance_nw = 0.001 #in mm #values of r within this tolerance of the needle width will be considered equal to the needle width

if args.needlewidth ==None: args.needlewidth = defaultnw
elif args.needlewidth<=0 :
    print "Needle width must be >0. Instead using default value of ",defaultnw
    args.needlewidth = defaultnw

if args.densitydiff ==None: args.densitydiff = defaultdd
elif args.densitydiff<=0 :
    print "Density difference must be >0. Instead using default value of ",defaultdd
    args.densitydiff = defaultdd

if args.dropdim ==None: args.dropdim = defaultR0
elif args.dropdim<=0 :
    print "Drop dimension must be >0. Instead using default value of ",defaultR0
    args.dropdim = defaultR0

if args.interfacialtension ==None: args.interfacialtension = defaultift
elif args.interfacialtension<=0 :
    print "Interfacial tension must be >0. Instead using default value of ",defaultift
    args.dropdim = defaultift

if args.gridwidth ==None : args.gridwidth = defaultN
elif args.gridwidth<=0 :
    print "N must be >0. Instead using default value of ",defaultN
    args.gridwidth = defaultN

if args.gridheight ==None : args.gridheight = defaultM
elif args.gridheight<=0 :
    print "M must be >0. Instead using default value of ",defaultM
    args.gridheight = defaultM

if args.scale ==None : args.scale = defaultScale
elif args.scale<=0 :
    print "Scale must be >0. Instead using default value of ",defaultScale
    args.scale = defaultScale

#---Calculate Bo
def Bo(dd, ift, R0):
    #convert dd from g/mL to kg/m^3
    #convert ift from mN/m to N/m
    #convert R0 from mm to m
    return ( (dd*1000.0*9.81*(R0/1000.0)*(R0/1000.0))/(ift/1000.0) )

#---Define function to take
#a) list of current y_i values
#b) current t
#c) any parameters needed to evaluate function
def func(y_i, t, params):
    #unpack parameters
    #y0 = phi, y1 = r_bar, y2 = z_bar, where bar denotes scaling by R0 to make it dimensionless
    Bo = params
    if y_i[1]==0.0: dy0_dt = 0
    else: dy0_dt = 2 - Bo*y_i[2] - math.sin(y_i[0])/y_i[1]
    dy1_dt = math.cos(y_i[0])
    dy2_dt = math.sin(y_i[0])
    return [dy0_dt, dy1_dt, dy2_dt]

#---Define function to calculate dimensionless drop volume (scaled by R0^3) from (r,z) profile, assuming symmetry about the z axis
def V_bar(r_bar,z_bar):
    f = math.pi * r_bar * r_bar
    result = np.trapz(f,z_bar)
    return result

#---Define function to calculate drop surface area (scaled by R0^2) from (r,z) profile
def A_bar(r_bar,z_bar):
    dr_dz = np.diff(r_bar)/np.diff(z_bar) #but this is one element shorter than r or z and centred at z[i+1]+z[i]/2
    r2 = (r_bar[:-1] + r_bar[1:]) / 2.0
    z2 = (z_bar[:-1] + z_bar[1:]) / 2.0
    f = r2*np.sqrt(1+dr_dz*dr_dz)
    result = 2.0*math.pi * np.trapz(f,z2)
    return result

#---Define initial conditions (phi,r,z)
y0 = [0.0, 0.0, 0.0]

#---Define timepoints to solve for
t = np.linspace(0, 10, 10000)

#---Integrate and solve
#Bondnum = Bo(args.dd, args.ift, args.R0)
Bondnum_calc = Bo(args.densitydiff, args.interfacialtension, args.dropdim)
sol = odeint(func, y0, t, args=(Bondnum_calc,))


#---Find z index where needle width equals drop profile position
r = args.dropdim * sol[:,1] #recall sol[:,1] is r_bar
i_rns = np.where(abs(r-(args.needlewidth/2.0)) < tolerance_nw)[0] #find r values within tolerance of half needle width
#difficult to narrow down which one to choose: it depends on the drop volume.
#Let's restrict ourselves to drops where more than a hemisphere of fluid is pushed out
#This means it's in the phi>pi/2 region
#Take the first point within tolerance in this region
i_phi_g_90 = np.where(sol[:,0]> (math.pi/2.0))[0] #indices in that phi region
i_rns = i_rns[np.where(i_rns > i_phi_g_90[0])[0]]
i_rn = None #initial value for the index of the point the needle will touch. If the needle is outside tolerance, this value will remain unchanged
if i_rns.size!=0: i_rn = i_rns[np.argmin(i_rns)]
else:
    print "Needle width is too small to have a drop which is more than a hemisphere. Such cases are currently not allowed."
    exit()

#---Form a pixel image of the drop and needle
#do one half, then mirror
image = np.zeros((args.gridheight,args.gridwidth))
image += 255 #make white
#for pixels within the profile, make black.
scale = args.scale
nheight = args.gridheight/10
#for each point in the solution, ceil to the nearest pixel
for point in sol[:i_rn,:]:
    r = point[1]*args.dropdim
    z = point[2]*args.dropdim
    #nearest pixel values
    r_px = int(r*scale)
    z_px = int(z*scale)
    #shift so that drop begins at gridheight/10
    z_top = int(sol[i_rn,2]*args.dropdim*scale)
    shift = (args.gridheight-z_top-nheight)
    #handle profile overrun
    if r_px > args.gridwidth/2 or shift<0:
        print "Profile exceeded frame. Adjust image scale or make a smaller drop."
        exit()
    #set image
    image[(z_px+shift)][args.gridwidth/2:((args.gridwidth/2)+r_px)] = 0

#add needle
for j in range(0,args.gridheight):
    if (j > args.gridheight-nheight):#sol[i_rn,2]*args.dropdim*scale): #if above needle contact, fill black to needle width
        needle_edge_px = int(scale*(args.needlewidth/2.0))
        image[j][args.gridwidth/2:(args.gridwidth/2 + needle_edge_px)] = 0

#mirror around centre line
for j in range(0,args.gridheight):
    for i in range(0, args.gridwidth/2):
        image[j][args.gridwidth/2 - i] = image[j][args.gridwidth/2 + i]

#---Add 100px buffer to drop bottom
#buffer = np.zeros((100,args.gridwidth)) + 255
#image = np.vstack([buffer,image])

#---Flip vertically
image = np.flip(image, 0)

#---Quick plot

V = V_bar(sol[:i_rn,1],sol[:i_rn,2])*(args.dropdim**3)
A = A_bar(sol[:i_rn,1],sol[:i_rn,2])*(args.dropdim**2)
today = str(datetime.date.today())
filename = "%s_gw%d_gh%d_s%.2f_nw%.2f_ift%.2f_dd%.4f_R0%.3f_V%.2f_A%.2f.png" % (\
today,\
args.gridwidth,\
args.gridheight,\
args.scale,\
args.needlewidth,\
args.interfacialtension,\
args.densitydiff,\
args.dropdim,\
V,\
A,)
int_image = image.astype(np.uint8)
imageio.imwrite(filename,int_image)
