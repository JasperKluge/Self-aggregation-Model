'''imports/modules'''

import numpy as np
import random as rn
import matplotlib.pylab as plt
import copy
import math 


'''basic variables'''

clch = False  #switch for the climate change model

R_0 = 5       #[km] Radius of the cells  

dl = 5        #[km] maximum displacement

N = 5         #number of starting cells

upbound = 100 #uppper boundry for the precipitation

colors = 'royalblue' # color 

if clch: 
    N = 39
    upbound= 186
    colors = 'navy'
    




r0 = np.sqrt(300)    #[km] reference for the wind

dr = 1 #[km] step length for the rad. prec.

L = 100 #[km] Domain length


'''less important starting variables'''

alph = 0.36    #direct wind constant

epsilon = 1e-10

global t_tot 
t_tot = 0   #time for calculation

dist_tot = 0

cs = []      #List of the cells

t_save = [] #cummulated time steps

t_saveinc =[] #Time steps

dist_save =[] #distanz per step

rad_prec = [] # save of the radial prec

aggrate= np.zeros(10**5) 

span = np.linspace(0, 10**5, 10**5)

end_inter = 10000


'''functions and classes'''

class ConvCells : # erstellt die Klasse Zelle die charakterisierende Variablen besitzt
    def __init__(self, P, rad, y, x, dy ,dx, mer,mer_w):
        self.P = P # Niederschlag
        self.rad = rad #Radius der Zelle
        self.y = y_zelle #Y-Position der Zelle 
        self.x = x_zelle #X-Position der Zelle
        self.dx = dx
        self.dy = dy
        self.mer = mer #gibt an ob die Zelle schon mal rezessiv verschmolzen ist 
        self.mer_w = mer_w #gibt an mit welchen Zellen die Zelle schon mal dominant verschmolzen ist 
        self.listx = [] #Liste der X-Positionen der Zellen
        self.listy = [] #Liste der Y-Positionen der Zellen
        self.listwind = [] #Liste der Windstärken jeder Zelle
        self.radprec_time = [] #radialer Niederschlag für L/2 mit abstand dr ueber Zeitschritte
        self.mergepoint = [] #punkte wo Zellen verschmelzen 
        self.listrad = []



def wind(t1, t2): 
   
    W_i = alph*cs[t2].P*cs[t2].rad**2/R_0**2*np.exp(-((dist(t1,t2)**2)/(2*r0**2)))
    return (cs[t2].x-cs[t1].x)*W_i/dist(t1,t2), (cs[t2].y-cs[t1].y)*W_i/dist(t1,t2)


def dist(t1,t2):
    return np.sqrt((cs[t1].x-cs[t2].x)**2+(cs[t1].y-cs[t2].y)**2)



def checkmerge(t1,t2):
    if t1 != t2:
        if cs[t1].mer and cs[t2].mer:
            if dist(t1,t2) <= (cs[t1].rad + cs[t2].rad):
                return True
    else:
        return False 

    
def number_active_cells():
    ac = 0
    for i in range(0,N):
        if cs[i].mer == True:
            ac += 1
    return ac 
            
def mergers():
    for o in range(0,N):
        for j in range(o+1,N):
            if checkmerge(o,j) is True:
                merging_cells.append([o,j])
    
    for i in range(0,len(merging_cells)):
        merge(merging_cells[i][0],merging_cells[i][1],n)
    
    
def merge(t1,t2,n):
    
    if t1 != t2:
        if t1>t2:
            t1,t2 = t2, t1
        cs[t1].P = (cs[t1].P*cs[t1].rad**2+cs[t2].P*cs[t2].rad**2)/(cs[t1].rad**2+cs[t2].rad**2)
        cs[t2].P = 0
        cs[t2].mer = False 
        cs[t1].mer_w.append(t2) 
        for s in cs[t2].mer_w:   
            cs[t1].mer_w.append(s)   
        cs[t1].x = (cs[t1].rad**2*cs[t1].x + cs[t2].rad**2*cs[t2].x)/(cs[t1].rad**2+cs[t2].rad**2)
        cs[t2].x = cs[t1].x
        cs[t1].y = (cs[t1].rad**2*cs[t1].y + cs[t2].rad**2*cs[t2].y)/(cs[t1].rad**2+cs[t2].rad**2)
        cs[t2].y = cs[t1].y
        cs[t1].rad = np.sqrt(cs[t1].rad**2+cs[t2].rad**2)
        cs[t2].rad = cs[t1].rad
        cs[t1].mergepoint.append([cs[t1].x, cs[t1].y])
        for inner in merging_cells:
            for index in range(len(inner)):
                if inner[index] == t2:
                    inner[index] = t1
    return

def calculating_displacement():
    for o in range(0,N):
        wind_dx =0.
        wind_dy =0.
        
        if cs[o].mer is True:
            for j in range(0,N):
                if cs[j].mer and o != j:
                    windx , windy  = wind(o,j)
                    wind_dx += windx
                    wind_dy += windy
       
        
        cs[o].dx = copy.deepcopy(wind_dx)
        cs[o].dy =  copy.deepcopy(wind_dy)
        wind_v.append(np.sqrt(wind_dy**2+wind_dx**2))
    

def find_dt(v,laenge):
    v_max = v[0]
    for i in range(len(v)):
        if v_max < v[i]:
            v_max =v[i]
        
    if abs(v_max) < epsilon: 
       
        return 0
    else: 
        return laenge/v_max

def displacement(dt, t_tot):
    for f in range(0,N):
        
        if cs[f].mer:
            
            cs[f].x += cs[f].dx*dt
            cs[f].y += cs[f].dy*dt
            
            for z in range(len(cs[f].mer_w)):
            
                cs[cs[f].mer_w[z]].x = cs[f].x
                cs[cs[f].mer_w[z]].y = cs[f].y
            
        cs[f].listx.append(cs[f].x)
        cs[f].listy.append(cs[f].y)
        cs[f].listwind.append(np.sqrt((cs[f].dx)**2+cs[f].dy**2))
   
    t_tot =dt + t_save[-1]
    
    t_save.append(t_tot)
    t_saveinc.append(dt)

def radial_prec(n):
    prec = []
    for q in range(1,int((L/2)+1)):
        r_p = 0
        for i in range(0,N):
            if cs[i].mer:
                r_p +=(check_overlap(i, q))
        prec.append(r_p/n)
    rad_prec.append(prec)
    

def check_overlap(t1,radi):
    dphi_P = 0
    for t2 in range(0,N):
        d = dist(t1,t2)
        if cs[t2].mer == False:
            continue
        else:
            if (d >= (radi + cs[t2].rad)) or d<= abs(radi-cs[t2].rad):
                continue


    
            else:
                a = (radi**2 - cs[t2].rad**2 + d**2) / (2 * d)
                h = np.sqrt(radi**2 - a**2)

                px = cs[t1].x + a * (cs[t2].x - cs[t1].x) / d
                py = cs[t1].y + a * (cs[t2].y - cs[t1].y) / d
    
                sx1 = px + h * (cs[t2].y - cs[t1].y) / d
                sy1 = py - h * (cs[t2].x - cs[t1].x) / d
                sx2 = px - h * (cs[t2].y - cs[t1].y) / d
                sy2 = py + h * (cs[t2].x - cs[t1].x) / d
                phi1 = math.atan2(sy1-cs[t1].y ,sx1-cs[t1].x)
                phi2 = math.atan2(sy2-cs[t1].y ,sx2-cs[t1].x)
                dphi = abs(phi1-phi2)
                
                if dphi > np.pi:
                    dphi = 2*np.pi - dphi
                
                dphi_P += dphi*cs[t2].P
    return dphi_P/(2*np.pi-schnittpunkte_domain(t1, radi))
    
    
            
def schnittpunkte_domain(t1,radi, domain_size=L):
    schnittpunkte = []
    
    for edge_x in [0, domain_size]:
        dy = np.sqrt(radi**2 - (edge_x - cs[t1].x)**2) if abs(edge_x - cs[t1].x) <= radi else None
        if dy is not None:
            schnittpunkte.append((edge_x, cs[t1].y + dy))
            schnittpunkte.append((edge_x, cs[t1].y - dy))
    
    for edge_y in [0, domain_size]:
        dx = np.sqrt(radi**2 - (edge_y - cs[t1].y)**2) if abs(edge_y - cs[t1].y) <= radi else None
        if dx is not None:
            schnittpunkte.append((cs[t1].x + dx, edge_y))
            schnittpunkte.append((cs[t1].x - dx, edge_y))
    
    schnittpunkte = [(sx, sy) for sx, sy in schnittpunkte if 0 <= sx <= domain_size and 0 <= sy <= domain_size]
    
    if len(schnittpunkte) < 2:
        return 0
    
    phi1 = math.atan2(schnittpunkte[0][1] - cs[t1].y, schnittpunkte[0][0] - cs[t1].x)
    phi2 = math.atan2(schnittpunkte[1][1] - cs[t1].y, schnittpunkte[1][0] - cs[t1].x)
    
    dqphi = abs(phi2 - phi1)
    if dqphi > np.pi:
        dqphi = 2 * np.pi - dqphi
    
    return dqphi         
            
    
def information_gathering():
    dist_tot = 0
    for i in range(0,N):
        if cs[i].mer:
            cs[i].listrad.append(cs[i].rad)
        else: cs[i].listrad.append(0)
        for g in range(0,N):
            
            dist_tot +=dist(i,g)/(N**2-N)
    dist_save.append(dist_tot)
    


def interpolation_aggr(t_save, dist_save, span): 
    xp = t_save
    yp = dist_save
    ep = np.interp(span, xp, yp)
    for i in range(0,len(ep)):
        aggrate[i] += ep[i]
        
        
    


'''initial conditions'''



for i in range(0,N):
    P = rn.uniform(20,upbound)
    rad = R_0
    x_zelle = rn.uniform(0,100)
    y_zelle = rn.uniform(0,100)
    dy = 0 
    dx = 0
    mer = True
    mer_w = []
    
   
    cell = ConvCells(P,rad, y_zelle, x_zelle, dy ,dx, mer, mer_w)
    cs.append(cell)
 
    cs[i].listx.append(cs[i].x)
    cs[i].listy.append(cs[i].y)
    



n = number_active_cells()
radial_prec(n)

for i in range(0,N):  
    cs[i].listrad.append(cs[i].rad)
    for g in range(0,N):
        dist_tot +=dist(i,g)/(N**2-N)
        
dist_save.append(dist_tot)
t_save.append(t_tot)



'''Simulation'''


while n > 1:
    
    wind_v = []
    merging_cells=[]
    
    
    mergers()
    
    calculating_displacement()
        
    dt = find_dt(wind_v,dl)
    
    displacement(dt,t_tot)
    
    n = number_active_cells()
    
    information_gathering()
    
    radial_prec(n)
    
interpolation_aggr(t_save, dist_save, span)


'''Plots'''


farben = plt.colormaps["hsv"]

#Aggregation rate

plt.figure(figsize=(10, 8))
plt.grid(True)
plt.title(f'mean distance of the cells throughout the simulation for $ r_0^2$={math.floor(r0**2)}')
plt.plot(span, aggrate, color='r')
plt.xlabel('time [min]')
plt.ylabel('distance [km]')
plt.show()








#trajectory plot

allemergingpointsx =[]
allemergingpointsy =[]

for s in range(N):
    for f in range(0,len(cs[s].mergepoint)):
        allemergingpointsx.append(cs[s].mergepoint[f][0])
        allemergingpointsy.append(cs[s].mergepoint[f][1])

for i in range(N):
    mergex = [] 
fig, ax = plt.subplots()
for i in range(N):
    ax.plot(cs[(N-i-1)].listx, cs[(N-i-1)].listy, 'o', color=farben((N-i-1)/N), markersize=3, label=f'Celle {N-i}')
    for xi, yi, si in zip(cs[(N-i-1)].listx, cs[(N-i-1)].listy, cs[N-i-1].listrad):
        circle = plt.Circle((xi, yi),si, color=farben((N-i-1)/N), alpha=0.05,)
        ax.add_patch(circle)
ax.plot(cs[0].mergepoint[-1][0],cs[0].mergepoint[-1][1], color="black", label="Last Point",marker='.')
ax.set_xlim(0,100)
ax.set_ylim(0,100)
ax.set_xlabel('x-coordinate [km]')
ax.set_ylabel('y-coordinate [km]')
ax.set_title('Trajectories of a single Simulation')
plt.grid(True)

plt.show()


#radial precipitation plot


for i in range(len(t_saveinc )):
    radial=np.linspace(1, len(rad_prec[0])+1)
    plt.bar(radial,rad_prec[i], edgecolor='black', alpha=0.7,color=colors )
    plt.xlabel('radius [km]')
    plt.ylabel('Precipitation [mm/h]')
    plt.title(f'radial Precipitation at time {round(t_save[i],2)}min')
    plt.grid(True)
    plt.show()





