"""
Reconstruccion COMPLETA de la geometria del exoesqueleto a partir de las
ecuaciones de CinematicaExoFinal.m / CinematicaExoModificada.m, para generar
un diagrama a escala con TODOS los eslabones, articulaciones y offsets.

Topologia (un solo motor):
  Motor -> engranes (reng) -> [5B#1 + 4B#1] -> falange proximal + soportes S1,S2
        -> [5B#2 + 4B#2] -> falange medial -> [4B#3 nuevo] -> falange distal

NOTA SOBRE EL MONTAJE DEL TERCER MECANISMO (DIP):
  El 4B#3 que acopla la falange distal a la cadena cinematica se monta del lado
  DORSAL del dedo (lado opuesto a la palma) por restricciones fisicas: las falanges
  apoyan el lado palmar contra los objetos manipulados, asi que los soportes y
  eslabones del exoesqueleto deben quedar por el dorso, fuera del area de prension.

  Por la simetria del 4 barras plano respecto a su bancada (la falange medial), un
  mecanismo "palmar" y su reflejo "dorsal" producen la MISMA cinemativa de salida
  (mismo theta_fd, misma flexion DIP). Por eso la cinematica numerica de
  CinematicaExoModificada.m / exo_18.py NO cambia: solo la posicion fisica de las
  barras intermedias (manivela, acoplador, balancin) se refleja al dorso para el
  CAD. Los parametros Lpc, Lpd, Lac, BETA1, BETA2 conservan los mismos valores.

  En este script se reflejan visualmente las posiciones CRK3 y ROK3 respecto a la
  recta IFP-IFD (la bancada del 4B#3) para mostrar el mecanismo donde se construira
  realmente.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------- Parametros nominales (dedo indice, mm/grados) -----------------
Bancada1=18.0; Bancada2=20.0
Link1=35.0; Link2=49.0; Link3=25.0; Link4=20.0; Link5=25.0
Link6=55.0; Link7=35.0; Link8=52.0; Link10=35.0
c2=46.01
TETHA1inicial=109.0; THETA14B=90.0
hsp=17.0; dsp=18.0
fp=49.0; fm=26.0; fd=24.0
THETAauxfm=51.39; THETAauxfd=38.78
reng=2.0
# Tercer mecanismo (acoplamiento DIP, montaje dorsal en CAD)
Lpc=8.0; Lpd=18.0; Lac=8.86; BETA1=40.0; BETA2=110.0

r1=Link4; r2=Link3; r3=Bancada1/2.0; r4=Link1; r5=Link2
a=Link4; b=Link5; c=np.sqrt(hsp**2+dsp**2); d=Bancada2
theta14B=np.deg2rad(THETA14B)
r1m2=fp-2*dsp; r2m2=Link7; r3m2=Link5/2.0; r4m2=Link3; r5m2=Link6
a2=Link7; b2=Link8; d2=np.sqrt(hsp**2+dsp**2)

def four_bar_open(a_, b_, c_, d_, th2_, th1_):
    k1=d_*np.cos(th1_)+a_*np.cos(th2_)
    k2=d_*np.sin(th1_)+a_*np.sin(th2_)
    k3=k1**2+k2**2+c_**2-b_**2
    A1=-k3-2*k1*c_; B1=4*k2*c_; C1=2*k1*c_-k3
    return 2*np.arctan((-B1-np.sqrt(B1**2-4*A1*C1))/(2*A1))

def reflect_across_line(P, A, u_hat):
    """Refleja el punto P respecto a la recta que pasa por A en direccion u_hat (unitario).
    Util para colocar el 4B#3 del lado dorsal del dedo manteniendo la cinematica."""
    v = P - A
    parallel = np.dot(v, u_hat) * u_hat
    perp = v - parallel
    return A + parallel - perp

def geom(THETA2):
    th2=np.deg2rad(THETA2)
    th1=np.deg2rad(THETA2/reng+TETHA1inicial)
    P={}
    # --- pivotes fijos ---
    A=np.array([-r3,0.0])          # pivote crank 5B#1 / crank 4B#1
    B=np.array([ r3,0.0])          # pivote crank r4 del 5B#1
    MCF=np.array([-r3,-d])         # articulacion metacarpofalangica (pivote rocker 4B#1)
    P['A']=A; P['B']=B; P['MCF']=MCF
    # --- 5B#1 ---
    M4=A+r1*np.array([np.cos(th1),np.sin(th1)])     # tip crank r1 (=Link4)
    J2=B+r4*np.array([np.cos(th2),np.sin(th2)])     # tip crank r4 (=Link1)
    e=(r1*np.sin(th1)-r4*np.sin(th2))/(r4*np.cos(th2)-r1*np.cos(th1)+2*r3)
    f=(2*(r1*r3*np.cos(th1)+r3*r4*np.cos(th2))-r1**2+r2**2+r4**2-r5**2)/(2*(r4*np.cos(th2)-r1*np.cos(th1)+2*r3))
    daux=e**2+1
    g=2*(e*f-e*r1*np.cos(th1)+e*r3-r1*np.sin(th1))
    h=f**2-2*f*(r1*np.cos(th1)-r3)-2*r1*r3*np.cos(th1)+r1**2+r3**2-r2**2
    pyP=(-g+np.sqrt(g**2-4*daux*h))/(2*daux); pxP=e*pyP+f
    Pp=np.array([pxP,pyP])
    P['M4']=M4; P['J2']=J2; P['P']=Pp
    # --- 4B#1 ---
    th4a=four_bar_open(a,b,c,d,th1,theta14B)
    TH4a=np.rad2deg(th4a)
    if TH4a<0: TH4a+=360
    th4a=np.deg2rad(TH4a)
    thfp=np.deg2rad(TH4a+np.rad2deg(np.arctan2(hsp,dsp)))
    IFP=MCF+fp*np.array([np.cos(thfp),np.sin(thfp)])
    S1=MCF+c*np.array([np.cos(th4a),np.sin(th4a)])
    thps2=thfp-np.arctan2(hsp,fp-dsp)
    rs2=np.sqrt(hsp**2+(fp-dsp)**2)
    S2=MCF+rs2*np.array([np.cos(thps2),np.sin(thps2)])
    P['IFP']=IFP; P['S1']=S1; P['S2']=S2; P['thfp']=thfp
    # --- 5B#2 -> P2 ---
    throll=np.arctan2(M4[1]-S1[1],M4[0]-S1[0])
    def angloc(pt_to,pt_from):
        ang=np.rad2deg(np.arctan2(pt_to[1]-pt_from[1],pt_to[0]-pt_from[0]))
        if ang<0: ang+=360
        return np.deg2rad(ang)-throll
    th1m2=angloc(S2,S1); th2m2=angloc(Pp,M4)
    em2=(r1m2*np.sin(th1m2)-r4m2*np.sin(th2m2))/(r4m2*np.cos(th2m2)-r1m2*np.cos(th1m2)+2*r3m2)
    fm2=(2*(r1m2*r3m2*np.cos(th1m2)+r3m2*r4m2*np.cos(th2m2))-r1m2**2+r2m2**2+r4m2**2-r5m2**2)/(2*(r4m2*np.cos(th2m2)-r1m2*np.cos(th1m2)+2*r3m2))
    dauxm2=em2**2+1
    gm2=2*(em2*fm2-em2*r1m2*np.cos(th1m2)+em2*r3m2-r1m2*np.sin(th1m2))
    hm2=fm2**2-2*fm2*(r1m2*np.cos(th1m2)-r3m2)-2*r1m2*r3m2*np.cos(th1m2)+r1m2**2+r3m2**2-r2m2**2
    pyP2=(-gm2+np.sqrt(gm2**2-4*dauxm2*hm2))/(2*dauxm2); pxP2=em2*pyP2+fm2
    p2=np.hypot(pxP2,pyP2); th2p2=np.arctan2(pyP2,pxP2)
    AUX=(S1+M4)/2.0
    P2=p2*np.array([np.cos(th2p2+throll),np.sin(th2p2+throll)])+AUX
    P['P2']=P2
    # --- 4B#2 -> medial ---
    th14B2=np.arctan2(S2[1]-IFP[1],S2[0]-IFP[0])
    a24=np.rad2deg(np.arctan2(P2[1]-S2[1],P2[0]-S2[0]))
    if a24<0: a24+=360
    th24B2=np.deg2rad(a24)
    th4am2=four_bar_open(a2,b2,c2,d2,th24B2,th14B2)
    TH4am2=np.rad2deg(th4am2)
    if TH4am2<0: TH4am2+=360
    P3=IFP+c2*np.array([np.cos(np.deg2rad(TH4am2)),np.sin(np.deg2rad(TH4am2))])
    thfm=np.deg2rad(TH4am2+THETAauxfm)
    IFD=IFP+fm*np.array([np.cos(thfm),np.sin(thfm)])
    P['P3']=P3; P['IFD']=IFD; P['thfm']=thfm
    # --- 4B#3 (DIP) -> distal ---
    # Cinematica calculada con la convencion original (palmar). theta_fd resultante
    # es identico al de un montaje dorsal por la simetria del 4 barras respecto a su
    # bancada (la falange medial). Despues reflejamos las posiciones CRK3 y ROK3 al
    # lado dorsal para el dibujo.
    thfp_=P['thfp']
    alpha1=thfp_+np.deg2rad(BETA1)-thfm
    Acrk_palmar=IFP+Lpc*np.array([np.cos(thfp_+np.deg2rad(BETA1)),np.sin(thfp_+np.deg2rad(BETA1))])
    Ax=Lpc*np.cos(alpha1); Ay=Lpc*np.sin(alpha1)
    Px=Ax-fm; Py=Ay; R=np.hypot(Px,Py)
    K=(Px**2+Py**2+Lpd**2-Lac**2)/(2*Lpd)
    phi=np.arctan2(Py,Px); alpha2=phi-np.arccos(K/R)
    thfd=thfm+(alpha2-np.deg2rad(BETA2))
    Brok_palmar=IFD+Lpd*np.array([np.cos(thfd+np.deg2rad(BETA2)),np.sin(thfd+np.deg2rad(BETA2))])
    # Reflexion al lado dorsal respecto a la recta IFP-IFD (la bancada del 4B#3)
    u_fm = (IFD - IFP) / fm
    Acrk = reflect_across_line(Acrk_palmar, IFP, u_fm)
    Brok = reflect_across_line(Brok_palmar, IFP, u_fm)
    TIP=IFD+fd*np.array([np.cos(thfd),np.sin(thfd)])
    P['CRK3']=Acrk; P['ROK3']=Brok; P['TIP']=TIP; P['thfd']=thfd
    P['CRK3_palmar']=Acrk_palmar; P['ROK3_palmar']=Brok_palmar  # disponible para debug
    return P

# ----------------- VERIFICACION de longitudes -----------------
P=geom(0.0)
def dist(a,b): return np.linalg.norm(P[a]-P[b])
checks=[
 ("|M4-A|=Link4",dist('M4','A'),Link4),
 ("|P-M4|=Link3",dist('P','M4'),Link3),
 ("|P-J2|=Link2",dist('P','J2'),Link2),
 ("|M4-S1|=Link5",dist('M4','S1'),Link5),
 ("|S1-MCF|=sqrt(hsp^2+dsp^2)",dist('S1','MCF'),np.sqrt(hsp**2+dsp**2)),
 ("|S2-S1|=fp-2dsp",dist('S2','S1'),fp-2*dsp),
 ("|P2-S2|=Link7",dist('P2','S2'),Link7),
 ("|P2-P|=Link6",dist('P2','P'),Link6),
 ("|P3-IFP|=c2",dist('P3','IFP'),c2),
 ("|P3-P2|=Link8",dist('P3','P2'),Link8),
 ("|IFD-IFP|=fm",dist('IFD','IFP'),fm),
 ("|TIP-IFD|=fd",dist('TIP','IFD'),fd),
 ("|CRK3-IFP|=Lpc (dorsal)",dist('CRK3','IFP'),Lpc),
 ("|ROK3-IFD|=Lpd (dorsal)",dist('ROK3','IFD'),Lpd),
 ("|ROK3-CRK3|=Lac (dorsal)",dist('ROK3','CRK3'),Lac),
]
print("=== VERIFICACION DE RECONSTRUCCION (mm) ===")
ok=True
for name,val,exp in checks:
    good=abs(val-exp)<1e-6
    ok=ok and good
    print(f"  [{'OK' if good else 'XX'}] {name:32s} calc={val:8.3f}  esp={exp:8.3f}")
print("TODAS CORRECTAS" if ok else "HAY DISCREPANCIAS")



# ============================ DIBUJO DEL DIAGRAMA ============================
P = geom(0.0)   # posicion de referencia (dedo extendido, THETA2 = 0)

def pt(name): return P[name]

def bar(ax, n1, n2, color, lw=2.5, label=None, ls='-', z=3,
        label_offset=(0,0)):
    p, q = pt(n1), pt(n2)
    ax.plot([p[0], q[0]], [p[1], q[1]], ls, color=color, lw=lw, zorder=z,
            solid_capstyle='round')
    if label:
        m = (p + q) / 2.0 + np.array(label_offset)
        ax.annotate(label, m, fontsize=8, color=color, ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.12', fc='white', ec=color, alpha=0.9),
                    zorder=z+1)

def joint(ax, name, txt=None, dx=2.2, dy=2.2, c='k', ms=7, ha='left', va='bottom'):
    p = pt(name)
    ax.plot(p[0], p[1], 'o', color=c, ms=ms, mfc='white', mew=1.8, zorder=6)
    if txt:
        ax.annotate(txt, (p[0]+dx, p[1]+dy), fontsize=8.5, fontweight='bold',
                    color=c, zorder=7, ha=ha, va=va)

# Figura mas ancha para acomodar la leyenda fuera del area de dibujo
fig, ax = plt.subplots(figsize=(17, 10))

# Colores por subsistema
C_FRAME='#444444'; C_5B1='#1f77b4'; C_4B1='#2ca02c'; C_5B2='#9467bd'
C_4B2='#ff7f0e'; C_4B3='#d62728'; C_PHAL='#111111'

# --- Bancada / marco fijo ---
bar(ax,'A','B',C_FRAME,lw=3,label='Bancada1 (2*r3=%.0f)'%Bancada1,z=2,label_offset=(0,2.5))
bar(ax,'A','MCF',C_FRAME,lw=3,label='Bancada2 (d=%.0f)'%Bancada2,z=2,label_offset=(-5,0))
ax.fill([pt('A')[0],pt('B')[0],pt('MCF')[0]],
        [pt('A')[1],pt('B')[1],pt('MCF')[1]],color=C_FRAME,alpha=0.06,zorder=1)

# --- Engranes (motor unico) ---
from matplotlib.patches import Circle
ax.add_patch(Circle(pt('B'),6.0,fill=False,ec=C_FRAME,lw=1.5,ls=':',zorder=2))
ax.add_patch(Circle(pt('A'),3.0,fill=False,ec=C_FRAME,lw=1.5,ls=':',zorder=2))
ax.annotate('engranes\n(reng=%.0f)\nMOTOR'%reng,(pt('B')[0]+8,pt('B')[1]+9),
            fontsize=8,color=C_FRAME,fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.15',fc='white',ec=C_FRAME,alpha=0.85))

# --- Primer mecanismo de 5 barras ---
bar(ax,'A','M4',C_5B1,label='Link4',z=4)
bar(ax,'M4','P',C_5B1,label='Link3',z=4)
bar(ax,'B','J2',C_5B1,label='Link1',z=4)
bar(ax,'J2','P',C_5B1,label='Link2',z=4)

# --- Primer mecanismo de 4 barras ---
bar(ax,'M4','S1',C_4B1,label='Link5',z=4)
bar(ax,'MCF','S1',C_4B1,label='c=v(hsp^2+dsp^2)',z=4,label_offset=(0,-1))

# --- Segundo mecanismo de 5 barras ---
bar(ax,'S1','S2',C_5B2,label='fp-2dsp',z=4,label_offset=(0,2))
bar(ax,'S2','P2',C_5B2,label='Link7',z=4)
bar(ax,'P','P2',C_5B2,label='Link6',z=4)

# --- Segundo mecanismo de 4 barras ---
bar(ax,'P2','P3',C_4B2,label='Link8',z=4)
bar(ax,'IFP','P3',C_4B2,label='c2=%.2f'%c2,z=4)
bar(ax,'IFP','S2',C_4B2,ls='--',lw=1.6,label='d2 (bancada)',z=3,label_offset=(0,-2))

# --- Tercer mecanismo de 4 barras (DIP, montaje DORSAL para CAD) ---
bar(ax,'IFP','CRK3',C_4B3,label='Lpc=%.0f (BETA1=%.0f)'%(Lpc,BETA1),z=5,label_offset=(0,1.2))
bar(ax,'CRK3','ROK3',C_4B3,label='Lac=%.2f'%Lac,z=5,label_offset=(0,1.2))
bar(ax,'IFD','ROK3',C_4B3,label='Lpd=%.0f (BETA2=%.0f)'%(Lpd,BETA2),z=5,label_offset=(0,1.2))

# --- Falanges (cuerpos rigidos) ---
bar(ax,'MCF','IFP',C_PHAL,lw=6,label='Falange proximal (fp=%.0f)'%fp,z=2,label_offset=(0,-3))
bar(ax,'IFP','IFD',C_PHAL,lw=6,label='Falange medial (fm=%.0f)'%fm,z=2,label_offset=(0,-3))
bar(ax,'IFD','TIP',C_PHAL,lw=6,label='Falange distal (fd=%.0f)'%fd,z=2,label_offset=(2,-3))

# --- Articulaciones etiquetadas (con offsets ajustados para evitar choques) ---
joint(ax,'A','A (-r3,0)',c=C_FRAME,dx=-12,dy=3.5,ha='right')
joint(ax,'B','B (r3,0)',c=C_FRAME,dx=3,dy=-1)
joint(ax,'MCF','MCF',c=C_PHAL,dx=2,dy=-5)
joint(ax,'M4','M4',c=C_5B1,dx=2,dy=2)
joint(ax,'P','P',c=C_5B1,dx=2,dy=2)
joint(ax,'J2','J2',c=C_5B1,dx=3,dy=-1)
joint(ax,'S1','S1',c=C_4B1,dx=-2,dy=4,ha='right')
joint(ax,'S2','S2',c=C_5B2,dx=2,dy=4)
joint(ax,'P2','P2',c=C_5B2,dx=2,dy=2)
joint(ax,'P3','P3',c=C_4B2,dx=2,dy=2)
joint(ax,'IFP','IFP (PIP)',c=C_PHAL,dx=2,dy=-5)
joint(ax,'IFD','IFD (DIP)',c=C_PHAL,dx=2,dy=-5)
joint(ax,'CRK3',None,c=C_4B3,ms=5)
joint(ax,'ROK3',None,c=C_4B3,ms=5)
joint(ax,'TIP','Punta',c=C_PHAL,dx=3,dy=-2)

# --- Leyenda FUERA del area de dibujo (a la derecha) ---
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],color=C_FRAME,lw=3,label='Bancada / marco fijo'),
     Line2D([0],[0],color=C_5B1,lw=2.5,label='Mecanismo 5 barras #1'),
     Line2D([0],[0],color=C_4B1,lw=2.5,label='Mecanismo 4 barras #1\n(falange proximal)'),
     Line2D([0],[0],color=C_5B2,lw=2.5,label='Mecanismo 5 barras #2'),
     Line2D([0],[0],color=C_4B2,lw=2.5,label='Mecanismo 4 barras #2\n(falange medial)'),
     Line2D([0],[0],color=C_4B3,lw=2.5,label='Mecanismo 4 barras #3\n(acoplamiento DIP, DORSAL)'),
     Line2D([0],[0],color=C_PHAL,lw=6,label='Falanges (lado palmar)')]
ax.legend(handles=leg,loc='center left',bbox_to_anchor=(1.01, 0.5),
          fontsize=9.5,framealpha=0.95,title='Subsistemas',title_fontsize=10)

# Anotacion explicativa del montaje dorsal del 4B#3
ax.annotate(
    'NOTA: el 4B#3 (rojo) se monta DORSAL por\n'
    'restriccion fisica (no puede invadir el lado\n'
    'palmar). Por simetria del 4 barras respecto a\n'
    'su bancada (falange medial), la cinematica de\n'
    'salida (theta_fd) es identica que en montaje\n'
    'palmar; los valores Lpc, Lpd, Lac, BETA1,\n'
    'BETA2 son los mismos.',
    xy=(0.01, 0.01), xycoords='axes fraction',
    fontsize=8, color=C_4B3,
    bbox=dict(boxstyle='round,pad=0.4', fc='white', ec=C_4B3, alpha=0.9),
    zorder=10)

ax.set_aspect('equal'); ax.grid(True,ls=':',alpha=0.5)
ax.set_title('Diagrama cinematico completo del exoesqueleto de dedo indice\n'
             '(posicion de referencia THETA2=0; cotas en mm) - conexiones para el CAD',
             fontsize=13,fontweight='bold')
ax.set_xlabel('X (mm) - DORSAL hacia arriba (+Y), PALMAR hacia abajo (-Y)')
ax.set_ylabel('Y (mm)')
plt.tight_layout()
plt.savefig('diagrama_mecanismo_completo.png',dpi=150,bbox_inches='tight')
plt.close()
print("\nGuardado diagrama_mecanismo_completo.png")
print("\n=== POSICIONES (mm) ===")
for k in ['A','B','MCF','M4','J2','P','S1','S2','P2','P3','IFP','IFD','CRK3','ROK3','TIP']:
    print(f"  {k:5s} = ({P[k][0]:8.2f}, {P[k][1]:8.2f})")
print(f"\n  thfp = {np.rad2deg(P['thfp']):7.2f} deg")
print(f"  thfm = {np.rad2deg(P['thfm']):7.2f} deg")
print(f"  thfd = {np.rad2deg(P['thfd']):7.2f} deg")
print(f"  DIP relativo (thfd - thfm) = {np.rad2deg(P['thfd']-P['thfm']):7.2f} deg")
