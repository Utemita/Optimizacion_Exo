"""
Genera trayectorias_referencia.csv con la posicion en el plano (X, Y, mm) de
TODOS los pivotes del exoesqueleto en cada paso del barrido del motor THETA2,
calculadas con las MISMAS ecuaciones de CinematicaExoModificada.m (re-implementadas
en numpy en este script).

Sirve para verificar el modelo CAD: al hacer el estudio de movimiento en el CAD
con omega2 = 10 RPM (o cualquier velocidad equivalente), las posiciones de
IFP, IFD y TIP en cada paso deben coincidir con las de este CSV.

Uso:
    python3 generar_trayectorias_referencia.py
Genera:
    trayectorias_referencia.csv     (133 filas x 33 columnas)
"""
import numpy as np
import csv
from diagrama_mecanismo import geom

# Mismo barrido del motor que CinematicaExoModificada.m: THETA2 = 0..132 deg en
# pasos de 1 deg (133 puntos).
THs = np.arange(0, 133, 1.0)

cols = ['THETA2_deg',
        'A_x','A_y','B_x','B_y','MCF_x','MCF_y',
        'M4_x','M4_y','J2_x','J2_y','P_x','P_y',
        'S1_x','S1_y','S2_x','S2_y',
        'P2_x','P2_y','P3_x','P3_y',
        'IFP_x','IFP_y','IFD_x','IFD_y','TIP_x','TIP_y',
        'CRK3_x','CRK3_y','ROK3_x','ROK3_y',
        'theta_fp_deg','theta_fm_deg','theta_fd_deg',
        'DIP_rel_deg']

with open('trayectorias_referencia.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(cols)
    for T in THs:
        P = geom(T)
        thfp = np.rad2deg(P['thfp'])
        thfm = np.rad2deg(P['thfm'])
        thfd = np.rad2deg(P['thfd'])
        dip_rel = thfd - thfm
        row = [T,
               *P['A'], *P['B'], *P['MCF'],
               *P['M4'], *P['J2'], *P['P'],
               *P['S1'], *P['S2'],
               *P['P2'], *P['P3'],
               *P['IFP'], *P['IFD'], *P['TIP'],
               *P['CRK3'], *P['ROK3'],
               thfp, thfm, thfd, dip_rel]
        w.writerow([f'{x:.4f}' if isinstance(x,(int,float,np.floating,np.integer)) else x for x in row])

# Resumen
P0 = geom(0.0)
PE = geom(132.0)
print('Generado trayectorias_referencia.csv')
print(f'  {len(THs)} pasos del motor (THETA2 = 0 a 132 deg, paso 1 deg)')
print()
print('=== POSICIONES CLAVE EN POSICION DE REFERENCIA (THETA2 = 0 deg) ===')
for k in ['A','B','MCF','M4','J2','P','S1','S2','P2','P3','IFP','IFD','CRK3','ROK3','TIP']:
    print(f'  {k:5s} = ({P0[k][0]:8.3f}, {P0[k][1]:8.3f}) mm')
print()
print(f'  theta_fp = {np.rad2deg(P0["thfp"]):7.3f} deg     theta_fm = {np.rad2deg(P0["thfm"]):7.3f} deg')
print(f'  theta_fd = {np.rad2deg(P0["thfd"]):7.3f} deg     DIP_rel  = {np.rad2deg(P0["thfd"]-P0["thfm"]):7.3f} deg')
print()
print('=== POSICIONES AL FINAL DEL CIERRE (THETA2 = 132 deg) ===')
for k in ['IFP','IFD','TIP']:
    print(f'  {k:5s} = ({PE[k][0]:8.3f}, {PE[k][1]:8.3f}) mm')
print(f'  theta_fp = {np.rad2deg(PE["thfp"]):7.3f} deg     theta_fm = {np.rad2deg(PE["thfm"]):7.3f} deg')
print(f'  theta_fd = {np.rad2deg(PE["thfd"]):7.3f} deg     DIP_rel  = {np.rad2deg(PE["thfd"]-PE["thfm"]):7.3f} deg')
